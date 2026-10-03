#!/usr/bin/env python3
"""Assign genome segments to curated TAHV metadata and group isolates.

Reads the curated metadata TSV plus one nextclade classification TSV per
segment reference (nextclade run against references/<seg>/, treeless), and
writes results/metadata.tsv with segment columns final:

  segment         NCBI label (M-RNA already normalized to M) when present;
                  otherwise the best-coverage segment, if it clears
                  --min-coverage; otherwise left empty (unassigned).
  segment_source  ncbi | nextclade-argmax | "" (unassigned)
  ncl_cov_s/m/l   nextclade coverage against each segment reference
                  (0 when nextclade excluded the sequence from that run)
  ncl_subs        totalSubstitutions on the assigned segment
  segment_label_conflict  true when the nextclade coverage argmax disagrees
                  with the NCBI segment label (flag only, never overridden)
  duplicate_isolate       true when an isolate has >= 2 records for the same
                  segment (e.g. the two Bardos-92 clones HM036208-213, or
                  XJ0625's GenBank + RefSeq trios) - its tips are not
                  independent isolates
  assembly        isolate grouping for tanglegram-style views: the pure
                  isolate name (`isolate` column) when >= 2 records of that
                  isolate carry >= 2 distinct non-empty segments; empty
                  otherwise.

Joins the nextclade TSVs on `strain`, which is the display id written into
the FASTA headers by ingest/bin/set_strain_display.py.

Labelled records keep their NCBI label even when nextclade disagrees - the
nextclade runs here are assignment-only and do not gate anything.
"""
import argparse
import csv
import sys

NCL_COV_COLS = {"s": "ncl_cov_s", "m": "ncl_cov_m", "l": "ncl_cov_l"}


def read_nextclade(path: str) -> dict:
    """seqName -> {'coverage': float, 'subs': int} for one segment run."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cov = row.get("coverage") or "0"
            subs = row.get("totalSubstitutions") or "0"
            out[row["seqName"]] = {
                "coverage": float(cov) if cov else 0.0,
                "subs": int(float(subs)) if subs else 0,
            }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--nextclade", nargs="+", required=True)
    ap.add_argument("--segments", nargs="+", required=True, help="e.g. s m l")
    ap.add_argument("--min-coverage", type=float, default=0.5)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if len(args.segments) != len(args.nextclade):
        sys.exit("--segments and --nextclade must have the same length")
    ncl = {seg: read_nextclade(p) for seg, p in zip(args.segments, args.nextclade)}

    with open(args.metadata, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        records = list(reader)

    for col in (NCL_COV_COLS[s] for s in args.segments):
        if col not in fieldnames:
            fieldnames.append(col)
    for col in ("segment_source", "ncl_subs", "segment_label_conflict", "duplicate_isolate", "assembly"):
        if col not in fieldnames:
            fieldnames.append(col)

    assigned = []  # (record, segment) with NCBI label or argmax filled in
    for rec in records:
        seq_id = rec["strain"]  # display id; matches nextclade seqName
        for seg in args.segments:
            rec[NCL_COV_COLS[seg]] = f"{ncl[seg].get(seq_id, {}).get('coverage', 0.0):.4f}"
        label = (rec.get("segment") or "").strip()
        if label:
            rec["segment_source"] = "ncbi"
            best = label
            # flag-only: does the coverage argmax disagree with the NCBI label?
            # (only meaningful where nextclade aligned the record at all)
            scores = {s: ncl[s].get(seq_id, {}).get("coverage", 0.0) for s in args.segments}
            argmax_seg = max(scores, key=scores.get).upper()
            best_cov = max(scores.values())
            rec["segment_label_conflict"] = "true" if best_cov > 0 and argmax_seg != best else ""
        else:
            scores = {s: ncl[s].get(seq_id, {}).get("coverage", 0.0) for s in args.segments}
            best_seg = max(scores, key=scores.get)
            if scores[best_seg] >= args.min_coverage:
                best = best_seg.upper()
                rec["segment_source"] = "nextclade-argmax"
            else:
                best = ""
                rec["segment_source"] = ""
        rec["segment"] = best
        rec["ncl_subs"] = str(ncl[best.lower()].get(seq_id, {}).get("subs", "")) if best else ""
        assigned.append(rec)

    # duplicate_isolate: an isolate with >=2 records for the same segment has
    # more than one genome sequenced (clones, GenBank+RefSeq) - its tips are
    # not independent isolates. Flag every record of such an isolate.
    per_isolate_segment = {}
    for rec in assigned:
        isolate = (rec.get("isolate") or "").strip()
        if isolate and rec["segment"]:
            per_isolate_segment.setdefault(isolate, {}).setdefault(rec["segment"], 0)
            per_isolate_segment[isolate][rec["segment"]] += 1
    for rec in assigned:
        isolate = (rec.get("isolate") or "").strip()
        seg_counts = per_isolate_segment.get(isolate, {})
        rec["duplicate_isolate"] = "true" if any(n >= 2 for n in seg_counts.values()) else ""

    # assembly: pure isolate name, >=2 records, >=2 distinct non-empty segments
    by_isolate = {}
    for rec in assigned:
        isolate = (rec.get("isolate") or "").strip()
        if isolate:
            by_isolate.setdefault(isolate, []).append(rec)
    for isolate, group in by_isolate.items():
        segments = {rec["segment"] for rec in group if rec["segment"]}
        if len(group) >= 2 and len(segments) >= 2:
            for rec in group:
                rec["assembly"] = isolate

    # Hard guard: tip names must be unique WITHIN each segment build (auspice
    # silently random-renames duplicates at render time, which breaks tangle
    # pairing). Aliased names repeat across segments by design; that is fine.
    per_segment_names = {}
    for rec in assigned:
        if rec["segment"]:
            per_segment_names.setdefault(rec["segment"], {}).setdefault(rec["strain"], []).append(
                rec["accession"]
            )
    for segment, names in sorted(per_segment_names.items()):
        dupes = {name: accs for name, accs in names.items() if len(accs) > 1}
        if dupes:
            for name, accs in sorted(dupes.items()):
                print(f"DUPLICATE tip name in segment {segment}: {name} <- {accs}", file=sys.stderr)
            sys.exit(f"duplicate strain names within segment {segment}; fix defaults/tip_name_aliases.tsv")

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(assigned)

    counts = {}
    for rec in assigned:
        key = (rec["segment_source"] or "unassigned", rec["segment"] or "-")
        counts[key] = counts.get(key, 0) + 1
    print("assign_segment: segment assignment counts:")
    for (source, seg), n in sorted(counts.items()):
        print(f"  {source or 'unassigned':18s} {seg or '-'}: {n}")
    n_asm = sum(1 for rec in assigned if rec.get("assembly"))
    print(f"assign_segment: {n_asm} records grouped into assemblies")


if __name__ == "__main__":
    main()
