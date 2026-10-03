#!/usr/bin/env python3
"""Build the record id ("strain + accession") in the curate stream.

Coordinator decision (2026-10-02): displayed tip names are isolate name +
accession; the pure isolate name must survive unmangled.

Tool reality: the strain value is the record id everywhere downstream (the
combined FASTA header, nextclade seqName, newick tip names, auspice nodes),
and the combined FASTA spans all segments - so ids must be globally unique.
Newick/FASTA ids are single tokens (whitespace/quotes/parens are newick syntax
or split FASTA ids), so the display form is one token:

    strain   (id)         XJ0625_EU622819   Prototype_92_Bardos_HM036208
    isolate  (pure name)  XJ0625            Prototype '92' Bardos

The tanglegram tip names (segment-independent biological IDs) are derived from
these in the PHYLO workfow (phylogenetic/bin/tangle_names.py + defaults/
tip_name_aliases.tsv), where per-segment naming makes byte-identical names
across segment trees safe. Do NOT alias ids here: two records of one isolate
would produce duplicate headers in the combined FASTA and mis-join every
downstream lookup.

Rules:
  - isolate := strain as it entered this step (pure lineage name; after
    transform-strain-name's accession backup, strain is never empty)
  - strain  := if strain is exactly the accession (with or without version),
    the accession.version; else <isolate>_<accession-without-version>, with
    characters unsafe for ids collapsed to "_". Collisions (same isolate +
    accession base, different version) get the full accession.version.

Usage: cat stream.ndjson | set_strain_display.py \
    --init-empty-field date_imputed --init-empty-field date_source > out.ndjson
"""
import argparse
import json
import re
import sys
from collections import Counter

UNSAFE = re.compile(r"[\s'\"()[\]{}:;,|#>]+")
VERSION_SUFFIX = re.compile(r"\.\d+$")


def token(text: str) -> str:
    return UNSAFE.sub("_", text).strip("_")


def build(rec: dict) -> str:
    accession = rec.get("accession") or ""
    isolate = (rec.get("strain") or "").strip()
    if isolate == accession or isolate == VERSION_SUFFIX.sub("", accession):
        # strain was backfilled from the accession (no isolate lineage):
        # the accession alone is already the isolate name.
        return accession
    return f"{token(isolate)}_{VERSION_SUFFIX.sub('', accession)}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--init-empty-field",
        action="append",
        default=[],
        help="create this field, empty, on every record (apply-record-annotations "
        "only overwrites existing fields, so annotation-set fields must pre-exist)",
    )
    args = ap.parse_args()
    records = [json.loads(line) for line in sys.stdin if line.strip()]
    for rec in records:
        for field in args.init_empty_field:
            rec[field] = ""
        # discrete collection year for auspice categorical coloring; "" when
        # undated at this point (imputed dates are applied later in the pipe)
        rec["year"] = (rec.get("date") or "")[:4].isdigit() and (rec.get("date") or "")[:4] or ""

    displays = Counter(build(rec) for rec in records)
    seen: dict = {}
    for rec in records:
        isolate = (rec.get("strain") or "").strip()
        display = build(rec)
        if displays[display] > 1:
            # collision on the unversioned accession: use the full versioned one
            display = f"{token(isolate)}_{rec['accession']}"
        if display in seen:
            sys.exit(f"strain display name collision: {display}")
        seen[display] = rec
        rec["isolate"] = isolate  # pure name, unmangled
        rec["strain"] = display   # display id used by every downstream stage

    for rec in records:
        print(json.dumps(rec, ensure_ascii=False))
    print(f"set_strain_display: {len(records)} records, {len(seen)} unique strain ids", file=sys.stderr)


if __name__ == "__main__":
    main()
