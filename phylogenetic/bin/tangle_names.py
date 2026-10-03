#!/usr/bin/env python3
"""Rename one segment's tips to segment-independent tangle names.

Auspice tanglegrams connect tips whose top-level `name` is byte-identical
across segment trees (no normalisation, no attribute lookup - verified in
investigations/2026-10-02-tahyna-nextstrain/06-auspice-tangle-condition.md).
Ingest ids are globally unique (<isolate>_<accession>), which never match
across segments; this script renames a segment's QC-passed metadata + aligned
FASTA to the biological tip ID:

    curated alias (defaults/tip_name_aliases.tsv, by accession.version or
    accession) -> else the pure isolate name, tokenized -> else the accession.

The alias table pairs multi-record isolates deliberately (XJ0625 GenBank vs
RefSeq trio, the two Bardos-92 clones, the 181-57 spelling mismatch); pairing
S-tip to M-tip is curator knowledge no local rule can derive. Names must be
unique WITHIN this segment - violations are a hard error naming the fix
(auspice would silently random-rename duplicates at render time instead).
"""
import argparse
import csv
import re
import sys

UNSAFE = re.compile(r"[\s'\"()[\]{}:;,|#>]+")


def token(text: str) -> str:
    return UNSAFE.sub("_", text).strip("_")


def load_aliases(path: str) -> dict:
    aliases = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].lstrip().startswith("#"):
                continue
            if len(row) != 2 or not row[1].strip():
                sys.exit(f"alias table malformed line: {row!r}")
            aliases[row[0].strip()] = row[1].strip()
    return aliases


def tangle_name(rec: dict, aliases: dict) -> str:
    accession = rec.get("accession") or ""
    if accession in aliases:
        return aliases[accession]
    isolate = (rec.get("isolate") or rec.get("strain") or "").strip()
    if not isolate or isolate == accession:
        return accession
    return token(isolate)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True, help="QC-passed metadata TSV")
    ap.add_argument("--fasta", required=True, help="aligned FASTA (accession-based ids)")
    ap.add_argument("--aliases", required=True)
    ap.add_argument("--output-metadata", required=True)
    ap.add_argument("--output-fasta", required=True)
    args = ap.parse_args()
    aliases = load_aliases(args.aliases)

    with open(args.metadata, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        records = list(reader)

    renames = {}
    names = {}
    for rec in records:
        new = tangle_name(rec, aliases)
        old = rec["strain"]
        renames[old] = new
        names.setdefault(new, []).append(old)
        rec["strain"] = new
    dupes = {n: olds for n, olds in names.items() if len(olds) > 1}
    if dupes:
        for n, olds in sorted(dupes.items()):
            print(f"DUPLICATE tangle name in this segment: {n} <- {olds}", file=sys.stderr)
        sys.exit(
            "duplicate tangle names; add disambiguating rows to "
            "phylogenetic/defaults/tip_name_aliases.tsv"
        )

    with open(args.output_metadata, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    with open(args.fasta) as src, open(args.output_fasta, "w") as dst:
        for line in src:
            if line.startswith(">"):
                old = line[1:].strip().split()[0]
                # headers not in the metadata (the OUTGROUP tip) pass through
                dst.write(f">{renames.get(old, old)}\n")
            else:
                dst.write(line)

    print(f"tangle_names: {len(renames)} tips renamed ({len(aliases)} alias rows available)")


if __name__ == "__main__":
    main()
