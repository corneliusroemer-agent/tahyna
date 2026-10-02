#!/usr/bin/env python3
"""Build the display tip name ("strain + accession") in the curate stream.

Coordinator decision (2026-10-02): displayed tip names are isolate name +
accession, e.g. XJ0625 (OP727994); the pure isolate name must survive
unmangled for tanglegram readability.

Tool reality: the strain value is the record id everywhere downstream (FASTA
headers, nextclade seqName, newick tip names, auspice nodes). Newick/FASTA
ids are single tokens - whitespace, quotes and parentheses either split the id
or are newick syntax. So the display form is the same name in one token:

    strain   (id/display)  XJ0625_EU622819   Prototype_92_Bardos_HM036208
    isolate  (pure name)   XJ0625            Prototype '92' Bardos

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

UNSAFE = re.compile(r"[\s'\"()[\]{}:;,|#>]+")
VERSION_SUFFIX = re.compile(r"\.\d+$")


def token(text: str) -> str:
    return UNSAFE.sub("_", text).strip("_")


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

    from collections import Counter

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


def build(rec: dict) -> str:
    accession = rec.get("accession") or ""
    isolate = (rec.get("strain") or "").strip()
    if isolate == accession or isolate == VERSION_SUFFIX.sub("", accession):
        # strain was backfilled from the accession (no isolate lineage):
        # the accession alone is already the isolate name.
        return accession
    return f"{token(isolate)}_{VERSION_SUFFIX.sub('', accession)}"


if __name__ == "__main__":
    main()
