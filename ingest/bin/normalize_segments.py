#!/usr/bin/env python3
"""Normalize NCBI Datasets segment labels in a `augur curate` NDJSON stream.

NCBI labels 19 TAHV M-segment records "M-RNA" (GQ386823-841); left alone they
would not group with M records downstream. Mechanical label fix, not curation -
the map is a config knob (ingest/defaults/config.yaml `segment_map`). Labels
not present in the map (S/M/L, empty) pass through unchanged.

Usage: cat records.ndjson | normalize_segments.py --map "M-RNA=M" > out.ndjson
"""
import argparse
import json
import sys


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--map",
        default="",
        help='comma-separated OLD=NEW pairs, e.g. "M-RNA=M"',
    )
    args = ap.parse_args()
    mapping = dict(pair.split("=", 1) for pair in args.map.split(",") if pair)
    n_changed = 0
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        seg = rec.get("segment") or ""
        if seg in mapping:
            rec["segment"] = mapping[seg]
            n_changed += 1
        print(json.dumps(rec, ensure_ascii=False))
    print(f"normalize_segments: rewrote {n_changed} segment labels", file=sys.stderr)


if __name__ == "__main__":
    main()
