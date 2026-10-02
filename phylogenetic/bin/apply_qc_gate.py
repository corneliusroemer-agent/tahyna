#!/usr/bin/env python3
"""Apply the approved nextclade QC gate to one segment's filtered metadata.

Approved criteria (docs/qc-review.md, user sign-off 2026-10-02): drop on ANY of
  - aligned coverage < --min-coverage
  - missing data > --missing-data-max
  - totalSubstitutions > --max-subs  (divergence band, fraction of reference)
  - frameshifts > 0
  - premature stops > 0
except records listed in --exceptions (kept, flagged qc_exception=true, their
failing metric stays visible in the report).

Writes the kept metadata (+ qc_exception column) and a dropped-strains report
with the exact failing metric(s), so nothing disappears silently.
"""
import argparse
import csv
import sys

COLS = {
    "coverage": "coverage",
    "missing": "qc.missingData.totalMissing",
    "subs": "totalSubstitutions",
    "frameshifts": "qc.frameShifts.totalFrameShifts",
    "stops": "qc.stopCodons.totalStopCodons",
}


def num(value):
    return float(value) if value not in (None, "") else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True, help="length-filtered metadata TSV")
    ap.add_argument("--nextclade", required=True, help="nextclade TSV for this segment")
    ap.add_argument("--min-coverage", type=float, required=True)
    ap.add_argument("--missing-data-max", type=int, required=True)
    ap.add_argument("--max-subs", type=int, required=True)
    ap.add_argument("--zero-frameshifts", action="store_true")
    ap.add_argument("--zero-stop-codons", action="store_true")
    ap.add_argument("--exceptions", nargs="*", default=[])
    ap.add_argument("--output", required=True, help="kept metadata TSV")
    ap.add_argument("--report", required=True, help="dropped-strains TSV")
    args = ap.parse_args()

    ncl = {}
    with open(args.nextclade, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ncl[row["seqName"]] = row

    with open(args.metadata, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        records = list(reader)
    for col in ("qc_exception", "coverage"):
        if col not in fieldnames:
            fieldnames.append(col)

    kept, dropped = [], []
    for rec in records:
        row = ncl.get(rec["strain"], {})
        cov = num(row.get(COLS["coverage"]))
        rec["coverage"] = f"{cov:.4f}"  # continuous coloring: aligned fraction of the reference
        missing = int(num(row.get(COLS["missing"])))
        subs = int(num(row.get(COLS["subs"])))
        frameshifts = int(num(row.get(COLS["frameshifts"])))
        stops = int(num(row.get(COLS["stops"])))
        reasons = []
        if cov < args.min_coverage:
            reasons.append(f"coverage<{args.min_coverage} (cov={cov:.3f})")
        if missing > args.missing_data_max:
            reasons.append(f"missing>{args.missing_data_max} (missing={missing})")
        if subs > args.max_subs:
            reasons.append(f"subs>{args.max_subs} (subs={subs})")
        if args.zero_frameshifts and frameshifts > 0:
            reasons.append(f"frameshifts={frameshifts}")
        if args.zero_stop_codons and stops > 0:
            reasons.append(f"stops={stops}")
        if reasons and rec["strain"] in args.exceptions:
            rec["qc_exception"] = "true"
            print(
                f"QC EXCEPTION kept: {rec['strain']} fails [{'; '.join(reasons)}] "
                "- kept via exceptions list, qc_exception=true",
                file=sys.stderr,
            )
            kept.append(rec)
        elif reasons:
            dropped.append((rec, reasons))
        else:
            rec["qc_exception"] = ""
            kept.append(rec)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept)
    with open(args.report, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["strain", "accession", "failing_criteria"])
        for rec, reasons in dropped:
            writer.writerow([rec["strain"], rec["accession"], "; ".join(reasons)])

    print(
        f"QC gate: {len(kept)} kept ({len([r for r in kept if r['qc_exception'] == 'true'])} by exception), "
        f"{len(dropped)} dropped"
    )
    for rec, reasons in dropped:
        print(f"  dropped: {rec['strain']} ({rec['accession']}): {'; '.join(reasons)}")


if __name__ == "__main__":
    main()
