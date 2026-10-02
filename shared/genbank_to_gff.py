#!/usr/bin/env python3
"""GenBank -> GFF3 for Nextclade (--input-annotation) AND augur translate (--reference-sequence).

Reads GenBank file(s), writes one GFF3 per record (per segment): both tools
want a single seqid per file.

Shape produced (satisfies both tools, verified 2026-10-02, nextclade 3.21.2 / augur 34.0.0):

    ##gff-version 3
    ##sequence-region NC_055206.1 1 978                      # augur needs this ('nuc' annotation)
    NC_055206.1 RefSeq gene 79 786 . + . ID=gene:N;Name=N;gene=N;locus_tag=KM532_sSgp1
    NC_055206.1 RefSeq CDS  79 786 . + 0 ID=CDS:N;Parent=gene:N;Name=N;gene=N;locus_tag=KM532_sSgp1

Why each part is there:
  nextclade: fundamental unit is the CDS; the CDS name comes from ID; same ID across
    rows = one joined (multi-exon) CDS. The gene row is optional (display only).
  augur: ignores CDS rows for gene discovery, reads type=gene rows and takes the gene
    name from attributes gene= / gene_name= / locus_tag=; ##sequence-region (or a
    region/source row) is mandatory or augur errors with "didn't define any
    information ... 'nuc' annotation".

Usage: genbank_to_gff.py [-o OUTDIR] FILE1.gb [FILE2.gb ...] [--combined OUT.gff]
"""
import argparse
from pathlib import Path
from Bio import SeqIO

NAME_KEYS = ("gene", "locus_tag", "protein_id", "product")  # first hit wins

def feature_name(ft, counter):
    for k in NAME_KEYS:
        if k in ft.qualifiers:
            return ft.qualifiers[k][0].replace(" ", "_")
    counter[0] += 1
    return f"CDS{counter[0]}"

def convert_record(rec, counter):
    rows = [f"##sequence-region {rec.id} 1 {len(rec.seq)}"]
    for ft in rec.features:
        if ft.type not in ("gene", "CDS"):
            continue
        name = feature_name(ft, counter)
        strand = "+" if ft.location.strand == 1 else "-"
        prefix = "gene" if ft.type == "gene" else "CDS"
        attr = f"ID={prefix}:{name};Name={name};gene={name}"
        if "locus_tag" in ft.qualifiers:
            attr += f";locus_tag={ft.qualifiers['locus_tag'][0]}"
        if ft.type == "gene":
            rows.append(f"{rec.id}\t{rec.id}\tgene\t{int(ft.location.start)+1}\t{int(ft.location.end)}\t.\t{strand}\t.\t{attr}")
            continue
        # CDS: one row per exon/part, same ID => joined into one CDS
        cum = 0
        for part in ft.location.parts:
            start, end = int(part.start)+1, int(part.end)  # GFF is 1-based inclusive
            phase = (3 - (cum % 3)) % 3
            rows.append(f"{rec.id}\t{rec.id}\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attr}")
            cum += end - start + 1
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genbank", nargs="+")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--combined", help="also write all records into one GFF (nextclade only; NOT augur)")
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    counter = [0]
    combined = []
    for gb in args.genbank:
        for rec in SeqIO.parse(gb, "genbank"):
            rows = convert_record(rec, counter)
            out = outdir / f"{rec.id}.gff"
            out.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
            print(f"{out} ({len(rows)-1} features)")
            combined.extend(rows)
    if args.combined:
        Path(args.combined).write_text("##gff-version 3\n" + "\n".join(combined) + "\n")

if __name__ == "__main__":
    main()
