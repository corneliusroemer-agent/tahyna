#!/usr/bin/env python3
"""Bootstrap references/{l,m,s}/ from fixture/refseq_segments.{gb,fna}.

Per-segment reference package for nextclade + augur:

    references/<seg>/reference.fasta   # RefSeq segment (nextclade -r, augur)
    references/<seg>/annotation.gff3   # both tools; see shared/genbank_to_gff.py
    references/<seg>/reference.gb      # GenBank copy (human reference)
    references/<seg>/pathogen.json     # NOT written here: curated thresholds,
                                       # checked in by hand, never overwritten

RefSeq trio: NC_055207.1 (L, 6979 nt), NC_055205.1 (M, 4489 nt),
NC_055206.1 (S, 978 nt, ambisense: N 79..786 + NSs 98..391, overlapping).

Run once from the repo root:
    python3 shared/build_references.py
Outputs are committed; rerunning overwrites fasta/gff/gb but not pathogen.json.
"""
import sys
from pathlib import Path

from Bio import SeqIO

REPO = Path(__file__).resolve().parent.parent
ACCESSION_TO_SEG = {
    "NC_055207.1": "l",  # RNA-dependent RNA polymerase, 6979 nt
    "NC_055205.1": "m",  # Gn/Gc polyprotein, 4489 nt
    "NC_055206.1": "s",  # ambisense N + NSs, 978 nt
}


def main() -> None:
    import subprocess

    gff_script = REPO / "shared" / "genbank_to_gff.py"
    # 1. GFF3 per record (single-seqid files; both tools require one seqid per file)

    res = subprocess.run(
        [sys.executable, str(gff_script), "-o", str(REPO / "references"),
         str(REPO / "fixture" / "refseq_segments.gb")],
        check=False,
    )
    if res.returncode != 0:
        sys.exit("genbank_to_gff.py failed")
    # outputs land as references/NC_05520X.1.gff

    # 2. Per-segment directories: fasta + gb + move the gff in
    for rec in SeqIO.parse(REPO / "fixture" / "refseq_segments.fna", "fasta"):
        seg = ACCESSION_TO_SEG[rec.id]
        d = REPO / "references" / seg
        d.mkdir(parents=True, exist_ok=True)
        SeqIO.write(rec, d / "reference.fasta", "fasta")
    for rec in SeqIO.parse(REPO / "fixture" / "refseq_segments.gb", "genbank"):
        seg = ACCESSION_TO_SEG[rec.id]
        d = REPO / "references" / seg
        SeqIO.write(rec, d / "reference.gb", "genbank")
    for acc, seg in ACCESSION_TO_SEG.items():
        src = REPO / "references" / f"{acc}.gff"
        dst = REPO / "references" / seg / "annotation.gff3"
        src.replace(dst)
    print("references/{l,m,s}/ built (pathogen.json untouched)")


if __name__ == "__main__":
    main()
