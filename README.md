# tahyna — a minimal Nextstrain workflow for Tahyna virus

Two self-contained Snakemake workflows (NCBI ingest + per-segment phylogenetic
builds) for **Tahyna virus** (TAHV, *Orthobunyavirus tahynaense*, California
serogroup). Built to the rotavirus/lassa two-directory archetype; no vendored
tooling, no CI, no upload machinery — a local-first repo.

- **Ingest**: NCBI Datasets (species taxid **3052445**, chosen over the virus
  node 45270 because 6 records are classified only at the species level) →
  `augur curate` chain → per-segment nextclade coverage assignment → one
  metadata TSV + sequences FASTA (all segments).
- **Phylogenetic**: per segment S / M / L — length filter → nextclade align +
  QC (custom reference/GFF/pathogen.json) → approved QC gate → iqtree →
  refine → ancestral → translate → traits → export v2.

## Run

```sh
# 1. ingest (uses the committed 121-record fixture by default)
cd ingest
snakemake -s Snakefile -j4 --dry-run
time snakemake -s Snakefile -j4

# 2. phylogenetic builds
cd ../phylogenetic
snakemake -s Snakefile -j4 --dry-run
time snakemake -s Snakefile -j4
```

Run each workflow from its own directory (`.snakemake/` lands there). Outputs:
`ingest/results/{metadata.tsv,sequences.fasta}`, `auspice/tahyna_{s,m,l}.json`
(repo root — the canonical location nextstrain.org/community serves from).

**Live refetch** (instead of the fixture): `cd ingest && snakemake -s Snakefile -j4 --config live_fetch=true`
— runs `datasets download virus genome taxon 3052445` against NCBI; everything
downstream is identical. The fixture is `fixture/ncbi_dataset.zip`
(121 records, captured 2026-10-02).

### View on nextstrain.org

Once this repo is **public** (keep Actions off), nextstrain.org/community
serves `auspice/<name>.json` from the default branch directly:

- S: https://nextstrain.org/community/corneliusroemer-agent/tahyna/tahyna_s
- M: https://nextstrain.org/community/corneliusroemer-agent/tahyna/tahyna_m
- L: https://nextstrain.org/community/corneliusroemer-agent/tahyna/tahyna_l
- Tangle, S vs M (shared isolate prefixes make the pairing readable):
  https://nextstrain.org/community/corneliusroemer-agent/tahyna/tahyna_s:community/corneliusroemer-agent/tahyna/tahyna_m

### View the trees locally

```sh
auspice view --datasetDir auspice     # from the repo root
# open http://localhost:4000 — datasets tahyna_s, tahyna_m, tahyna_l
```

`auspice` 2.73 is in the `bioinfo` env. The JSONs are self-contained (v2
schema, root sequences inline), so any static file server works too.

## Layout

```
ingest/        Snakefile + defaults/ (config.yaml, date_imputation.tsv) + bin/
phylogenetic/  Snakefile + defaults/ (config.yaml, auspice_config.json) + bin/
auspice/       tahyna_{s,m,l}.json — the served deliverable (nextstrain.org/community)
references/    l/ m/ s/: reference.fasta, annotation.gff3, reference.gb, pathogen.json
shared/        genbank_to_gff.py, build_references.py (bootstrap, already run)
fixture/       committed NCBI Datasets package + RefSeq GenBank/FASTA
docs/          qc-review.md, skill-draft/, this README
```

Regenerate `references/` with `python3 shared/build_references.py` (never
touches the curated `pathogen.json` files).

## Config knobs

Ingest (`ingest/defaults/config.yaml`):

| knob | default | note |
|---|---|---|
| `ncbi_taxon_id` | `"3052445"` | species node; `45270` loses 6 records |
| `live_fetch` / `cached_package` | `false` / fixture | live-vs-cached switch |
| `ncbi_datasets_fields` | 18 fields | lean by design; Title Case drives the rename map |
| `segment_map` | `M-RNA: M` | 19 records NCBI labels `M-RNA` |
| `nextclade.assign_min_coverage` | `0.5` | blank-label records only: argmax segment kept if coverage ≥ floor |
| `curate.annotations` | `defaults/date_imputation.tsv` | curated 27-record date imputation (see "What to review") |

Phylogenetic (`phylogenetic/defaults/config.yaml`):

| knob | default | note |
|---|---|---|
| `segments.<seg>.min_length` | S 490 / M 2245 / L 3490 | ~50% of each RefSeq segment |
| `qc_gate.*` | coverage ≥ 0.5; missing < 500/2250/3500; subs ≤ 25% of ref; 0 frameshifts; 0 stops | approved criteria, `docs/qc-review.md` |
| `qc_gate.exceptions` | `PJ01_OP727996` | kept despite 1 frameshift (`qc_exception=true`) |
| `refine.*` | root best, coalescent opt, iqd 10, no stochastic-resolve | unpinned; polytomies kept |
| `refine.overrides.s.timetree` | `false` | S is a divergence tree — see below |
| `ancestral.inference` | `joint` | |
| `traits.columns` | country, region | |

Tool versions (validated; already in the devcontainer `bioinfo` env — do not
create a new env; `environment.yaml` documents the pins): augur 34.0.0,
nextclade 3.21.2, snakemake 9.21.0, datasets/dataformat 18.33.1, iqtree 3.1.2
(augur calls `iqtree3`), treetime 0.12.1, auspice 2.73.0, biopython for the
bootstrap. `augur tree --method iqtree` passes augur's defaults verbatim:
**fixed GTR** + `--ninit 2 -n 2 --epsilon 0.05 -T AUTO --redo` (verified in
augur source; matches nextstrain/ncov practice, which overrides no model).

## Validated run (the committed fixture)

Ingest: 121 records in → 121 curated, unique ids, FASTA↔metadata exact match.
Segments 54 S + 31 M + 25 L by NCBI label, +3 by coverage argmax, 8 unassigned
(≤294 nt fragments aligning to nothing). 0 undated (27 imputed, flagged).
15 lab-host kept, 12 duplicate-isolate flagged, 0 segment-label conflicts.
≈4 s wall.

Phylo: S 39 / M 33 / L 24 tips after floor + QC gate (PJ01 L kept by
exception). Wall ~56 s total, dominated by iqtree (~10 s/segment) and refine
(~7 s/segment); nextclade/translate/ancestral < 1 s each.

Timetrees: M rate **5.8e-5** subs/site/yr (σ 2.6e-5), L **5.3e-5** (σ 1.8e-5)
— both in the plausible 1e-5–1e-4 band. Root ages are poorly constrained (M
95% CI spans ~3900 years) — small sample, reported as-is. **S has no timetree:
** treetime's rate estimate is negative on every rooting (deep Asia/Europe
split; the European lineage is temporally flat 1957→2021, matching the
literature's "no clear temporal clustering" for TAHV S), so S ships as a
divergence tree via `refine.overrides.s` — evidence in
`docs/qc-review.md` §6.

Literature checks (numbers in `investigations/2026-10-02-tahyna-nextstrain/05-build.md`):
Bardos-92 prototype nests inside shallow European diversity on all three
segments; the Asian lineage is deeply split (3–30x within-Europe distances);
XJ0708's M is anomalously distant (0.232 to nearest neighbour, vs 0.012–0.020
for its S/L) — the known reassortment signal; PJ01 is deeply divergent in M
and L with S comparatively close to the Inner Mongolia strains.

## What to review (self-serve)

1. **Date imputation** — `ingest/defaults/date_imputation.tsv`: 27 records,
   3-column format (accession, field, value), method + sources documented in
   the header. Per-source breakdown: Camp et al. 2021 Suppl. Table S5 isolate
   years (19 GQ3868xx records, 1958–1968), strain-name years verified in the
   literature (181/57 → 1957, Růžek 2020; prototype "92" Bardos → 1958),
   prototype-lineage inference (EU185046.1, lower confidence), release-year
   fallback (AX230490.1 only). Imputed records carry `date_imputed=true`,
   `date_source=<method tag>`.
2. **QC gate + S timetree** — `docs/qc-review.md`: approved criteria, the
   per-segment outcome distribution, the PJ01 exception, and the evidence
   behind the S divergence-tree override.

## Caveats

- **Bardos-92 clones**: HM036208-213 are two biologically cloned genomes of
  the same 1958 isolate (Camp 2021 S5), and XJ0625 contributes both a GenBank
  and a RefSeq trio. All 12 records are kept but flagged
  `duplicate_isolate=true` — for phylodynamic use, de-duplicate or thin.
- **Outgroup rooting is deferred**: v1 roots auto (`--root best`). A planned
  refinement roots on a California-serogroup relative (Snowshoe hare / Inkoo /
  Jamestown Canyon) or on the 181/57 complete genomes EU277663-665 — which,
  notably, are absent from the NCBI Datasets species-level pull.
- The 8 unassigned records ("OccaBUN"/"AeveBUN"/"Carynthia"/"NMGBT007") align
  to no TAHV reference; check whether they are TAHV before ever forcing them
  into a build.
- `Czechoslovakia` is kept as a country value (historically accurate; augur's
  built-in geolocation rules do not map it).
- Private-mutations QC gating becomes possible only after a first tree exists
  (build → nextclade with tree → re-gate → rebuild); divergence-band gating is
  the treeless proxy in use.
