---
name: nextstrain-pathogen-basic
description: Building a basic Nextstrain pathogen workflow end-to-end in this devcontainer (NCBI ingest + augur/nextclade phylogenetic builds, rotavirus/lassa archetype). Covers the two-Snakefile shape, the augur curate NDJSON pipe, per-segment build matrices with custom nextclade references, QC gating, and which output validations actually catch broken builds. Use whenever creating or debugging a Snakefile under a pathogen repo (ingest/Snakefile, phylogenetic/Snakefile), wiring `augur curate`, `nextclade run --input-ref`, `augur refine/export v2`, or diagnosing symptoms like "must contain a single '=' character", "Skipping annotation for field ... does not exist in record", "rate estimate is negative", or an auspice JSON with missing annotations.
---

# Basic Nextstrain pathogen workflow (ingest + phylo), the working shape

Fast path — a small virus (1-11 segments, < a few hundred records) is this,
and nothing more:

```
ingest/Snakefile:     cached-or-live ncbi_dataset.zip -> dataformat TSV + unzip FASTA
                      -> augur curate passthru (joins -> NDJSON)
                      -> one curate pipe -> results/{metadata.tsv,sequences.fasta}
                      -> nextclade per segment (assignment only) -> final segment column
phylogenetic/Snakefile: copy -> augur filter (segment + min_length) -> nextclade align+QC
                      -> QC gate -> tree -> refine -> ancestral -> translate -> traits
                      -> augur export v2 -> auspice/<prefix>_<segment>.json
```

Worked example (every claim below verified 2026-10-02/03, augur 34.0.0,
nextclade 3.21.2, snakemake 9.21.0): `scratch/tahyna/` — read its Snakefiles
before writing new ones. Reference material: the lassa ingest rules are human
authored and current; agent-written repos (e.g. scratch/rotavirus) are
secondary prior art.

## Archetype shape

- Two self-contained Snakefiles, all rules inline. Two, not one: ingest is a
  rare refetch, phylo is the iterate-fast loop; shared config would retrigger
  ingest on every phylo tweak. Skip rules/*.smk, vendored helpers, profiles,
  nextstrain-pathogen.yaml unless publishing to nextstrain.org.
- Each Snakefile: `workdir: workflow.current_basedir` (invocable from
  anywhere) + `configfile: "defaults/config.yaml"`. But `.snakemake/` still
  lands in the invocation cwd — run from the workflow dir.
- References are checked in: `references/<seg>/{reference.fasta,annotation.gff3,pathogen.json}`
  (plus reference.gb for humans). Bootstrap script under shared/, run once.
  pathogen.json is hand-curated and never overwritten by the bootstrap.

## augur curate pipe (NDJSON in, NDJSON out; only the last step writes files)

- Quoting: `--field-map`, `--date-fields`, `--expected-date-formats`,
  `--backup-fields`, `--titlecase-fields` are nargs-style — they MUST
  word-split. In snakemake, `{params.x:q}` on a joined string produces one
  argument and dies with `must contain a single '=' character` or
  `Expected date field '...' not found in record 0`. Emit the quotes inside
  the param (field-map pairs like `"Source database"=database`) and pass
  unquoted.
- `dataformat tsv virus-genome` headers are Title Case (`Source database`,
  `BioSample accession`, `Is Lab Host`); the rename map keys must match
  verbatim. Verify against the actual TSV, never from the mnemonics.
- `augur curate apply-record-annotations` only OVERWRITES fields that already
  exist on the record; new fields are skipped with `WARNING: Skipping
  annotation for field 'X' that does not exist in record`. Pre-create the
  fields empty in a tiny custom NDJSON step.
- `transform-strain-name --strain-regex '.+' --backup-fields accession` is the
  "never drop" idiom. Tip ids must be single tokens: whitespace/parens are
  newick syntax and split FASTA headers. Keep the pure isolate name in a
  second column (`isolate`) when the id needs mangling.
- Species-level taxids matter: NCBI classifies some records directly under the
  species node where the virus node cannot see them (TAHV: 121 vs 115).

## Per-segment builds

- `wildcard_constraints: segment = "s|m|l"`; per-segment values live in config
  (`label`, `min_length`, `ref_length`). Length floors ~50% of reference are a
  defensible default; set from the real length distribution.
- Custom nextclade references: `nextclade run --input-ref X --input-annotation
  Y --input-pathogen-json Z` (no `--input-dataset`, no server). One GFF per
  segment (both tools reject multi-seqid); GFF needs `##sequence-region` +
  gene rows with `gene=`/`Name=` attributes to also serve `augur translate`.
  pathogen.json needs `schemaVersion: "3.0.0"` and explicit numbers — bare
  `{"enabled": true}` gives degenerate QC (threshold 0).
- Treeless nextclade is fine for alignment+QC: clade/private-mutation columns
  come out null/zero. Private-mutation gating is impossible until a tree
  exists (proxy: substitutions-vs-reference band + coverage + frameshifts +
  stops). A QC gate = small python script joining the nextclade TSV to the
  metadata, writing kept-metadata + a dropped-strains report (never drop
  silently), then `augur filter --metadata kept --sequences aligned.fasta
  --output-sequences`.
- Segment labels from NCBI are dirty (`M-RNA`); normalize mechanically in the
  pipe, and let nextclade coverage-argmax assign blank labels above a floor.
  Flag (never override) conflicts between label and argmax.

## augur 34 specifics that bite

- `augur tree --method iqtree` passes a FIXED GTR (`-m GTR`) + `--ninit 2 -n 2
  --epsilon 0.05 -T AUTO --redo` — verified in augur/tree.py. nextstrain/ncov
  overrides only search intensity. Use defaults verbatim for small trees.
- `augur refine` fails with `rate estimate is negative` when the root-to-tip
  regression has no signal (deep between-lineage splits, temporally flat
  lineages). No rooting choice fixes a negative unrooted correlation. Ship the
  segment as a divergence tree via a per-segment config override rather than
  pinning an invented clock.
- refine and ancestral both want to write nt_muts.json — give them separate
  files (refine_<seg>.json keeps num_date/rate for export) and pass all node
  data files to export.
- export v2 config schema rejects `"version"` and `"metadata"`; the columns
  key is `metadata_columns`. Gene annotations land at
  `meta.genome_annotations` in the output (not top-level).
- `augur translate --reference-sequence` takes the same GFF3 nextclade used.

## Validation worth doing (catches real breakage)

1. Curated ids: record count in == out, unique ids, set(FASTA headers) ==
   set(metadata ids). Catches every join bug.
2. Per build: tips == QC-pass metadata rows == filter report; dropped-strains
   report non-empty ONLY with reasons.
3. Dates: 0 tips without dates if refine needs them; imputed flagged.
4. Timetrees: rates in 1e-5..1e-4 subs/site/yr for RNA viruses (>=1e-3 = red
   flag); report root-age CIs — with <50 tips they are often huge; say so.
5. Topology sanity against prior literature via patristic nearest-neighbours
   (Bio.Phylo): known lineages nested/shallow, known deep splits deep, known
   reassortment (one segment's neighbour disagrees with the others).
6. Auspice JSONs: keys tree/meta/version, tips == expected, genome_annotations
   populated, colorings count, root_sequence length == reference length.
