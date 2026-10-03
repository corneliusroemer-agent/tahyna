# Ingest notes — tahyna

Milestone notes for `ingest/`, written when the ingest first validated end-to-end
(2026-10-02). Future fixes get noted here, not built.

## What runs

`ingest/Snakefile` (self-contained, all rules inline): cached-or-live NCBI
Datasets package → `dataformat` TSV + `unzip` FASTA → `augur curate passthru`
joins them into NDJSON → one streaming curate pipe (rename → normalize-strings
→ transform-strain-name → set_strain_display.py → format-dates →
normalize_segments.py → parse-genbank-location → titlecase → abbreviate-authors
→ apply-geolocation-rules → apply-record-annotations) → assignment-only
nextclade per segment reference (treeless) → `assign_segment.py` writes
`results/metadata.tsv` (+ `results/sequences.fasta` from curate).

Config knobs live in `ingest/defaults/config.yaml`; the taxon is the species
node 3052445 (not the virus node 45270 — 6 records visible only at the species
level). `live_fetch: false` uses `fixture/ncbi_dataset.zip` (121 records,
captured 2026-10-02); `--config live_fetch=true` refetches from NCBI.

## Validated-run numbers (fixture, 121 records)

- 121 records in = 121 curated; 121 unique `strain` ids; FASTA ids == metadata
  ids exactly.
- Segment: 54 S + 31 M (12 `M` + 19 `M-RNA` normalized) + 25 L by NCBI label;
  +3 assigned by nextclade coverage argmax (all 1.0000 coverage: 2 M
  [AF123485.1, AF229129.1], 1 S [U47142.1]); 8 left unassigned (all ≤ 294 nt
  fragments that align to no TAHV reference at coverage ≥ 0.5 — see below).
- Dates: 0 records without a date; 94 real (already `YYYY-MM-DD` or masked
  `YYYY-XX-XX`), 27 imputed (`date_imputed=true`) to 1957/1958/1962/1963/1964/
  1966/1968/2001. Span 1957 → 2024-07-18.
- Geo: country filled 115/121 (Czech Republic 37, Italy 22, China 22,
  Czechoslovakia 20, Austria 8, Slovakia 3, France 3, empty 6); region
  Europe 93 / Asia 22 / empty 6; division+location only 8 records have them
  (Austria 6, Inner Mongolia 2) — NCBI geo strings are mostly bare country.
- `is_lab_host=true`: 15 records, kept.
- `duplicate_isolate=true`: 12 records (2 isolates with more than one genome
  sequenced: Prototype '92' Bardos clones HM036208-213, XJ0625 GenBank+RefSeq
  trios) — their tips are not independent isolates.
- `segment_label_conflict`: 0 records — nextclade's coverage argmax never
  disagrees with an NCBI segment label in this dataset (the flag is wired and
  will fire if that changes with refetches).
- `year`: collection year as a discrete column (from `date`; imputed years
  included, `date_imputed` says which) for auspice categorical coloring.
- `abbr_authors` filled 121/121.
- Assemblies (exact `isolate` string, ≥2 records with ≥2 distinct segments):
  30 assemblies covering 86 records, incl. XJ0625 (6 records: GenBank + RefSeq
  trios) and Prototype '92' Bardos (6 records: HM036208-213).
- Wall time: full ingest ≈ 4 s on this container (3 nextclade runs dominate).

## Date imputation (the curated mapping)

`ingest/defaults/date_imputation.tsv` (3-column, headerless, comments allowed;
consumed by `augur curate apply-record-annotations`). Priority order:

1. Camp et al. 2021 (PMC8556901) Suppl. Table S5 isolate years, mapped at
   accession level — resolves the 19 GQ3868xx M segments to per-isolate years
   1958-1968. They are NOT one 2009 study year; 2009 is only the sequencing
   study's release year.
2. Strain-name year suffix verified in literature: 181/57 → 1957 (Růžek et al.
   2020, PMID 33020166 — isolated from a sick child, Czechoslovakia, 1957);
   prototype "92" Bardos → 1958 (Camp S5). EU185046.1 "Bardos 92-like" imputed
   1958 at lower confidence (a real 92-like field isolate would have its own
   unknown year).
3. Release-year last resort: only AX230490.1 (294 nt patent fragment, no
   isolate, no geo; released 2001). Release-year was NOT used for the old
   records — they were released 1996-2011, decades off.

Records with real dates keep `date_imputed` empty (NCBI boolean convention:
true when present, empty otherwise).

## Tip naming (display id)

`strain` is the display tip name = isolate name + accession, one token:
`XJ0625_EU622819`, `Prototype_92_Bardos_HM036208`, `AX230490.1` (no isolate →
bare accession). The pure, unmangled isolate name is kept in `isolate`. The
accession without version suffixes the name, so ids are unique by construction
(a version collision would fall back to the full versioned accession — not
triggered in this dataset; the script errors if it ever is).

Deviation 1 from the literal `XJ0625 (OP727994)` form: spaces/parentheses are
newick syntax and whitespace splits FASTA ids, so the same content is carried
with `_` instead. `isolate` preserves the readable form for tanglegrams and
export colorings.

## Lab-host divergence check (records kept, flagged, not dropped)

From the ingest nextclade runs: substitutions-vs-reference and coverage on each
record's assigned segment.

| segment | group | n | subs mean | subs max | cov mean | cov min |
|---|---|---|---|---|---|---|
| L | non-lab | 21 | 1135.4 | 1207 | 0.952 | 0.138 |
| L | lab | 4 | 100.8 | 228 | 1.000 | 1.000 |
| M | non-lab | 29 | 798.2 | 808 | 0.999 | 0.989 |
| M | lab | 4 | 250.8 | 838 | 1.000 | 1.000 |
| S | non-lab | 48 | 49.9 | 77 | 0.720 | 0.197 |
| S | lab | 7 | 29.0 | 67 | 0.937 | 0.561 |

Reading: lab-host records are not extra-divergent — the high non-lab L/M means
are the *European* strains' distance from the Chinese XJ0625 reference, while
every lab-host record happens to be a Chinese isolate (XJ0625/XJ0708/XJ0710,
NM08003/10). The one lab outlier, XJ0708's M (838 subs), is the known
anomalously divergent M segment (Lu 2011), not a lab artifact. Subs-vs-ref
mostly encodes geography here; judging "adaptation mutations" properly needs a
tree and ancestral states (later). Lab records ARE among the best-sequenced
(coverage 1.0 on L/M), consistent with lab-passaged high-quality isolates.

## Observations on the data (not blocking)

- The 8 unassigned records: 7 tiny N-gene-ish fragments (226-294 nt) plus
  GQ480358/HM068013-15 (the "OccaBUN"/"AeveBUN" mosquito-pool isolates) and
  KJ575081-2 ("Carynthia12a/b"), PX470108 ("NMGBT007") — all align to no TAHV
  segment reference (coverage 0.0 in all three runs). Worth checking whether
  they are TAHV at all or mislabelled/untyped orthobunyaviruses grouped under
  the species node; until then they stay out of every build (unassigned →
  no segment → filtered out), which is the safe direction.
- The actual 181/57 prototype complete genomes (EU277663-665, Růžek 2020) are
  absent from the datasets species-level pull — NCBI Virus portal gap, noted
  for a future fetch via Entrez if the prototype is wanted as an outgroup.
- `Czechoslovakia` is left as the country value (augur's built-in geolocation
  rules do not map it; mapping it to "Czech Republic" would be historically
  wrong for some isolates). Add a custom rule in
  `ingest/defaults/config.yaml` territory if ever needed.
- Year-precision dates appear in two spellings: real NCBI year-only dates are
  masked by format-dates to `2006-XX-XX`, imputed ones are bare `1958` (the
  annotation step runs after format-dates). augur's date handling treats both
  as year-precision ranges, so this is cosmetic only.

## Gotchas hit while building (for the next person)

- `augur curate rename --field-map` and other nargs-style flags
  (`--date-fields`, `--expected-date-formats`, `--backup-fields`,
  `--titlecase-fields`) must be word-split by the shell — quoting the joined
  string (`:q` in snakemake) turns it into one argument and dies with "must
  contain a single '=' character" / "not found in record".
- `augur curate apply-record-annotations` only OVERWRITES fields that already
  exist on the record; it silently skips (warning) fields that don't. Hence
  `set_strain_display.py --init-empty-field date_imputed --init-empty-field
  date_source` creates them empty before the annotations are applied.
- nextclade 3.21.2 writes a TSV row for sequences it FAILED to align (with
  empty metrics) in addition to warning "will not be included in the results" —
  contrary to what we saw in the 03 investigation (sequences absent from
  outputs). Empty coverage reads as 0, which is exactly what the argmax wants;
  just don't validate by row count here.
- `dataformat` Title Case headers verified: `Source database`,
  `BioSample accession`, `Is Lab Host`, `Submitter Names` — the rename map keys
  must match them verbatim.
- Snakemake + `workdir: workflow.current_basedir` moves the DAG's working
  directory, but `.snakemake/` (logs of the run itself) still lands in the
  invocation cwd — run from `ingest/` (`cd ingest && snakemake -s Snakefile`)
  to keep it in one place.
- nextclade leaves an `ncbi_sequences.fasta.fxi` index file next to its input;
  gitignored as junk.

## Not built (future fixes, noted not built)

- Private-mutations QC gating (needs a first tree; planned later refinement).
- Outgroup rooting with a California-serogroup relative (needs Entrez fetch of
  e.g. Snowshoe hare / Inkoo / Jamestown Canyon sequences, or EU277663-665).
- Custom geolocation rules (e.g. Czechoslovakia → modern names) if wanted.
- Fetching the 181/57 complete genomes (EU277663-665) into the dataset.
