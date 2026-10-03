# QC review — tahyna phylogenetic builds (CHECKPOINT)

> **VERDICT (user, 2026-10-02, applied):** the gate is approved exactly as
> proposed below and is wired in (`phylogenetic/bin/apply_qc_gate.py`,
> `phylogenetic/defaults/config.yaml` `qc_gate:`). PJ01_OP727996's L is kept
> via the config `exceptions` list (flags `qc_exception=true`; its frameshift
> stays visible in the tree annotations so artifact-vs-real can be judged from
> branch mutations). One deviation surfaced during tree building: **the S
> build ships as a divergence tree** — treetime's root-to-tip rate estimate is
> negative on every rooting tried (best / min_dev / skyline / rooting on each
> Asian tip; evidence in the "S timetree" section at the bottom). M and L
> timetrees built unpinned as locked.

The workflow is built and run through the nextclade align+QC stage
(`phylogenetic/Snakefile`: copy ingest outputs → `augur filter` per segment →
`nextclade run` with the custom reference/GFF/pathogen.json). Everything below
is observed on the real filtered sets (S n=39, M n=33, L n=24).

## 1. Proposed QC gating criterion (treeless-computable only)

Drop a record from a segment build if ANY of:

| # | criterion | exact rule | why |
|---|---|---|---|
| 1 | aligned coverage | `coverage < 0.5` | half-genome of aligned reference is the same bar the length floor sets |
| 2 | missing data | `qc.missingData.totalMissing > 500 (S) / 2250 (M) / 3500 (L)` | the pathogen.json `missingData` threshold (~50% of reference); belt-and-braces with #1 — catches N-heavy records that align |
| 3 | divergence outlier | `totalSubstitutions > 25% of reference length` (S > 244, M > 1122, L > 1745) | treeless proxy for the "private muts < N" gate: a mis-segmented or foreign sequence sits at ~25-40% divergence; observed maxima are ≤ 17.3% (XJ0708 M, real biology, kept) |
| 4 | frameshifts | `qc.frameShifts.totalFrameShifts > 0` | CDS-consistency artifact signal |
| 5 | premature stops | `qc.stopCodons.totalStopCodons > 0` | ditto |

Implementation shape (to be wired after review): a `filter_qc` rule reading
`results/nextclade_{seg}.tsv` + `results/filtered_{seg}_metadata.tsv`, writing
`results/qc_pass_{seg}_metadata.tsv` + `results/dropped_{seg}_qc.tsv` (strain,
accession, the failing metric(s) — nothing disappears silently), then an
`augur filter` pass restricting the aligned FASTA to the kept metadata. All
thresholds live in `phylogenetic/defaults/config.yaml` (`qc:` block), so
strictness is a config edit.

**Explicit proxy note (user-confirmed):** the usual `private muts < N` gate is
NOT implementable here — nextclade without a reference tree reports
`privateNucMutations.*` as null/zero (verified in
investigations/.../03-nextclade-custom.md). Criterion #3 (divergence band vs
reference) is the treeless stand-in. Once a first tree exists, a refinement
pass can do it properly: build tree → nextclade run WITH that tree → re-gate
on real private mutations → rebuild. Deferred until there is ground truth to
eyeball artifacts against.

## 2. Observed QC outcome distribution (real data, this run)

| segment | n (length floor passed) | coverage min–max | missing max | subs min–max | frameshifts > 0 | stops > 0 | overallStatus | would FAIL proposed gate |
|---|---|---|---|---|---|---|---|---|
| S (978 nt, floor 490) | 39 | 0.561–1.000 | 0 | 0–77 | 0 | 0 | 39 good | **0** |
| M (4489 nt, floor 2245) | 33 | 0.989–1.000 | 0 | 0–838 | 0 | 0 | 33 good | **0** |
| L (6979 nt, floor 3490) | 24 | 0.866–1.000 | 868 | 0–1207 | 1 | 0 | 23 good, 1 mediocre | **1** |

Length floors dropped: 16 S records (196-294 nt N-gene fragments), 0 M
(all ≥ 4437), 1 L (EU185046.1, 963 nt). The only proposed-gate failure:

- **PJ01_OP727996 (OP727996.1, L segment, complete 7022 nt, China 2016, Cai
  2023): 1 frameshift detected** (nextclade `overallStatus=mediocre`, frameshift
  score flagged). PJ01 was passaged 1x C6/36 + 2x BHK-21 before sequencing
  (literature note); the frameshift may be real (degraded ORF), a sequencing
  error, or an alignment edge. Reviewer decides: gate drops it (default), or
  we add an explicit keep-exception. Everything else passes every criterion.

Notes: `overallStatus` alone would pass the same 95/96 and drop PJ01_L too —
the hard-criteria gate is preferred because each drop carries a named reason
per the user's recipe. S coverage min 0.561 = the three 549-nt 1958-70s
partials, which the floors deliberately keep. L missing max 868 = HM243136
(prototype Bardos L, China CDC sequencing) — passes under the 50% threshold.

## 3. pathogen.json (verbatim, one per segment under references/)

All three share this shape; only `missingData.missingDataThreshold` and
`privateMutations.typical/cutoff` are per-segment (S/M/L). `schemaVersion:
"3.0.0"` is mandatory (parsing fails without it). privateMutations and
snpClusters are treeless no-ops today (kept, with real numbers, for the later
tree-based pass).

S (references/s/pathogen.json):

```json
{
  "qc": {
    "missingData": {
      "enabled": true,
      "missingDataThreshold": 500,
      "scoreBias": 50
    },
    "mixedSites": {
      "enabled": true,
      "mixedSitesThreshold": 10
    },
    "frameShifts": {
      "enabled": true
    },
    "stopCodons": {
      "enabled": true
    },
    "privateMutations": {
      "enabled": true,
      "typical": 25,
      "cutoff": 100,
      "weightLabeledSubstitutions": 4,
      "weightReversionSubstitutions": 6,
      "weightUnlabeledSubstitutions": 1,
      "weightLabeledDeletions": 1,
      "weightReversionDeletions": 1,
      "weightUnlabeledDeletions": 1
    },
    "snpClusters": {
      "enabled": true,
      "windowSize": 100,
      "clusterCutOff": 6,
      "scoreWeight": 50
    }
  },
  "schemaVersion": "3.0.0"
}
```

M (references/m/pathogen.json): identical except
`"missingDataThreshold": 2250`, `"typical": 150`, `"cutoff": 800`.

L (references/l/pathogen.json): identical except
`"missingDataThreshold": 3500`, `"typical": 200`, `"cutoff": 600`.

Threshold rationale: missingData = ~50% of each reference (the policy the
length floors implement); privateMutations typical/cutoff sized to observed
max divergence (S ≤ 7.9%, M ≤ 18.7% incl. XJ0708, L ≤ 17.3%) so real divergent
lineages are "mediocre" at most, never silently "bad"; mixedSites 10;
frameshifts/stopCodons bare-enabled (verified functional treeless in
03-nextclade-custom.md).

## 4. IQ-TREE model (verified, per request)

`augur 34.0.0 tree --method iqtree` passes a **fixed GTR** by default —
`augur/tree.py`: `DEFAULT_SUBSTITUTION_MODEL = "GTR"`, emitted as `-m GTR`,
plus `DEFAULT_ARGS["iqtree"] = "--ninit 2 -n 2 --epsilon 0.05 -T AUTO --redo"`.
ModelFinder only runs with `--substitution-model auto`. Canonical precedent
checked (not from memory): nextstrain/ncov `defaults/parameters.yaml` sets
only `tree-builder-args: "'-ninit 10 -n 4'"` (a heavier search for its big
tree) and overrides no model — i.e. it rides augur's fixed GTR too. **We use
augur's defaults verbatim** (fixed GTR, light search — right for ≤ 39-tip
alignments), recorded here as the matched precedent. No `--tree-builder-args`.

## 5. Locked settings this chain will use (for the record)

- refine: `--timetree --date-confidence --coalescent opt --date-inference
  marginal --root best --clock-filter-iqd 10` (IQD loosened from the default 4
  to keep historic/imputed tips in the clock; no `--stochastic-resolve` —
  polytomies kept).
- ancestral: `--inference joint`. aa annotation: `augur translate` with the
  same GFF3 (verified compatible), isolated behind one rule.
- traits: `--columns country region --confidence`. export v2 per segment with
  `--include-root-sequence-inline`, default colors.
- imputed dates: single timetree run per segment (no sensitivity pair).

## 6. S timetree — why the S build is a divergence tree (post-review finding)

With the locked unpinned refine, S fails: treetime raises
`calc_rate_susceptibility: rate estimate is negative` for `--root best`,
`--root min_dev`, `--coalescent skyline`, and rooting directly on each Asian
tip (XJ0625_EU622820, XJ0710_HM243142). Tip-name mangling was ruled out (all
39 newick tip names match the metadata exactly). The cause is in the data:

- midpoint-rooted root-to-tip regression slope on the S ML tree:
  **−1.2e-4 subs/site/yr** (negative);
- the S tree has a deep Asia-Europe split (Europe-Asia mean patristic
  distance 3.1x the within-Europe mean), and the European lineage is
  temporally flat: strains from 1957, 1958, 1963, 1966, 1984, 2019 and 2021
  all sit ~0.050-0.059 from the root — 64 years with almost no divergence
  accumulation in N;
- this matches the literature finding already recorded for TAHV: "no clear
  temporal clustering despite samples from 1958-2019".

M and L carry usable clocks on the same settings (see 05-build.md for the
numbers), so the override is per-segment (`refine.overrides.s.timetree:
false`). Revisit options for S later: external outgroup rooting, a pinned
clock at a literature rate, or more S sampling.

## 7. Outgroup rooting (2026-10-03, replaces --root best)

Lumbo virus (TAHV's closest relative, ~89% nt on S; the Calzolari 2022
precedent) RefSeq segments NC_043631.1 (S) / NC_043630.1 (M) / NC_043632.1 (L),
strain SAAr 1881, isolated **1962** in Argentina. Committed at
`fixture/outgroup/<segment>.fasta` with controlled tip name OUTGROUP.

Mechanics: `mafft --add --keeplength` joins the outgroup to the aligned
ingroup; iqtree + refine root on the OUTGROUP tip; export hides the tip via a
`{"nodes": {OUTGROUP: {"hidden": "always"}}}` node-data file instead of
pruning the tree (pruning would merge the outgroup's long branch into the
ingroup base and distort divergence display). Config switch:
`outgroup.enabled` in `phylogenetic/defaults/config.yaml`.

Observed:

- **S still does not timetree.** The rate estimate stays negative with the
  outgroup rooted - the European lineage is temporally flat (1957-2021 at
  near-constant divergence), and the 1962-dated outgroup is not older than the
  old TAHV tips, so it cannot anchor the regression (a 1962 tip at ~11%
  divergence would actively steepen the slope the wrong way). S ships as an
  OUTGROUP-ROOTED divergence tree (`refine.overrides.s.timetree: false` ->
  `--root OUTGROUP` only).
- M and L timetrees run with the outgroup present; **treetime prunes the
  undated outgroup tip itself** ("pruning leaf OUTGROUP" in the logs), so the
  M/L trees display ingroup-only without extra work. The outgroup still shaped
  those trees: rates shifted up from the no-outgroup build
  (M 5.8e-5 -> 1.0e-4, at the top of the plausible 1e-5..1e-4 band; L
  5.3e-5 -> 7.6e-5) because the outgroup-inclusive alignment changes the tree
  the clock is fit on. Both remain in-band; the shift is reported, not hidden.
- Root placement verified: the S tree root's children are [OUTGROUP, ingroup];
  M/L roots sit at the outgroup attachment point (ingroup MRCA).

Note for anyone retrying an S timetree: pinning a clock (`--clock-rate`) is
excluded by the standing "unpinned" decision; an outgroup OLDER than 1957
would be the other lever, and none is available (Lumbo 1962 is the closest
relative with complete genomes).
