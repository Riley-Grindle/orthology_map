# ORTHOLOGER: similarity thresholds in effect

The `ORTHOLOGER` process (`modules/local/orthologer.nf`) runs OrthoDB's
[OrthoLoger](https://orthologer.ezlab.org/) pipeline (v3.0.5, `quay.io/rgrindle/orthologer`)
with `ext.args = ''` (see `conf/modules.config`), i.e. no overrides — it runs entirely on
the pipeline's built-in defaults. This document records what those defaults actually are
and where they live, since none of it is visible from our own repo.

OrthoLoger has a single mode switch, `ODB_RUN_MODE`, with values `DEFAULT` / `XFAST` /
`SENSITIVE` / `USER`. We don't set it, so it falls back to `DEFAULT`. That mode fixes the
following parameters (from `/usr/local/ORTHOLOGER-3.0.5/bin/common_odb.sh`,
`set_run_mode()`, and confirmed against the generated `orthologer_conf.sh` for a fresh
project):

| Variable | Default value (DEFAULT mode) | What it controls |
|---|---|---|
| `MIN_OVERLAP` | `30` | Minimum alignment overlap / sequence length (residues) considered at all. |
| `MIN_COVERAGE` | `0.8` | Minimum fraction of a sequence's length that must be covered by the alignment for a hit to count. |
| `MAX_EVALUE` | `1.0e-3` | Max e-value allowed for a homology match. Passed straight to the aligner, e.g. `alignment_diamond.sh` calls `diamond ... -e $MAX_EVALUE`. |
| `SELECT_PID` | `97` | Min %identity used in the SELECT step, which collapses near-identical sequences within a species down to one representative before alignment. Intentionally close to 100% — this is a dedup step, not a similarity filter between orthologs. |
| `MMSEQS_SENSITIVITY` | `3.0` | mmseqs2 prefilter sensitivity. **Not actually in effect for us** — see below. |
| `MERGEINPAR_FCLUSTMIN` | `0.5` | Inparalog cluster merge: min fraction of a smaller cluster's members that must be linked to the bigger cluster to merge them. |
| `MERGEINPAR_FSPECMIN` | `0.4` | Inparalog cluster merge: min species overlap required between two clusters before merging. |
| `MAKEBRH_NALT` | `2` | "Fuzzy" best-reciprocal-hit: keep up to this many near-best candidates per gene, not just the single best. |
| `MAKEBRH_FSEPMAX` | `0.03` | Max score separation (fraction, 0–1) a fuzzy-BRH candidate may have from the best hit and still be kept. |

Two more clustering parameters are fixed at their hardcoded defaults (not touched by
`ODB_RUN_MODE`, but still worth knowing since they directly gate how aggressively
`BRHCLUS` merges clusters):

| Variable | Default | What it controls |
|---|---|---|
| `BRHCLUS_MAX_SPRATIO` | `0.5` | Max allowed species-overlap ratio between two clusters for a merge to be considered (`0.0` = no overlap allowed, `1.0` = overlap ignored). |
| `BRHCLUS_MERGE_FNLINKS` | `0.3` | Min fraction of extra cross-links (beyond the seed BRH pair) required to merge two clusters. |

## Which aligner is actually used

`orthologer_conf.sh`'s defaults (independent of `ODB_RUN_MODE`) are:

```
PREPROC_LABEL=CHECKFASTA
ALIGNMENT_LABEL=DIAMOND
MASKER_LABEL=SEGMASKER
SELECT_LABEL=CDHIT
CLUSTER_LABEL=BRHCLUS
```

So our run uses **DIAMOND** for the all-vs-all alignment step, not mmseqs2 or BLAST+.
That's why `MMSEQS_SENSITIVITY` (listed above) is a dead parameter for us — it only
applies when `ALIGNMENT_LABEL=MMSEQS`. The `DEFAULT`/`XFAST`/`SENSITIVE` run-mode
distinction otherwise only changes that one knob, so for a DIAMOND-based run all three
modes are actually identical except for name — the real cutoffs (`MIN_OVERLAP`,
`MIN_COVERAGE`, `MAX_EVALUE`, `SELECT_PID`, the merge parameters) are the same in all of
them.

## Net effect, in plain terms

A candidate hit between two sequences only becomes a usable match if:
- the aligned region is at least 30 residues long,
- it covers at least 80% of (one of?) the sequence(s)' length,
- DIAMOND's e-value for it is ≤ 1e-3.

Best-reciprocal-hit pairs are then clustered (`BRHCLUS`) with fairly permissive merge
settings (species-overlap ratio up to 0.5, only ~30% extra links required), and
in-paralogs are folded back into those clusters (`MERGEINPAR`) under similarly permissive
thresholds (cluster-link fraction ≥ 0.5, species-overlap ≥ 0.4). None of this is an
"ortholog similarity %" in the colloquial sense (e.g. "80% identical") — there is no
global percent-identity cutoff on orthologs themselves; `SELECT_PID=97` is a
same-species deduplication step, and the real ortholog-defining thresholds are the
e-value/coverage/overlap triad plus the cluster-merge ratios above.

## How to change these, if ever needed

Set `ODB_RUN_MODE=USER` plus the desired variables via `ext.args` in
`conf/modules.config`'s `withName: 'ORTHOLOGER'` block (OrthoLoger reads these as shell
env vars inside the container, so they'd need to be passed as `-e VAR=value`-style
environment, not CLI args to `orthologer -c run` — current module has no mechanism for
this and would need a small change to actually honor overrides).
