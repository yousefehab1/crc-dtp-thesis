# CRC DTP Pipeline

A single R project merging three previously separate analyses of drug-tolerant
persister (DTP) signatures in colorectal cancer. They share one core (loaders,
scoring, statistics, plotting) but stay as distinct modules because two of them
are the same analysis class and one is not.

- **`crc_survival`** (was `Final.R`) — ssGSEA → 3-year OS/RFS in two CRC cohorts
  (GSE39582 microarray, TCGA-COAD RNA-seq).
- **`pancan_survival`** (was `Panncan.R`) — the same ssGSEA → survival logic
  across all eligible TCGA cohorts, with pan-cancer batch correction.
- **`mets_de`** (was `Mets.R`) — a different class: paired tissue-stage
  differential expression (GSE50760 Normal/Primary/Metastasis), GSEA, PCA, and
  single-sample ROC. No survival.

`crc_survival` and `pancan_survival` are one analysis class and overlap on
TCGA-COAD; they now score it identically. `mets_de` shares infrastructure only.

## Layout

```
CRC_DTP_Pipeline/
  main.R                  orchestrator (sources core, runs the 3 modules)
  core/
    config.R             all shared constants / decisions
    io.R                 seed, timestamped run root, caching, sessionInfo
    signatures.R         single sig.csv loader
    expression.R         per-platform normalisation -> one ssGSEA input convention
    scoring.R            ssGSEA wrapper + score-naming convention
    clinical.R           CDR loader + censoring-aware endpoint derivation
    stats.R              Wilcoxon/KM/Cox, effect size, gating, FDR, runner
    plotting.R           one theme/palette, violin/KM/forest, 3-way routing
  modules/
    crc_survival.R
    pancan_survival.R
    mets_de.R
  subtyping/
    crc_subtyping.R      CMS/CRIS/PDS — PULLED OUT, not run by main.R (see #11)
```

## Inputs (place in the working directory)

- **`sig.csv`** — the master signature panel: one column per signature, gene
  symbols down the rows, blanks where columns differ in length. Every module
  reads this one file and selects the columns it needs (decision #1).
  Pan-cancer requires columns named `Up` and `Down` (it derives `Composite =
  Up − Down` after scoring). The Mets PCA overlay uses `METS_CORE_SIG`
  (default `UP`).
- **`TCGA-CDR.csv`** — TCGA Clinical Data Resource (Liu et al., 2018).

TCGA expression is downloaded via `TCGAbiolinks`; GEO series via `GEOquery`.
Downloads are cached under `cache/` so re-runs resume (decision #16).

## Run

```r
# from inside CRC_DTP_Pipeline/
Rscript main.R
```

Each run creates `CRC_DTP_YYYYMMDD_HHMM/` containing one subfolder per module,
plus `sessionInfo.txt`. Subtyping is run separately and on demand:

```r
source("core/config.R"); source("core/io.R"); source("core/expression.R")
source("subtyping/crc_subtyping.R")
run_crc_subtyping(init_run("subtyping"))
```

## Decisions applied (audit items 1–17)

| # | Decision | Where |
|---|----------|-------|
| 1 | One `sig.csv` holds all signatures; each module selects by name | `core/signatures.R` |
| 2 | TCGA arm uses `tpm_unstrand`; ssGSEA fed log2(TPM+1); back-transform dropped | `core/expression.R` |
| 3 | TCGA duplicates collapsed in log space (`avereps`); Mets FPKM still summed (additive sub-features) | `core/expression.R`, `modules/mets_de.R` |
| 4 | One shared CDR recurrence endpoint (`CDR_RFS_TYPE = "PFI"`) for both arms | `core/config.R`, `core/clinical.R` |
| 5 | Outcome labels `Dead_3yr` / `Alive_3yr` everywhere | `core/clinical.R` |
| 6 | Censoring-aware `case_when` endpoint derivation, else `NA` | `core/clinical.R` |
| 7 | Cox on the continuous score; median split kept for KM visualisation only | `core/stats.R` |
| 8 | Rank-biserial effect size, `Is_Testable`, 3-way routing (significant / non_significant / not_tested) | `core/stats.R`, `core/plotting.R` |
| 9 | Centralised gating thresholds applied identically | `core/config.R`, `core/stats.R` |
| 10 | **FDR strategy — FINALIZED (per-module grouping; see below)** | `core/stats.R` |
| 11 | CMS/CRIS/PDS subtyping pulled to `subtyping/crc_subtyping.R` | `subtyping/` |
| 12 | One score-naming convention `<Signature>_ssGSEA` | `core/scoring.R` |
| 13 | One ggplot theme + palette | `core/plotting.R` |
| 14 | Timestamped run root for all modules | `core/io.R` |
| 15 | One global `set.seed(42)` + `sessionInfo` capture | `core/io.R` |
| 16 | Shared caching / checkpoint discipline | `core/io.R`, `modules/` |
| 17 | Primary-tumour only in both CRC and pan-cancer (metastatic fallback removed) | `core/config.R`, `core/expression.R` |

### #10 — FDR strategy (finalized)

`apply_fdr()` takes an explicit `by` grouping. Each module keeps its own
grouping so that each preserves its original statistical behaviour; the two
groupings are intentionally different because the modules ask different
questions (single-cohort CRC vs. cross-cohort pan-cancer):

- `crc_survival`: `by = c("Dataset","Project","Test","Metric")`
- `pancan_survival`: `by = c("Family","Test")`

Trade-off of the correction family: per-cohort × test × metric (more power,
more tests overall) versus a single pooled family across the whole project
(stricter, fewer false positives). The adopted approach is the former for CRC
and Family × Test for pan-cancer, matching each module's original design.

## Open items (still need your decision)

**#11 — Subtyping.** Extracted and standalone. Before relying on its calls,
resolve the three notes at the top of `subtyping/crc_subtyping.R`: the
`ematAdjust` centering asymmetry between cohorts, the CMScaller RNA-seq flag vs.
input-scale double-transform risk, and the undocumented PDS normalisation.

## Required packages

CRAN: `dplyr`, `tidyr`, `tibble`, `stringr`, `readr`, `ggplot2`, `survival`,
`survminer`, `ROCR`, `msigdbr`.
Bioconductor: `limma`, `GSVA`, `GSEABase`, `GEOquery`, `Biobase`,
`SummarizedExperiment`, `AnnotationDbi`, `hgu133plus2.db`, `TCGAbiolinks`,
`clusterProfiler`, `enrichplot`, `BiocParallel`.
Subtyping only: `CMScaller`.
