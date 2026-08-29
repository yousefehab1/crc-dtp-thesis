# CRC DTP Pipeline

A single R project merging three previously separate analyses of drug-tolerant
persister (DTP) signatures in colorectal cancer. They share one core (loaders,
scoring, statistics, plotting) but stay as distinct modules because two of them
are the same analysis class and one is not.

- **`crc_survival`** (was `Final.R`) — ssGSEA → 3-year OS/RFS in two CRC cohorts
  (GSE39582 microarray, TCGA-COAD RNA-seq), plus CMS/PDS molecular subtyping and
  the confounding/effect-modification analysis built on top of it.
- **`pancan_survival`** (was `Panncan.R`) — the same ssGSEA → survival logic
  across all eligible TCGA cohorts, with pan-cancer batch correction.
- **`mets_de`** (was `Mets.R`) — a different class: paired tissue-stage
  differential expression (GSE50760 Normal/Primary/Metastasis), GSEA, PCA, and
  single-sample ROC. No survival.

`crc_survival` and `pancan_survival` are one analysis class and overlap on
TCGA-COAD; they now score it identically. `mets_de` shares infrastructure only.

## Layout

```
Thesis/
  main.R                  orchestrator (sources core, runs the 3 modules + composites)
  core/
    config.R             all shared constants / decisions
    io.R                 seed, timestamped run root, caching, sessionInfo
    id_conversion.R      symbol <-> Ensembl ID mapping (ID_TYPE)
    signatures.R         signature-panel loader (SIG_FILE)
    expression.R         per-platform normalisation -> one ssGSEA input convention
    scoring.R            ssGSEA wrapper + score-naming convention
    clinical.R           CDR loader + censoring-aware endpoint derivation
    stats.R              Wilcoxon/KM/Cox, effect size, gating, FDR, runner
    plotting.R           one theme/palette, violin/KM/forest, 3-way routing
    subtyping.R          CMS + PDS callers, shared by crc_survival and subtyping/
    gsea_plots.R         GSEA enrichment plot primitive (used by mets_de)
  modules/
    crc_survival.R
    pancan_survival.R
    pancan_treated.R      standalone pan-cancer treated-only section (§3.5-treated)
    mets_de.R
    composites.R          composite/publication figures for all modules
  subtyping/
    crc_subtyping.R       standalone CMS/PDS driver (subtype calls only, no survival)
  scripts/
    verify_signature_sources.R   ad-hoc MSigDB overlap check for the panel's gene sets
```

## Inputs (place in the working directory)

- **Signature panel** (`SIG_FILE` in `core/config.R`) — one column per
  signature, gene symbols/IDs down the rows, blanks where columns differ in
  length. Every module reads this one file and selects the columns it needs
  (decision #1). Pan-cancer requires columns named `Up` and `Down` (it derives
  `Composite = Up − Down` after scoring). The Mets PCA overlay uses
  `METS_CORE_SIG` (default `Up`).
  **This file is per-sample confidential data and is not committed to this
  repository** (see `.gitignore`) — obtain it separately and place it at the
  path `SIG_FILE` points to.
- **`Data/TCGA-CDR.csv`** — TCGA Clinical Data Resource (Liu et al., 2018).
  Included in the repo.

TCGA expression is downloaded via `TCGAbiolinks`; GEO series via `GEOquery`.
Downloads are cached under `cache/` so re-runs resume (decision #16).

## Run

```r
# from inside Thesis/
Rscript main.R
```

Each run creates `CRC_DTP_YYYYMMDD_HHMM/` containing one subfolder per module
(including molecular subtyping, run as part of `crc_survival`), the composite
publication figures, and `sessionInfo.txt`. To get subtype calls alone,
without the survival analysis, run the standalone driver:

```r
source("core/config.R"); source("core/io.R"); source("core/expression.R")
source("core/scoring.R"); source("core/subtyping.R")
source("subtyping/crc_subtyping.R")
run_crc_subtyping(init_run("subtyping"))
```

## Decisions applied (audit items 1–17)

| # | Decision | Where |
|---|----------|-------|
| 1 | One signature panel (`SIG_FILE`) holds all signatures; each module selects by name | `core/signatures.R` |
| 2 | TCGA arm uses `tpm_unstrand`; ssGSEA fed log2(TPM+1); back-transform dropped | `core/expression.R` |
| 3 | TCGA duplicates collapsed in log space (`avereps`); Mets FPKM still summed (additive sub-features) | `core/expression.R`, `modules/mets_de.R` |
| 4 | One shared CDR recurrence endpoint (`CDR_RFS_TYPE = "PFI"`) for both arms | `core/config.R`, `core/clinical.R` |
| 5 | Outcome labels `Dead_3yr` / `Alive_3yr` everywhere | `core/clinical.R` |
| 6 | Censoring-aware `case_when` endpoint derivation, else `NA` | `core/clinical.R` |
| 7 | Cox on the continuous score; median split kept for KM visualisation only | `core/stats.R` |
| 8 | Rank-biserial effect size, `Is_Testable`, 3-way routing (significant / non_significant / not_tested) | `core/stats.R`, `core/plotting.R` |
| 9 | Centralised gating thresholds applied identically | `core/config.R`, `core/stats.R` |
| 10 | **FDR strategy — finalized.** Each module keeps its own per-analysis grouping (narrow families maximise power per stratum); there is no single project-wide FDR number, so the dissertation's primary confirmatory claim is defined narrowly (Up/Down/Composite vs. 3-yr OS/RFS in the main CRC cohorts) and everything built on subtype calls is reported as exploratory. | `core/stats.R` |
| 11 | **CMS/PDS subtyping — integrated.** Shared callers in `core/subtyping.R` (both cohorts classified as already-log2, `RNAseq=FALSE`, removing the double-transform/centering-asymmetry risk); used inline by `crc_survival` for per-subtype survival cohorts and the confounding/interaction Cox models, and by the standalone `subtyping/crc_subtyping.R` driver. CRIS dropped in favour of PDS. | `core/subtyping.R`, `modules/crc_survival.R`, `subtyping/crc_subtyping.R` |
| 12 | One score-naming convention `<Signature>_ssGSEA` | `core/scoring.R` |
| 13 | One ggplot theme + palette | `core/plotting.R` |
| 14 | Timestamped run root for all modules | `core/io.R` |
| 15 | One global `set.seed(42)` + `sessionInfo` capture | `core/io.R` |
| 16 | Shared caching / checkpoint discipline (cache keys that depend on `ID_TYPE` embed it explicitly, e.g. `tpm_<project>_<ID_TYPE>`) | `core/io.R`, `modules/` |
| 17 | Primary-tumour only in both CRC and pan-cancer (metastatic fallback removed) | `core/config.R`, `core/expression.R` |

Full methodological rationale, assumptions, and limitations (landmark analysis
trade-offs, FDR framing, PDS retained as exploratory, etc.) are written up in
the thesis itself — see `writing/Section_2_Methods.md` onward.

## Required packages

CRAN: `dplyr`, `tidyr`, `tibble`, `stringr`, `readr`, `ggplot2`, `survival`,
`survminer`, `ROCR`, `msigdbr`.
Bioconductor: `limma`, `GSVA`, `GSEABase`, `GEOquery`, `Biobase`,
`SummarizedExperiment`, `AnnotationDbi`, `hgu133plus2.db`, `TCGAbiolinks`,
`clusterProfiler`, `enrichplot`, `BiocParallel`.
Subtyping: `CMScaller`, `PDSclassifier` (GitHub-only; pinned commit — see
`subtyping/crc_subtyping.R`). If `PDSclassifier` is absent, PDS is skipped and
CMS + survival still run.
