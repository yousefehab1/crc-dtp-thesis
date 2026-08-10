# CRC DTP Pipeline: Code Quality & Architectural Review Improvement Plan

## 1. Executive Summary

This comprehensive code review evaluates the **CRC DTP Pipeline** (~5,000 lines across 19 R files and project documentation). The codebase demonstrates a strong architectural core: shared infrastructure (`core/`) centralises configuration, data I/O, ssGSEA scoring, clinical endpoint derivation, statistical gating, and publication-ready plotting. 

However, as the project expanded to accommodate pan-cancer analyses, subtyping, and extensive composite figures, several architectural debt patterns emerged. The overall codebase is in a functional and well-structured state, but suffers from significant presentation-layer bloat (most notably in `modules/composites.R`), subtle statistical inconsistencies between modules evaluating identical cohorts, and cache-invalidation risks.

### Top Highest-Impact Opportunities
1. **Harmonise Cross-Module Statistical & FDR Schemes**: Eliminate discrepancy in TCGA-COAD reporting between `modules/crc_survival.R` and `modules/pancan_survival.R` caused by differing FDR adjustment grouping scopes (`by` parameters in `apply_fdr()`).
2. **Deconstruct & Decouple `modules/composites.R`**: Restructure the 2,002-line monolithic presentation file into modular component scripts under `modules/composites/`, and shift per-SD hazard ratio pre-computation into `core/stats.R` to eliminate redundant raw data loading and SD re-estimations.
3. **Robustify Cache Keys for Global Configuration Changes**: Incorporate `ID_TYPE` and signature panel versions into `cache_rds()` keys to prevent silent loading of stale RDS files when switching gene ID namespaces (`symbol` vs `ensembl`).
4. **Fix Edge-Case Hazard Ratio Squishing**: Prevent `scales::squish()` in forest plot generation from converting `0` or `Inf` hazard ratios (resulting from non-convergent Cox fits in small subgroups) into valid boundary values (0.1 or 10.0).

---

## 2. Redundancy Findings

### Theme A: Duplicated Dataset, Endpoint, and Palette Dictionaries
- **Observation**: `DATASET_FULL`, `DATASET_SHORT`, and `ENDPOINT_LABEL` maps are defined in `core/plotting.R` (lines 32-38), but are re-declared or shadowed in `modules/composites.R` (lines 38, 786, 1274). `END_LAB <- c(OS = "3-yr OS", RFS = "3-yr RFS")` in `modules/composites.R` (line 786) duplicates `ENDPOINT_LABEL`.
- **Judgment**: Maintaining duplicate mapping vectors across files risks desynchronization if a dataset or endpoint label is renamed.
- **Verification**: Verified directly in `core/plotting.R:32-38` and `modules/composites.R:38, 786, 1274`.

### Theme B: Re-implementation of Statistical Helpers and Gating Thresholds
- **Observation**: `metric_cols()` is defined in `core/stats.R` (lines 7-10), but re-implemented locally in `modules/composites.R` at lines 321-322 (`metric_cols <- function(m) ...`) and lines 758-759.
- **Observation**: `MIN_KM_GROUP_N <- 5` and `MIN_EVENTS <- 5` are defined in `core/config.R` (lines 15-16), but redefined as local constants in `modules/composites.R` (lines 46-47).
- **Observation**: Standard deviation calculation for scaling continuous Cox log-HRs to per-SD HRs (`sd_of()` / `sd_vec`) is re-implemented 3 separate times in `modules/composites.R` (lines 339-351, 521-527, 761-768).
- **Judgment**: Redefining constants and helper functions leads to maintenance drift if thresholds in `core/config.R` are adjusted.
- **Verification**: Verified directly in `core/config.R:15-16`, `core/stats.R:7-10`, and `modules/composites.R:46-47, 321-322, 339-351, 521-527, 758-768`.

### Theme C: Duplicated Subtype Attachment and Data Loading Boilerplate
- **Observation**: Subtype attachment logic (`call_cms`, `call_pds`, `attach_subtypes`) is executed in `modules/crc_survival.R` (lines 90-98, 191-201) and duplicated in `subtyping/crc_subtyping.R` (lines 28-54).
- **Observation**: Loading and preprocessing TCGA STAR TPM matrices via `TCGAbiolinks` and `prep_tcga_tpm()` is implemented separately in `modules/crc_survival.R` (lines 133-149) and `modules/pancan_survival.R` (lines 239-258).
- **Judgment**: While `subtyping/crc_subtyping.R` is designed as a standalone entry point, sharing a high-level wrapper would prevent code drift between standalone subtyping and survival-integrated subtyping.
- **Verification**: Verified in `modules/crc_survival.R:90-98, 133-149`, `modules/pancan_survival.R:239-258`, and `subtyping/crc_subtyping.R:28-54`.

---

## 3. Unnecessary Complexity Findings

### Theme A: Monolithic Presentation Module (`modules/composites.R`)
- **Observation**: `modules/composites.R` is 2,002 lines long and handles 5 distinct CRC figure groups, pan-cancer figures (Fig 6, Fig 7, Appendix Forest, Table 3.5), and metastasis figures (GSEA and Liver Purity composites).
- **Observation**: The script implements `.dplyr_local()` (lines 88-96) to dynamically re-assign 14 `dplyr` functions into local environments to prevent S4 generic masking by `AnnotationDbi`/`GSEABase`.
- **Judgment**: A 2,000-line presentation script that re-parses raw CSVs, re-calculates standard deviations, re-fits batch corrections (`limma::removeBatchEffect` at lines 1380, 1471), and builds complex patchwork objects violates the Single Responsibility Principle.
- **Verification**: Verified directly in `modules/composites.R:1-2002`.

### Theme B: Deep Monolithic Execution Bodies in Analysis Modules
- **Observation**: `run_crc_survival()` in `modules/crc_survival.R` spans lines 12-360 (~350 lines), and `run_pancan_survival()` in `modules/pancan_survival.R` spans lines 166-414 (~250 lines).
- **Observation**: In `modules/crc_survival.R` (lines 35-50), positional field verification for GSE39582 characteristics (`characteristics_ch1.2` through `26`) uses a hardcoded loop inline within the main orchestrator function.
- **Judgment**: Long sequential functions with inline string parsing and data cleaning obscure the core analytical workflow and make unit testing or step-wise debugging difficult.
- **Verification**: Verified in `modules/crc_survival.R:12-360` and `modules/pancan_survival.R:166-414`.

### Theme C: Dual Rendering Passes for Titled vs. Publication Figures
- **Observation**: Every figure function in `modules/composites.R` builds `make_fig(TRUE)` (detailed titled version) and `make_fig(FALSE)` / `.strip_titles(make_fig(FALSE))` (clean publication version), triggering two full `ggplot` object construction passes per figure.
- **Judgment**: Executing dual plotting passes increases execution time without adding structural value; theme modification or title stripping can be applied to the pre-rendered plot object.
- **Verification**: Verified in `modules/composites.R:198-217, 267-285, 395-414, 456-471, 567-586, 628-638, 840-867, 1105-1138`.

---

## 4. Efficacy / Correctness Risk Findings

### Risk A: TCGA-COAD FDR P-Value Inconsistency Between Survival Modules
- **Fact**: `README.md` (lines 16-17) states that `crc_survival` and `pancan_survival` overlap on TCGA-COAD and score it identically. However:
  - In `modules/crc_survival.R` (line 232), FDR correction is applied using `by = c("Dataset", "Project", "Test", "Metric")`.
  - In `modules/pancan_survival.R` (line 386), FDR correction is applied using `by = c("Family", "Test")` where `Family` is `"PerCohort"`.
- **Judgment**: The reported `FDR_P` and `Is_Significant` values for TCGA-COAD differ between `modules/crc_survival.R` and `modules/pancan_survival.R` because the p-value adjustment family in `pancan_survival` includes all 24 TCGA cohorts (a larger pool of tests), whereas `crc_survival` adjusts only within TCGA-COAD and GSE39582. This creates conflicting statistical outputs for the same cohort in thesis tables.
- **Verification**: Confirmed in `modules/crc_survival.R:231-232` and `modules/pancan_survival.R:384-386`.

### Risk B: Silent Boundary Conversion of Failed or Infinite Cox Hazard Ratios
- **Fact**: In `modules/composites.R` (lines 346-348, 538-540, 1284-1288), continuous Cox hazard ratios are converted to per-SD hazard ratios via `exp(log(HR) * sd_sc)` and then passed to `scales::squish(hr_sd, FOREST_XLIM)`.
- **Judgment**: If a Cox model fails to converge or exhibits quasi-complete separation in a small subgroup (yielding `HR = 0` or `HR = Inf`), `log(HR)` produces `-Inf` or `Inf`. `scales::squish()` converts `-Inf` to `0.1` and `Inf` to `10.0`, silently plotting unstable or non-convergent estimates as valid extreme hazard ratios on the forest plot boundaries.
- **Verification**: Confirmed in `modules/composites.R:346-348, 538-540, 1284-1288`.

### Risk C: Probe Multi-Mapping Selection in Microarray Normalization
- **Fact**: `prep_microarray_symbols()` in `core/expression.R` (lines 20-28) maps Affymetrix probe IDs to gene symbols/Ensembl IDs using `AnnotationDbi::mapIds(..., multiVals = "first")`.
- **Judgment**: Using `multiVals = "first"` selects the first target ID returned by the annotation package for multi-mapping probes (~2.8% of probes on HG-U133 Plus 2.0). If probe-to-gene mapping order changes across AnnotationDbi releases, probe assignment can silently shift prior to duplicate collapsing with `limma::avereps()`.
- **Verification**: Confirmed in `core/expression.R:20-28`.

### Risk D: Missing Cache Invalidation on Namespace or Parameter Changes
- **Fact**: `cache_rds()` in `core/io.R` (lines 24-34) generates cache keys solely from the string key passed by the caller (e.g., `cache_rds("GSE39582_eset", ...)` or `cache_rds("tpm_TCGA-COAD", ...)`).
- **Judgment**: `config.R` defines `ID_TYPE <- "ensembl"` (line 46). If a user switches `ID_TYPE` to `"symbol"`, existing RDS files under `cache/` (which were built under `"ensembl"`) will be re-loaded without being invalidated, leading to a silent namespace mismatch between the cached expression matrix and the signature panel.
- **Verification**: Confirmed in `core/config.R:46`, `core/io.R:24-34`, and `core/expression.R:15-74`.

### Risk E: Inconsistent PDS Subtype Sample Filtering
- **Fact**: `harmonize_crc_modifiers()` in `core/stats.R` (lines 105-112) excludes `"Mixed"` PDS samples (`pds[!pds %in% paste0("PDS", 1:3)] <- NA`) for adjusted and interaction Cox models. However, `Subtype_Score_Stats.csv` (generated via `get_kruskal_stats()` in `modules/crc_survival.R:252`) includes `"Mixed"` samples when computing Kruskal-Wallis statistics across PDS.
- **Judgment**: In `modules/composites.R` (lines 927-929), the Kruskal-Wallis p-value displayed on the PDS violin plot (which omits `"Mixed"`) reflects a 4-group comparison including `"Mixed"`, while the plot displays 3 groups. This distinction is noted in line 1030, but represents an underlying inconsistency in dataset filtering across statistical tests.
- **Verification**: Confirmed in `core/stats.R:105-112`, `modules/crc_survival.R:252`, and `modules/composites.R:927-929, 1030`.

---

## 5. Performance Findings

### Performance Issue A: Repeated Disk Reads and Re-processing in `modules/composites.R`
- **Fact**: `modules/composites.R` reads `GSE39582_clinical.csv` and `TCGA_COAD_clinical.csv` from disk 5 separate times (lines 138, 323, 513, 756, 917).
- **Fact**: In lines 1380 and 1471, `modules/composites.R` re-evaluates `limma::removeBatchEffect()` on the 24 TCGA cohorts to compute corrected scores for pan-cancer violins and aggregate tables, duplicating work already performed during `run_pancan_survival()`.
- **Judgment**: Reading CSVs and re-running batch corrections inside the figure generation layer slows down composite figure regeneration and creates unnecessary I/O overhead.
- **Verification**: Confirmed in `modules/composites.R:138, 323, 513, 756, 917, 1380, 1471`.

### Performance Issue B: Uncached Full Organism Annotation Fetching
- **Fact**: `core/id_conversion.R` (lines 29-40) builds `.get_full_id_map` for all ~60,000 Entrez/Ensembl/Symbol keys in `org.Hs.eg.db` and caches it as `full_id_map_symbol_to_ensembl.rds`.
- **Judgment**: Building the full genome-wide map on the initial run takes 10-15 seconds and allocates significant memory, whereas querying `org.Hs.eg.db` for only the ~20,000 genes in the active dataset would be substantially faster.
- **Verification**: Confirmed in `core/id_conversion.R:29-40`.

### Performance Issue C: Full MSigDB Querying in `modules/mets_de.R`
- **Fact**: `modules/mets_de.R` (lines 101-107) queries `msigdbr::msigdbr(species = "Homo sapiens")` to extract `HSIAO_LIVER_SPECIFIC_GENES`.
- **Judgment**: Querying `msigdbr` fetches a ~100,000-row data frame on every run of `run_mets_de()`. Wrapping this lookup in `cache_rds()` or storing the liver gene set in `core/` eliminates redundant package operations.
- **Verification**: Confirmed in `modules/mets_de.R:101-107`.

---

## 6. `composites.R` Deep-Dive & Proposed Restructuring

### Architectural Assessment
`modules/composites.R` is the largest file in the repository (2,002 lines, 122 KB). It serves as a presentation layer, but currently performs data transformation, statistical calculation, and plot rendering in a single tightly-coupled file.

```
Current State (Monolithic):
modules/composites.R (2002 lines)
├── Part 1: CRC Composites (Groups 1-5, lines 42-1204)
├── Part 2: Pan-Cancer Composites (Fig 6, 7, Table 3.5, lines 1205-1741)
└── Part 3: Metastasis Composites (GSEA & Liver Purity, lines 1742-1955)
```

### Proposed Modular Architecture
Instead of maintaining a 2,000-line script that re-calculates standard deviations and re-reads CSVs, `composites.R` should be split into domain-focused figure builders under a dedicated `modules/composites/` directory:

```
Proposed State (Modular):
modules/composites.R (Orchestrator wrapper, ~60 lines)
└── modules/composites/
    ├── crc_composites.R        (CRC survival figure groups 1-5)
    ├── pancan_composites.R     (Pan-cancer forests, violins, Table 3.5)
    ├── mets_composites.R       (Metastasis GSEA & Liver purity figures)
    └── composite_helpers.R     (Shared layout, formatting, & theme utilities)
```

### Data Flow Optimization
1. **Shift per-SD HR calculation to `core/stats.R`**: Modify `run_all_stats()` and `get_cox_stats()` to record the cohort/subgroup score standard deviation (`Score_SD`) and per-SD hazard ratios (`HR_perSD`, `HR_lower_perSD`, `HR_upper_perSD`) directly in `CRC_Statistical_Summary.csv` and `FDR_Stats_Summary.csv`.
2. **Eliminate Raw Clinical Re-loading**: By embedding per-SD metrics in the statistical summary CSVs, figure builders in `pancan_composites.R` and `crc_composites.R` can draw forest plots directly from summary tables without re-loading `GSE39582_clinical.csv` or `TCGA_COAD_clinical.csv`.

---

## 7. Relationship to README.md ("Decisions Applied" & "Open Items")

### Confirmation & Alignment
- **Decision #1 (Single `sig.csv` panel)**: Confirmed in `core/signatures.R`. (Note: `scripts/verify_signature_sources.R:15` hardcodes `Data/Fatemeh.csv`; updating it to `SIG_FILE` ensures consistency).
- **Decisions #2 & #3 (TCGA log2 TPM & geometric mean `avereps`)**: Confirmed in `core/expression.R`.
- **Decisions #4 & #6 (CDR endpoint PFI & censoring-aware `derive_endpoints`)**: Confirmed in `core/clinical.R`.
- **Decision #7 (Continuous Cox primary, median-split KM visual)**: Confirmed in `core/stats.R`.
- **Decision #12 (Score naming `<Sig>_ssGSEA`)**: Confirmed in `core/scoring.R`.
- **Decision #17 (Primary-only TCGA sample filtering)**: Confirmed in `core/config.R` and `core/expression.R`.

### Clarifications & Additions to Documented Items

#### Open Item #10 (FDR Strategy)
- **Documented Status**: Marked as `PROVISIONAL` with separate call sites in `crc_survival` (`by = c("Dataset","Project","Test","Metric")`) and `pancan_survival` (`by = c("Family","Test")`).
- **Review Finding**: This split causes conflicting p-values for TCGA-COAD when compared across modules. 
- **Recommendation**: Formalise a two-tier FDR specification in `core/config.R`:
  1. `FDR_BY_COHORT <- c("Dataset", "Project", "Test", "Metric")` for primary single-cohort analyses.
  2. `FDR_BY_PANCAN <- c("Family", "Test")` for pan-cancer cross-cohort comparisons.
  Explicitly document in thesis methods that TCGA-COAD p-values in the pan-cancer section reflect multi-cohort adjustment.

#### Open Item #11 (Subtyping Integration)
- **Documented Status**: Extracted to `subtyping/crc_subtyping.R` with notes on `CMScaller` `RNAseq=FALSE` and `ematAdjust`.
- **Review Finding**: `core/subtyping.R` correctly implements `RNAseq=FALSE` to avoid double log-transformation. However, `modules/crc_survival.R` still invokes `call_cms()` and `call_pds()` inline (lines 92-93, 196-197) for subgroup Cox analysis. The documentation should clarify that subtyping code is shared via `core/subtyping.R`, while `subtyping/crc_subtyping.R` serves as the standalone driver.

---

## 8. Proposed Action Plan

The following prioritized steps provide a clear roadmap for refactoring and enhancing the codebase.

| Step | Action Item | Impact / Area | Risk | Effort | Target Files |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | **Harmonise TCGA-COAD FDR & Document Strategy**<br>Standardise FDR grouping keys or document two-tier FDR explicitly in `core/config.R`. | Correctness & Thesis Defensibility | Medium (changes COAD FDR values) | Small | `core/config.R`<br>`core/stats.R`<br>`modules/crc_survival.R`<br>`modules/pancan_survival.R` |
| **2** | **Add Namespace Cache Invalidation**<br>Incorporate `ID_TYPE` into `cache_rds()` key generation to prevent stale cache loading when switching gene ID namespaces. | Correctness & Data Integrity | Low | Small | `core/io.R`<br>`core/expression.R`<br>`core/id_conversion.R` |
| **3** | **Guard Forest Plot HR Squishing Against Infinite Values**<br>Check for `is.infinite(hr_sd)` or non-convergent Cox fits before calling `scales::squish()` to prevent plotting failed fits at axis boundaries. | Correctness & Visual Accuracy | Low | Small | `modules/composites.R` |
| **4** | **Export Per-SD Hazard Ratios to Statistical Summaries**<br>Compute `Score_SD` and `HR_perSD` in `core/stats.R` and write to CSV summaries, eliminating clinical CSV re-loading in composites. | Performance & Architecture | Low | Medium | `core/stats.R`<br>`modules/crc_survival.R`<br>`modules/pancan_survival.R`<br>`modules/composites.R` |
| **5** | **Restructure `modules/composites.R` into Sub-modules**<br>Split 2,002-line file into `modules/composites/` (`crc_composites.R`, `pancan_composites.R`, `mets_composites.R`, `composite_helpers.R`). | Maintainability & Code Quality | Low | Large | `modules/composites.R`<br>`modules/composites/*.R` |
| **6** | **Consolidate Duplicate Constants & Dictionaries**<br>Move `DATASET_FULL`, `ENDPOINT_LABEL`, and gating thresholds to `core/config.R` / `core/plotting.R`; remove redundant local declarations. | Redundancy & Maintainability | Low | Small | `core/config.R`<br>`core/plotting.R`<br>`modules/composites.R` |
| **7** | **Cache MSigDB Queries in `modules/mets_de.R`**<br>Wrap `msigdbr` fetch of liver-specific genes in `cache_rds()` to avoid downloading/loading full MSigDB data frames repeatedly. | Performance | Low | Small | `modules/mets_de.R` |
| **8** | **Update Signature Verification Script Input Path**<br>Update `scripts/verify_signature_sources.R` to read `SIG_FILE` (`Data/final.csv`) instead of hardcoded `Data/Fatemeh.csv`. | Code Hygiene | Low | Small | `scripts/verify_signature_sources.R` |
