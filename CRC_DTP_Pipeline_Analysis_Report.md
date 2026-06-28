# CRC DTP Pipeline: Analysis Report and Methodological Rationale

## 1. Scientific framing

The pipeline interrogates a single hypothesis from three angles: that a transcriptomic
drug-tolerant persister (DTP) program carries prognostic and biological signal in
colorectal cancer (CRC). The DTP state is the slow-cycling, reversible phenotype a
subpopulation of tumour cells adopts under cytotoxic pressure such as 5-FU, and it is a
candidate substrate for relapse. If a DTP signature is real and clinically meaningful,
three independent predictions should hold, and each is tested by one module.

First, in primary CRC the signature should track patient outcome: tumours expressing more
of the program should relapse or die sooner. The `crc_survival` module tests this in two
independent CRC cohorts (GSE39582 microarray and TCGA-COAD RNA-seq). Second, if the signal
is a genuine cancer-stress program rather than a CRC-specific artefact, a directionally
consistent association should appear across many cancer types. The `pancan_survival` module
tests this across all eligible TCGA cohorts. Third, the program should be enriched in the
tissue compartment where persister biology is most active, namely metastatic deposits that
have survived selection. The `mets_de` module tests this in a paired Normal / Primary /
Metastasis design (GSE50760).

The first two modules are the same analysis class, signature score to survival, and they
overlap on TCGA-COAD. The third is a different class, tissue-stage differential expression
with no survival endpoint. This is why the project is one shared analytical core with three
modules rather than one statistical frame: forcing the metastasis differential-expression
analysis into a survival template, or scoring the shared COAD cohort two different ways,
would manufacture inconsistencies rather than resolve them.

## 2. The shared core

The core exists so that any choice made once is made everywhere. The single most important
consequence is that TCGA-COAD, which both survival modules touch, is now processed by
identical code: the same expression normalisation, the same clinical endpoint derivation,
the same statistical tests and gates. Before the merge, the two arms scored that cohort
differently, which is the main material inconsistency the audit surfaced.

`config.R` holds every shared constant (landmark cutpoints, gating thresholds, the tumour
sample filter, the recurrence endpoint, the random seed). `signatures.R` loads one master
`sig.csv` and lets each module name the columns it needs. `expression.R` converts each
platform to one convention, a log2-scale, gene-symbol-indexed matrix. `scoring.R` runs
ssGSEA and enforces one score-naming convention. `clinical.R` loads the TCGA Clinical Data
Resource and derives censoring-aware endpoints. `stats.R` holds the three tests, the effect
size, the gating logic and the false-discovery correction. `plotting.R` holds one theme,
one palette, and the three-way routing of figures by significance.

## 3. Module A: crc_survival

The module runs two cohorts through the same downstream logic and differs only where the
platform forces it. For GSE39582, Affymetrix probe intensities are mapped to gene symbols,
collapsed to one value per symbol, and left on the log2 scale the array delivers. The GEO
clinical table is parsed defensively: because GEO stores characteristics as positional
fields, the code first verifies that the field labels are identical across all samples
before trusting the hardcoded column map, and it reports any value that fails to parse to a
number rather than silently coercing it to missing.

For TCGA-COAD, primary tumours are queried, the STAR `tpm_unstrand` assay is taken,
transformed to log2(TPM+1), and collapsed to one value per symbol in log space. Clinical
data comes from the CDR for survival, from the GDC clinical record for treatment status, and
from the published colData for microsatellite-instability status.

Both cohorts then pass through the same endpoint derivation and the same battery of cohort
strata. The strata exist to isolate the populations where a persister signal would actually
be expected: treated patients, microsatellite-stable tumours, and specific stage groupings,
because in untreated or MSI-high disease the biological premise (selection under chemotherapy
in a stable genome) does not apply cleanly. Each stratum is scored for both overall survival
and recurrence with all three statistical tests, the results are corrected for multiple
testing, and figures are written into significant, non-significant or not-tested folders.
Subtyping (CMS/CRIS/PDS) has been removed from this module and lives standalone.

## 4. Module B: pancan_survival

This module generalises the CRC survival logic to the whole atlas in five phases. Phase one
loads the CDR and inventories every project, recording how many patients have usable overall
survival and recurrence data, and excludes the non-solid malignancies (LAML, DLBC) where the
DTP premise does not apply. Phase two processes each cohort independently and checkpoints to
disk, so a cohort that has already been scored is skipped on re-run; for each it downloads
expression, applies the identical TPM normalisation, checks that enough of the signature
genes are actually present before scoring, runs ssGSEA, derives the composite score, joins
clinical data, and only keeps the cohort if enough matched patients survive the join. Phase
three compiles the master table, derives the shared endpoints, and applies batch correction
across cancer types. Phase four runs statistics at two levels, and phase five renders the
figures and the cross-cohort forest plots.

The batch-correction logic deserves its own note because it is the one place the pan-cancer
analysis legitimately diverges from the single-cohort analysis. Cancer type is a genuine
batch for a pan-cancer score: baseline ssGSEA values differ between tissues for reasons that
have nothing to do with persister biology. So the aggregate pan-cancer association is computed
on batch-corrected scores, where cancer type is removed as a nuisance factor. The per-cohort
associations, by contrast, are computed on the uncorrected scores, because within a single
cohort there is no between-tissue batch to remove and correcting would only distort the
within-cohort variance. The output keeps both, clearly labelled.

## 5. Module C: mets_de

This module shares the core infrastructure (the signature loader, the ssGSEA wrapper and
naming, the theme, the seed, the cache) but keeps its own analytical core because it answers
a different question. The design is the 18-patient, three-tissue GSE50760 series, parsed and
validated to confirm it really is a clean Normal / Primary / Metastasis design before
proceeding.

Differential expression uses limma on log2(FPKM+1) with a paired model that includes patient
as a blocking factor, so each patient acts as their own control and inter-patient variation
does not leak into the tissue contrast. The central confound is addressed explicitly: a liver
metastasis sample is physically contaminated with hepatocytes, so any "metastasis-up" gene
might simply be a liver gene. The module quantifies hepatocyte content per sample with a
liver-specific gene-set score, validates that this score behaves as a contamination proxy
should (higher in liver metastases, tested with a paired non-parametric test), then refits
the differential-expression model with that score as a covariate and reports which genes
survive the adjustment. It also builds a purity-corrected expression matrix used downstream
for the single-sample diagnostics. Beyond differential expression it runs gene-set enrichment
on the ranked statistic, principal-component analysis to show separation by tissue, and a
single-sample ROC analysis asking whether each signature score alone can separate primary
from metastatic tissue in a paired comparison.

## 6. Review table: data sources

| Module | Cohort | Platform | What it tests | Why this design |
|---|---|---|---|---|
| crc_survival | GSE39582 | Affymetrix microarray | Signature vs OS/recurrence in CRC | Large, richly annotated CRC cohort with treatment and MMR status, enabling the stratified analyses where a persister signal is biologically expected |
| crc_survival | TCGA-COAD | RNA-seq (STAR) | Same, in an independent CRC cohort on a different platform | Cross-platform replication; processed identically to the pan-cancer COAD arm to remove the prior inconsistency |
| pancan_survival | 30+ TCGA projects | RNA-seq (STAR) | Whether the prognostic signal generalises across cancer types | Distinguishes a genuine pan-cancer stress program from a CRC-specific artefact |
| mets_de | GSE50760 | RNA-seq (FPKM) | Whether the program is enriched in metastasis vs primary | Paired tissue-stage design isolates the compartment where persister selection has occurred |

## 7. Review table: preprocessing and scoring

| Step | Method used | Logic | Alternative considered and why not |
|---|---|---|---|
| TCGA quantification | `tpm_unstrand` (log2(TPM+1)) | TPM is length-normalised, so ranks reflect relative transcript abundance, which is what a rank-based score needs | CPM is not length-normalised; because ssGSEA ranks genes within a sample, gene length reorders those ranks and changes the score. This was the material fix in the audit |
| ssGSEA input scale | log2 throughout; pan-cancer back-transform dropped | ssGSEA is a within-sample rank statistic, so monotonic transforms (log2, linear) give identical scores; standardising on log2 makes the inputs literally one convention and removes a pointless back-transform | Mixing linear and log2 inputs across modules was harmless numerically but obscured that the inputs were the same |
| Duplicate gene collapse (TCGA, array) | log-space `avereps` (geometric mean) | One value per symbol, computed consistently in the scale fed to ssGSEA; matches what the pan-cancer arm already did | Summing duplicates inflates highly-multi-mapped genes and is wrong for intensity/abundance data |
| Duplicate gene collapse (Mets FPKM) | summation | Cufflinks reports additive sub-features per gene id; their FPKM values are additive, so summing reconstructs the gene-level value | Averaging would understate genuine gene-level abundance for these additive sub-features |
| Composite score | Up minus Down, after scoring | A directional score is a score-level operation, not a gene set; computing it post-scoring keeps the signature file clean and the operation explicit | Encoding it as a third gene set conflates "what genes" with "how scores combine" |
| Signature coverage gate | require a minimum percent of signature genes present before scoring | A score computed from a fraction of its genes is not the same score; the gate prevents silently degraded scores in cohorts with poor coverage | Scoring regardless would produce non-comparable values across cohorts |

## 8. Review table: statistical methods (primary rationale table)

| Method | Where | What it answers | Why this method | Key assumption / caveat |
|---|---|---|---|---|
| ssGSEA | all modules | Per-sample signature activity | Rank-based and single-sample, so it needs no reference population and is robust to platform-specific scaling | Within-sample ranking; comparable across samples only after consistent normalisation |
| Mann-Whitney (Wilcoxon) | survival modules | Does score differ between outcome groups | Non-parametric; ssGSEA scores are not guaranteed normal and groups are unequal in size | Tests distributional shift, not a survival model; used as a screen alongside, not instead of, time-to-event analysis |
| Rank-biserial r | survival modules | Effect size for the Wilcoxon test | A p-value alone does not convey magnitude; r is the natural effect size for a rank test and is reported with sign relative to the first outcome level | Sign convention must be read against the reference level (Dead/Recurred) |
| Kaplan-Meier + log-rank | survival modules | Visual and test of separation over time | Standard, assumption-light way to display survival separation; the median split makes the figure legible | Median split is for visualisation only; it discards information and is not the inferential model |
| Cox proportional hazards, continuous score | survival modules | Effect of the score on hazard, per unit | Uses the full continuous score, avoids an arbitrary cutpoint, and retains statistical power; gives an interpretable hazard ratio with confidence interval and a concordance index | Assumes proportional hazards; landmarked at 36 months to keep the assumption reasonable |
| Landmark at 36 months, censoring-aware | survival modules | 3-year outcome with correct handling of incomplete follow-up | A patient censored before the landmark cannot be classified as an event or a non-event; the case_when logic returns missing rather than guessing, so no patient is silently misclassified | Patients with insufficient follow-up are excluded from that endpoint by design |
| Unified gating thresholds | survival modules | Whether a test is even admissible | Minimum group sizes and event counts prevent unstable estimates from tiny strata; applying them identically everywhere makes "not tested" mean the same thing in every module | Stricter gates reduce the number of reported tests |
| `removeBatchEffect` | pan-cancer aggregate only | Pan-cancer association net of tissue | Cancer type is a real batch for a cross-tissue score; removing it isolates the persister signal from baseline tissue differences | Applied only to the aggregate; per-cohort analyses use uncorrected scores because there is no between-tissue batch within a cohort |
| limma, paired design | mets_de | Differential expression Metastasis vs Primary | Designed for the small-n, many-gene setting with empirical-Bayes variance shrinkage; the paired term uses each patient as their own control | log2(FPKM+1) is not raw counts, so a count-based model is inappropriate here |
| Liver-purity covariate + negative control | mets_de | Is a metastasis-up gene real or hepatocyte contamination | Quantifies contamination directly and adjusts for it, then reports which genes survive; the negative-control validation confirms the proxy behaves correctly before it is trusted as a covariate | Adjustment can only correct for measured contamination, not unmeasured confounders |
| GSEA on ranked statistic | mets_de | Whether signatures are coordinately shifted | Uses the whole ranked gene list rather than an arbitrary significance threshold, so it detects coordinated modest shifts a cutoff would miss | Tie-breaking is the only stochastic step; fixed by the global seed |
| ROC / AUC, paired | mets_de | Can a single score separate tissue states | A threshold-free measure of separability that does not assume a model form | Paired structure respected; AUC summarises ranking, not calibration |
| Benjamini-Hochberg FDR | all modules | Control false discoveries across many tests | Controls the expected proportion of false positives, appropriate when many related hypotheses are tested and some signal is expected | The correction family (which tests are pooled) is an open decision; see section 11 |

## 9. Review table: cross-cutting design decisions

| Decision | Choice | Logic |
|---|---|---|
| One signature file | `sig.csv` holds all signatures; modules select by name | One source of truth; a signature edit propagates everywhere and cannot drift between modules |
| Shared recurrence endpoint | PFI for TCGA in both arms | PFI has the broadest coverage across TCGA and is the CDR-recommended progression endpoint; using one endpoint in both arms makes the shared COAD cohort directly comparable |
| Outcome labels | Dead_3yr / Alive_3yr, Recurred / Recurrence-Free | Self-documenting and identical across modules, so a figure or table reads the same everywhere |
| Primary tumour only | both survival arms | Defines the COAD overlap identically and removes the metastatic-sample fallback that previously differed between arms; keeps the survival populations biologically comparable |
| Continuous score for inference, split for display | Cox on continuous, median split for KM only | Preserves power and avoids an arbitrary threshold in the model while keeping the figure interpretable |
| Timestamped runs, global seed, sessionInfo | applied to every module | Reproducibility: nothing is overwritten, the only stochastic step is seeded, and the exact package versions are captured per run |
| Shared caching | GEO/GDC objects and per-cohort results | Slow downloads happen once and re-runs resume, which matters most for the 30-plus pan-cancer cohorts |
| Subtyping extracted | standalone file, not run by `main.R` | Its fate is undecided and three methodological questions are unresolved; isolating it prevents unreviewed calls from entering the analysis |

## 10. Assumptions and limitations

The survival analyses are landmark analyses at three years and therefore exclude patients
with shorter follow-up from the relevant endpoint; this is a deliberate trade of sample size
for correct classification. The Cox models assume proportional hazards, which the landmark
helps but does not guarantee. ssGSEA produces sample-relative scores, so all comparability
rests on consistent upstream normalisation, which is exactly why the normalisation was
unified. GSEA enrichment depends on the ranking metric and gene-set definitions. The
metastasis purity adjustment corrects only for measured hepatocyte contamination. Across all
modules the multiple-testing correction reduces but does not eliminate false positives, and
its grouping is not yet finalised.

## 11. Open methodological questions

Two decisions are deliberately left open and flagged in the code rather than resolved
silently.

The false-discovery correction family is undecided. The correction is currently applied with
the grouping that reproduces each module's original behaviour, so no result has changed: the
CRC arm corrects within dataset, cohort, test and metric; the pan-cancer arm corrects within
family and test. The real choice is whether to keep correcting within those narrow strata,
which gives more power per stratum, or to pool into a single project-wide family, which is
stricter and more defensible if the headline claim spans modules. The right answer depends on
how many independent tests actually underpin the dissertation's primary claims, which is the
sensible starting point for that conversation.

The molecular subtyping (CMS/CRIS/PDS) is extracted and not run. Before its calls can be
trusted, three issues need resolving: a centering asymmetry where one cohort was median-
centered for the classifier and the other was not; the risk of a double transform from
passing an already-log-transformed matrix to a classifier that applies its own RNA-seq
transform; and an undocumented normalisation expectation for the PDS classifier. These are
noted at the top of the subtyping file.
