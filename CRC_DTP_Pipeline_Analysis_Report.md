# CRC DTP Pipeline: Analysis Report and Methodological Rationale

*Last updated: 2026-07-03 — molecular subtyping (CMS + PDS) is now integrated into the
`crc_survival` module, and the module additionally tests whether CMS / PDS / Stage / MSI
confound or modify the DTP signature's survival associations (adjusted + interaction Cox).
A pre-launch methodological review settled seven open decisions: the FDR family (narrow,
per-analysis families kept; primary claim vs exploratory tests defined), Stage and MSI
binning (MSI-L folded into MSS), the recurrence endpoint (PFI, not DFI), the
events-per-variable gate (raised to 10), the confounding decision rule (descriptive), and
PDS (retained as exploratory). See sections 3, 8, 9, 10, 11.*

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
one palette, and the three-way routing of figures by significance. `subtyping.R` holds the
download-free molecular subtype callers (CMS and PDS) so the survival module and the
standalone subtyping script share one implementation rather than duplicating it.

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
Violin figures show every patient rather than a subsample: the violins are width-scaled so
the shapes stay comparable across unequal group sizes, and the per-group n is annotated on
each panel. An earlier balanced-subsampling step (downsampling the majority group to the
minority size for display) was removed because it discarded data and made the density
estimates noisier without improving the comparison — width-scaling already handles the
imbalance, and the statistics were always computed on the full data regardless.

Molecular subtyping is now folded into this module rather than living standalone. Both
cohorts are classified into consensus molecular subtypes (CMS1–4) and pathway-derived
subtypes (PDS1–3, plus a threshold-aware "Mixed" call) from the same expression matrices the
module already builds, so no data is re-downloaded. The calls are attached to each patient
and used two ways. First, per-subtype survival cohorts (for example `GSE_CMS1`, `TCGA_PDS2`)
are added to the existing strata and pass through the identical three-test, gated,
FDR-corrected machinery, asking whether the DTP signal predicts outcome within a subtype.
Second, a score-across-subtype comparison tests whether the DTP score itself differs between
subtypes, using a Kruskal-Wallis test with an epsilon-squared effect size and a dedicated
multi-group violin. The classifier input scale is handled explicitly to avoid a
double-transform (see section 11): both cohorts are classified as already-log2 expression so
the CMS classifier applies identical centering to each. CRIS was dropped in favour of PDS.

Building on those subtype calls, the module then asks whether the signature's survival and
recurrence associations are confounded or modified by CMS subtype, PDS subtype, disease stage,
or MSI status. Two model families are fitted on the core DTP scores (Up, Down, Composite).
Adjusted Cox models add each factor to the score and report the score's adjusted hazard ratio,
its change from the unadjusted estimate on the same patients, and a likelihood-ratio test of
whether the score still adds prognostic value beyond clinicopathology — this is the confounding
question. The change-in-estimate is reported descriptively (the absolute shift in the score's
log-hazard-ratio and a floored percent version), with the conventional ~10% change read as a
guideline for confounding rather than an automatic flag (Greenland 1989; Rothman). Interaction
Cox models test `score × factor` against `score + factor` with a likelihood-ratio test, one
factor at a time, and render the per-subgroup score hazard ratios as a forest — this is the
effect-modification question. Both families are gated by an events-per-variable rule (≥10
events per model parameter; Peduzzi 1996) so underpowered subgroup models — chiefly the
higher-parameter score×CMS interactions in the smaller cohort — are marked not-testable rather
than fitted unstably. Results are written to `Adjusted_Cox_Summary.csv` and
`Interaction_Cox_Summary.csv`.

The four modifiers are harmonised into clean, cohort-consistent factors before modelling.
Stage is dichotomised into Early (I/II) versus Late (III/IV), the standard node-negative
versus node-positive/metastatic prognostic split. MSI status is coded as MSI-high (dMMR)
versus everything else: MSI-low is grouped with MSS rather than with MSI-high, because MSI-low
lacks the hypermutator and immune phenotype that defines the distinct MSI-high entity and
behaves like MSS (Rantanen 2023 pooled MSS/MSI-low against MSI-high); this also makes
GSE39582, which only carries dMMR/pMMR, consistent with TCGA. MSI is a comparatively weak
modifier in TCGA-COAD, where it is missing for roughly 56% of patients. These confounder and
effect-modification analyses, together with the per-subtype survival cohorts and the
score-across-subtype comparison, are reported as **exploratory / hypothesis-generating**
rather than confirmatory (see sections 10–11).

## 4. Module B: pancan_survival

This module generalises the CRC survival logic to the whole atlas in five phases. Phase one
loads the CDR and inventories every project, recording how many patients have usable overall
survival and recurrence data, and excludes the non-solid malignancies (LAML, DLBC) where the
DTP premise does not apply. Both the pre-queue exclusions (with their reason: manually
excluded or below the patient threshold) and the per-cohort run outcomes are reported to the
console and written to CSV, so it is always explicit which cancer types entered the analysis
and why the rest did not. Phase two processes each cohort independently and checkpoints to
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
| Kruskal-Wallis + epsilon-squared | crc_survival (score-across-subtype) | Does the DTP score differ across molecular subtypes | Non-parametric extension of Wilcoxon to the 3–4 subtype groups (CMS1–4, PDS1–3); epsilon-squared gives a rank-based effect size | Tests distributional shift across groups, not which pair differs; underpowered subtypes are gated out per level |
| Kaplan-Meier + log-rank | survival modules | Visual and test of separation over time | Standard, assumption-light way to display survival separation; the median split makes the figure legible | Median split is for visualisation only; it discards information and is not the inferential model |
| Cox proportional hazards, continuous score | survival modules | Effect of the score on hazard, per unit | Uses the full continuous score, avoids an arbitrary cutpoint, and retains statistical power; gives an interpretable hazard ratio with confidence interval and a concordance index | Assumes proportional hazards; landmarked at 36 months to keep the assumption reasonable |
| Adjusted (multivariable) Cox | crc_survival (confounding) | Does the score's HR survive adjustment for CMS / PDS / Stage / MSI, and does the score add signal beyond clinicopathology | Separate adjustment models (clinicopathology = Stage+MSI; then +CMS; then +PDS) rather than one collinear model; reports the score's adjusted HR, the change-in-estimate vs the unadjusted HR on the same subset, and a likelihood-ratio test of score-added value | CMS is collinear with MSI and PDS, so they are not co-adjusted; models gated by events-per-variable |
| Interaction Cox + LRT | crc_survival (effect modification) | Does the score's effect differ across subgroups of CMS / PDS / Stage / MSI | Likelihood-ratio test of `score * modifier` vs `score + modifier`, one modifier at a time; complemented by per-subgroup score HRs shown as a forest | Interaction tests are low-powered; underpowered models are gated out (events-per-variable) rather than reported as null |
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
| Shared recurrence endpoint | PFI for TCGA in both arms | PFI has the broadest coverage across TCGA and is the CDR-recommended progression endpoint; using one endpoint in both arms makes the shared COAD cohort directly comparable. DFI was considered and rejected: only 190 COAD patients have it and it carries just 24 events, which would collapse the stratified, subtype and interaction models. The trade-off is a mild cross-cohort concept mismatch — GSE39582 uses true relapse-free survival, whereas TCGA PFI counts progression, new tumour, or death — accepted so the TCGA recurrence arm retains usable power |
| Outcome labels | Dead_3yr / Alive_3yr, Recurred / Recurrence-Free | Self-documenting and identical across modules, so a figure or table reads the same everywhere |
| Primary tumour only | both survival arms | Defines the COAD overlap identically and removes the metastatic-sample fallback that previously differed between arms; keeps the survival populations biologically comparable |
| Continuous score for inference, split for display | Cox on continuous, median split for KM only | Preserves power and avoids an arbitrary threshold in the model while keeping the figure interpretable |
| Timestamped runs, global seed, sessionInfo | applied to every module | Reproducibility: nothing is overwritten, the only stochastic step is seeded, and the exact package versions are captured per run |
| Shared caching | GEO/GDC objects and per-cohort results | Slow downloads happen once and re-runs resume, which matters most for the 30-plus pan-cancer cohorts |
| Subtyping integrated | CMS + PDS folded into `crc_survival`, shared callers in `core/subtyping.R` | The classifier scale questions are resolved (log2 input, `RNAseq=FALSE`, symbol-keyed matrices); subtype calls now feed both per-subtype survival cohorts and a score-across-subtype comparison. A standalone runner is kept for subtype calls without survival |
| Subtyping input namespace | dedicated log2 gene-symbol matrix, independent of the ssGSEA `ID_TYPE` | CMScaller and PDS expect gene symbols, whereas ssGSEA runs on the configured ID namespace (ensembl); building a separate symbol matrix keeps each tool on its required identifiers without changing the scoring path |
| Confounding vs modification kept separate | adjusted Cox for confounding, interaction Cox for effect modification | They are different questions and conflating them misleads; an adjusted HR answers "does the signal survive control for this factor", an interaction answers "does the signal differ by this factor" |
| Separate adjustment models, not one | Stage+MSI, then +CMS, then +PDS individually | CMS is collinear with MSI (CMS1↔MSI) and with PDS (PDS3⊂CMS2); a single all-in model would be unstable, so each factor's confounding contribution is isolated |
| Events-per-variable gate | require events ≥ `MIN_EPV` (=10) × parameters for adjusted/interaction Cox | Multivariable and interaction models overfit when events are scarce; 10 events per parameter is the conventional Cox threshold (Peduzzi 1996). The confirmatory whole-cohort models clear it comfortably; the gate mainly marks the higher-parameter score×CMS interactions in the smaller cohort not-testable rather than reporting unstable estimates |
| MSI-L grouped with MSS | MSI = MSI-high/dMMR only; MSS = MSS/pMMR/MSI-low | MSI-low lacks the hypermutator/immune phenotype that defines the distinct MSI-high entity and behaves like MSS (Rantanen 2023 pooled MSS/MSI-low vs MSI-high); also makes GSE39582 (dMMR/pMMR only) consistent with TCGA |
| Confounding read descriptively | report adjusted vs unadjusted HR and the change-in-estimate; no automatic flag | Change-in-estimate is a heuristic without an inferential test and the percent metric is unstable when the unadjusted HR is near 1; reporting the absolute log-HR shift with the ~10% rule as a stated guideline (Greenland 1989) is more honest than a binary flag |
| Subtype/confounder tests exploratory | per-subtype survival, Kruskal-Wallis, adjusted and interaction Cox labelled hypothesis-generating | The primary confirmatory claim is the DTP scores→OS/RFS association in the main CRC cohorts; the secondary tests are numerous and underpowered, so they are framed as exploratory rather than pooled into the confirmatory error rate |

## 10. Assumptions and limitations

The survival analyses are landmark analyses at three years and therefore exclude patients
with shorter follow-up from the relevant endpoint; this is a deliberate trade of sample size
for correct classification. The Cox models assume proportional hazards, which the landmark
helps but does not guarantee. ssGSEA produces sample-relative scores, so all comparability
rests on consistent upstream normalisation, which is exactly why the normalisation was
unified. GSEA enrichment depends on the ranking metric and gene-set definitions. The
metastasis purity adjustment corrects only for measured hepatocyte contamination. The
confounding and effect-modification analyses in `crc_survival` are constrained by power: the
interaction (effect-modification) tests in particular need many events to detect a subgroup
difference, so a non-significant interaction is not evidence of a uniform effect, and the
events-per-variable gate will legitimately leave some subtype models untested in the smaller
cohort. Adjustment removes only the measured factors (CMS, PDS, Stage, MSI), not unmeasured
confounders. The false-discovery correction is applied within narrow per-analysis families
(see section 11), so significance is a per-table statement rather than a single study-wide
error rate; the subtype, confounder and interaction analyses are read as exploratory for
exactly this reason.

## 11. Open methodological questions

A pre-launch methodological review resolved the decisions that were previously left open;
what remains is a small number of framing choices that are now settled and documented here.

The false-discovery correction keeps its narrow, per-analysis families. Each module corrects
within the grouping that reproduces its original behaviour — the CRC arm within dataset,
cohort, test and metric; the pan-cancer arm within family and test — which maximises power
per stratum. The consequence, accepted deliberately, is that there is no single project-wide
false-discovery number: significance means "within this table". To keep that honest, the
dissertation's **primary confirmatory claim** is defined narrowly — the Up, Down and
Composite DTP scores against 3-year OS and RFS in the main CRC cohorts — and everything built
on top of the subtype calls (per-subtype survival cohorts, the score-across-subtype
Kruskal-Wallis, and the adjusted and interaction Cox models) is reported as **exploratory /
hypothesis-generating** rather than folded into the confirmatory error rate. This ties the
correction to the claim structure instead of pooling heterogeneous tests.

The molecular subtyping is now integrated (CMS + PDS; CRIS dropped) and its former blockers
are resolved. Both cohorts are classified as already-log2 expression with the CMS
classifier's own RNA-seq transform switched off, which removes the double-transform risk and,
because the same centering is applied to each cohort, removes the centering asymmetry as well.
This was verified directly: on TCGA-COAD the correct setting versus the old double-transform
changes 21 of 481 CMS calls, and the correct setting yields the expected CMS2-dominant
distribution. The earlier concern about the PDS classifier's input scale is resolved
empirically: `PDSpredict` produces a clean, non-degenerate spread across PDS1–3 in both
cohorts (roughly 20% Mixed and a balanced remainder in each), which is the signature of a
correctly scaled input — a mis-scaled matrix collapses the calls onto one class. The "Mixed"
label uses the `PDSpredict` package default confidence threshold of 0.6 (a sample is assigned
to a pure subtype only when that subtype's posterior probability exceeds 0.6), retained rather
than tuned, and the package is pinned to a specific commit for reproducibility.

PDS is kept as an exploratory modifier because it earns its place empirically rather than by
assumption: the score×PDS interaction on OS in GSE39582 is significant for all three core
scores and, for the Down score, is captured by PDS (FDR ≈ 0.01) where CMS entirely misses it
(FDR ≈ 0.86) — a non-redundant signal. The usual caveats apply and are stated as such: it is
single-cohort, does not replicate in the smaller TCGA cohort, and lives in the exploratory
family, so it is a hypothesis to pursue rather than a confirmed effect.
