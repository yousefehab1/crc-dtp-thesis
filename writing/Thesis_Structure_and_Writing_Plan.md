# Thesis Structure and Writing Plan

**Module:** Research Project and Dissertation (SCM8053), MSc Bioinformatics and Computational Genomics, Queen's University Belfast
**Working title:** A transcriptomic drug-tolerant persister programme as a prognostic and biological signal in colorectal cancer: evidence from survival, pan-cancer, and metastasis analyses
**Target length:** approximately 15,000 words excluding references, appendices, tables and figure legends
**Submission:** electronic via Canvas, first half of September 2026 (indicative 10 September 2026, to be confirmed)
**Marking:** two internal examiners, nine graded criteria (see mapping below)

---

## 1. Guiding decisions for this plan

Three modules are presented as one connected investigation of a single hypothesis rather than three separate studies. The unifying thread is that a transcriptomic drug-tolerant persister (DTP) programme, if it is a genuine cancer-stress state rather than a colorectal-specific artefact, should produce three independent and directionally consistent predictions, one tested by each module:

1. In primary colorectal cancer the signature should track patient outcome (crc_survival).
2. The association should generalise across many cancer types (pancan_survival).
3. The programme should be enriched in the compartment where persister selection has occurred, namely metastasis (mets_de).

This framing is stated once in the Introduction, restated as the study design in Methods, carried through as three results strands, and drawn back together in the Discussion. Every section should point back to it so the thesis reads as one argument.

The thesis follows the standard structure named in the study guide: Abstract, Introduction, Materials and Methods, Results, Discussion, Conclusions and Future Directions, References, Appendices.

---

## 2. Word budget

The allocation below sums to roughly 15,000 words of main text. Shorter is not penalised, so treat these as ceilings, not quotas. The abstract sits on its own page and is not part of the running total.

| Section | Target words | Share of body |
|---|---|---|
| Abstract | 400 (250 to 500) | not counted in body |
| 1. Introduction | 3,800 | 25% |
| 2. Materials and Methods | 3,200 | 21% |
| 3. Results | 4,800 | 32% |
| 4. Discussion | 2,800 | 19% |
| 5. Conclusions and Future Directions | 900 | 6% |
| Body total | approx. 15,500 | 100% |

Results is the largest section because it carries three modules. If the total runs long, trim Introduction background and Methods description first, and protect Results and Discussion.

---

## 3. How sections map to the nine examiner criteria

The Appendix 1 assessment form grades nine items. Keeping them in view while writing keeps effort aligned with marks.

| Examiner criterion | Where it is earned |
|---|---|
| 1. Abstract | Abstract |
| 2. Introduction | Introduction, especially the problem statement and aims |
| 3. Materials and Methods | Materials and Methods (reproducibility is the test) |
| 4. Results: handling of results and data analysis | Results plus the rigour of the statistics reported |
| 5. Discussion | Discussion |
| 6. Conclusions | Conclusions and Future Directions |
| 7. Ideas for continuation or modification | Future Directions and the open methodological questions |
| 8. Quality and clarity of English | throughout; budget editing time |
| 9. Organisation and presentation | overall structure, figures, tables, referencing consistency |

---

## 4. Section-by-section outline

### Front matter

Title page with author full name, degrees held, dissertation title as agreed with supervisor, degree offered, faculty, and date. A declaration that the work is the student's own original work except where acknowledged, including any use of AI, with a clearly stated word count. These are required by the study guide and are quick to prepare, but should not be left to the final day.

### Abstract (approx. 400 words, one page)

Structured in four moves as the guide suggests: introduction (background and purpose), methods, results, conclusions. Written last, after Results and Discussion are stable, because it must reflect the final findings. It should convey what was done, why, and what the findings imply for the field, in language a non-specialist scientist can follow.

### 1. Introduction (approx. 3,800 words)

This is where the January literature review is reused and extended. It should read as an almost complete introduction already, tightened and updated.

1.1 Colorectal cancer burden, treatment, and the problem of relapse. Establish why relapse after chemotherapy is the clinical problem the thesis speaks to.

1.2 The drug-tolerant persister state. Define the slow-cycling, reversible phenotype that a subpopulation of tumour cells adopts under cytotoxic pressure such as 5-FU, and its candidacy as a substrate for relapse. Review the persister literature across cancer types.

1.3 Transcriptomic signatures and single-sample scoring. Introduce signature-based approaches and why a rank-based, single-sample method is suitable. Set up ssGSEA conceptually without method detail.

1.4 Prognostic signatures in colorectal cancer and molecular subtypes. Review consensus molecular subtypes (CMS) and pathway-derived subtypes (PDS), microsatellite instability, and stage, since these are later used as confounders and effect modifiers.

1.5 Gap, hypothesis, and aims. State the gap the thesis addresses, the single hypothesis, and the three predictions with their corresponding aims. This subsection is what examiner criterion 2 rewards most, so it must be explicit and well justified.

### 2. Materials and Methods (approx. 3,200 words)

The test here is reproducibility by an independent scientist. Written mostly in the past tense and passive where natural. A shared-core-then-modules structure mirrors the actual pipeline and avoids repetition.

2.1 Data sources. GSE39582 (Affymetrix microarray), TCGA-COAD (RNA-seq, STAR), the wider TCGA atlas for pan-cancer, and GSE50760 (paired Normal, Primary, Metastasis RNA-seq). One paragraph plus the data-source table.

2.2 Signature definition. The master signature panel, the Up, Down and Composite (Up minus Down) construction, and the coverage gate that requires a minimum fraction of signature genes present before scoring.

2.3 Expression processing. Per-platform normalisation to one convention: log2-scale, gene-symbol-indexed matrices. TCGA uses tpm_unstrand transformed to log2(TPM+1); duplicates collapsed in log space by averaging for TCGA and array data, summed for the additive sub-features in the Mets FPKM data. State the reasoning briefly, since a choice with a rationale reads as deliberate.

2.4 Signature scoring. ssGSEA as a within-sample rank statistic, the score-naming convention, and why monotonic transforms leave scores unchanged.

2.5 Clinical endpoints. TCGA Clinical Data Resource, the shared PFI recurrence endpoint, censoring-aware landmark derivation at 36 months, and the Dead_3yr / Alive_3yr and Recurred / Recurrence-Free labels.

2.6 Statistical analysis. Mann-Whitney with rank-biserial effect size, Kaplan-Meier with log-rank for display, Cox proportional hazards on the continuous score for inference, the median split used for visualisation only, the gating thresholds, and Benjamini-Hochberg correction. Note that the correction family is a considered choice, cross-referenced to the open questions.

2.7 Molecular subtyping and adjustment. CMS and PDS calling from the same expression matrices, the log2 input handling that avoids a double transform, per-subtype survival cohorts, the score-across-subtype Kruskal-Wallis comparison, and the two Cox families: adjusted Cox for confounding and interaction Cox for effect modification, both gated by an events-per-variable rule.

2.8 Pan-cancer batch correction. Cancer type as a genuine batch for a cross-tissue score, removeBatchEffect applied to the aggregate only, per-cohort associations left uncorrected.

2.9 Metastasis differential expression. limma paired design with patient as a blocking factor on log2(FPKM+1), the liver-contamination covariate and its negative-control validation, the purity-corrected matrix, GSEA on the ranked statistic, PCA, and paired single-sample ROC.

2.10 Reproducibility. Global seed, timestamped run directories, caching, and captured session information.

### 3. Results (approx. 4,800 words) — write this first

Organised so that the colorectal strand builds from a primary univariable result to a subgroup breakdown to an independence test, before the two generalisation strands (pan-cancer and metastasis) and a closing synthesis. Each subsection states the question, the finding, and points to the figure or table. Interpretation is kept light here and saved for the Discussion, but the qualified reading of each result belongs here. The boundary between subsections is the type of analysis: 3.2 and 3.3 are univariable, 3.4 is multivariable.

3.1 Cohort and data overview. What entered each analysis and why, including the pan-cancer inclusion and exclusion accounting. A short orientation before the findings.

3.2 Primary association: the signature and colorectal cancer outcome (univariable). Association of the DTP score with overall survival and recurrence in GSE39582 and TCGA-COAD, in the whole cohort and the adjuvant-treated subset only. Report Wilcoxon screens with rank-biserial effect size, Kaplan-Meier separation, and continuous-score Cox hazard ratios with confidence intervals and concordance. State cross-platform consistency plainly, and be explicit that whole-cohort overall survival is borderline while recurrence and the adjuvant-treated subset are the cleaner signals. Figures 1 and 2 live here, supported by a results table carrying the exact effect sizes, hazard ratios, intervals and FDR p-values.

3.3 Signal across molecular and clinical subgroups (univariable). The same DTP-to-outcome test run within CMS subtypes, PDS subtypes, stage groups, and MSI or microsatellite-stable status, plus the score-across-subtype question (does the DTP score itself differ across CMS or PDS, Kruskal-Wallis with epsilon-squared). Presented as one compact table showing the events-per-variable gate result per subgroup and at most one summary figure. The subsection states up front that most subgroups are underpowered, so a non-significant or not-testable result is not evidence of absence. Framed this way the section supports the specificity of the primary signal rather than diluting it.

3.4 Independence from clinicopathology (multivariable). Adjusted Cox results showing whether the score's hazard ratio survives adjustment for stage, MSI, CMS and PDS, the change-in-estimate on the same patients, and the likelihood-ratio test of added value beyond clinicopathology. Interaction Cox results and the per-subgroup forest. This is the "the signal is independent of, and not merely explained by, the standard taxonomy" section.

3.5 Generalisation across cancer types (pancan_survival). The aggregate batch-corrected pan-cancer association and the per-cohort forest on uncorrected scores. Whether the direction is consistent with a shared stress programme.

3.6 Enrichment in metastasis (mets_de). Differential expression Metastasis versus Primary before and after purity adjustment, the contamination proxy validation, GSEA enrichment, PCA separation by tissue, and paired ROC for signature separability.

3.7 Synthesis of results. One short paragraph drawing the three strands together at the level of evidence, without yet moving to interpretation.

### 4. Discussion (approx. 2,800 words)

4.1 Principal findings. Restate what the three strands collectively show about the DTP programme in colorectal cancer, in relation to the hypothesis.

4.2 Interpretation in the context of the literature. Place the findings against persister biology, prognostic signatures, and molecular subtypes. Critically evaluate significance and potential importance, which is what examiner criterion 5 rewards.

4.3 Strengths. Cross-platform replication, identical processing of the shared COAD cohort, biologically motivated stratification, explicit confounding and effect-modification analysis, and direct handling of the metastasis contamination confound.

4.4 Limitations. Landmark analysis trades sample size for correct classification; proportional-hazards assumption; sample-relative ssGSEA scores depend on consistent normalisation; interaction tests are low-powered and a null is not evidence of uniformity; adjustment removes only measured confounders; multiple-testing correction reduces but does not eliminate false positives and its family is not yet finalised; the PDS input-scale check remains open.

4.5 Open methodological questions. The false-discovery correction family (narrow per-stratum versus pooled project-wide) and the PDS classifier input scale, presented as considered, flagged decisions rather than oversights.

### 5. Conclusions and Future Directions (approx. 900 words)

A concise statement of what can and cannot be concluded from the work, followed by concrete next steps that build on it: resolving the correction family against the number of independent primary tests, validating the PDS calls, external validation in an independent treated colorectal cohort, and experimental follow-up of the persister hypothesis. This section carries examiner criteria 6 and 7, so the future directions should be specific and grounded in the results, not generic.

### References

Alphabetical, full details, in the Watson and Crick style specified by the guide: one or two authors named in text, three or more as first author et al. Every reference cited must have been read and every reference listed must be cited. Choose a reference manager early and keep the style consistent.

### Appendices

Supporting material that would break the flow of the main text: full result tables (Adjusted_Cox_Summary, Interaction_Cox_Summary, per-cohort pan-cancer results), supplementary figures, the signature panel, and session information for reproducibility. Appendices are excluded from the word count.

---

## 5. Suggested writing order

Results first, as chosen, because the analysis is complete and documented and writing around the existing outputs is the fastest way to build momentum. A workable order from there:

1. Results (all subsections in the three-strand order above).
2. Materials and Methods (already largely specified by the pipeline; write while Results is fresh so the two stay consistent).
3. Introduction (extend the January literature review to a full introduction, ending on the hypothesis and aims).
4. Discussion (needs stable Results and a complete Introduction to reference).
5. Conclusions and Future Directions.
6. Abstract (last, so it reflects the final findings).
7. Front matter, reference tidy, and a full English and consistency edit for criteria 8 and 9.

---

## 6. Open items to confirm with supervisor

The study guide says the final format should be agreed with the primary supervisor, so a few points are worth settling before drafting far:

- Confirmation of the submission date and the exact word-count expectation.
- Agreement that all three modules belong in the main narrative rather than one being moved to an appendix.
- The referencing style, since the guide allows variations with supervisor approval.
- Structure decided: 3.2 is the primary univariable association (whole cohort and adjuvant-treated only); 3.3 is the univariable subgroup breakdown (CMS, PDS, stage, MSI or MSS, plus score-across-subtype); 3.4 is the multivariable independence analysis (adjusted and interaction Cox). Confirm with supervisor whether the full 3.3 subgroup table stays in the main text or moves to an appendix given its length.
- Resolution, or at least a stated position, on the two open methodological questions before the Discussion is finalised.
