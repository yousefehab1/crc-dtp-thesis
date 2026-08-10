# Appendix: supplementary result tables

The tables below are referenced from the Results chapter and are excluded from the word count. To keep the appendix readable, the tables are not printed here; each is provided as a downloadable supplementary data file and linked below, so the full statistics behind every figure are one click away without breaking the flow of the text. All values come from run CRC_DTP_20260806_1651.

## Appendix Table A1. Univariable DTP-to-outcome association within molecular and clinical subgroups

Referenced from Section 3.3 (Figure 6). Downloadable supplementary file: [AppendixTableA1_subgroup_survival.csv](tables/AppendixTableA1_subgroup_survival.csv).

One row per Dataset, subgroup Category (CMS, PDS, Stage, MSI), Subgroup level, Score (DTP up, down, composite) and Endpoint (overall survival, recurrence). Columns: sample size (N), event count (N_events), testability flag (Is_Testable), rank-biserial effect size (Effect_r) with its Wilcoxon FDR p-value, hazard ratio per one standard deviation of the score (HR_perSD) with 95% confidence interval (HR_lower, HR_upper), concordance index (C_index), and the Benjamini-Hochberg corrected Cox and log-rank p-values (Cox_FDR_P, KM_logrank_FDR_P). The table holds 168 rows, of which 159 are testable; the 9 not-testable rows are the stage I strata below the events-per-variable gate.

## Appendix Table A2. Adjusted Cox models: DTP score added value beyond clinicopathology

Referenced from Section 3.4 (Figure 8). Downloadable supplementary file: [AppendixTableA2_adjusted_cox.csv](tables/AppendixTableA2_adjusted_cox.csv).

One row per Dataset, Endpoint, Score and adjustment Model (clinicopathology, CMS-adjusted, PDS-adjusted). Columns: adjusted hazard ratio (HR) with 95% confidence interval, raw p-value, concordance index (C_index), likelihood-ratio p-value for the score's added fit beyond the covariates (LRT_score_P), proportional-hazards test p-value (PH_P), the unadjusted hazard ratio (Unadj_HR), the change in the score's log hazard ratio after adjustment as an absolute value and a percentage (Delta_logHR, Delta_logHR_pct), sample size and events (N, N_events), testability flag, and the Benjamini-Hochberg corrected p-value with its significance flag.

## Appendix Table A3. Interaction Cox models: effect modification of the DTP score

Referenced from Section 3.4 (Figure 9). Downloadable supplementary file: [AppendixTableA3_interaction_cox.csv](tables/AppendixTableA3_interaction_cox.csv).

One row per Dataset, Endpoint, Score and Modifier (CMS, PDS, stage, microsatellite status). Columns: interaction likelihood-ratio raw p-value (Raw_P), the number of modifier levels (K), sample size and events (N, N_events), testability flag, and the Benjamini-Hochberg corrected p-value with its significance flag. The corrected test compares a score-by-modifier Cox model against the additive score-plus-modifier model.

## Appendix Table A4. Pan-cancer per-cohort and pooled DTP-score survival association

Referenced from Section 3.5 (Figures 10 and 11). Downloadable supplementary file: [AppendixTableA4_pancan_survival.csv](tables/AppendixTableA4_pancan_survival.csv).

One row per cohort per score per endpoint across the 25 solid-tumour TCGA cohorts, plus the batch-corrected pooled aggregate rows. Columns: Project, Score, Endpoint, sample size and events (N, N_events), testability flag, per-SD hazard ratio with 95 percent confidence interval (HR_perSD, HR_lower, HR_upper), concordance index (C_index), and the Benjamini-Hochberg corrected Cox p-value with its significance flag. Accompanying appendix figures for Section 3.5 are the batch-corrected pooled hazard-ratio forest (`Fig_Pancan_aggregate_forest`) and the pooled Kaplan-Meier of the batch-corrected score, both of which show the near-chance separation of the pooled association.

## Appendix Table A5. Univariable association of the DTP score with colorectal cancer outcome

Referenced from Section 3.2 (Figures 4 and 5). One row per Cohort, Subset (whole cohort, adjuvant-treated), Score (DTP up, down, composite) and Endpoint (overall survival, recurrence). Columns: sample size and events, the Mann-Whitney rank-biserial effect size with its Wilcoxon FDR, the continuous-score Cox hazard ratio per one standard deviation of the score with 95% confidence interval, the concordance index, and the Benjamini-Hochberg corrected Cox and log-rank p-values. Hazard ratios are per standard deviation so they are comparable across rows; the adjuvant-treated subset of TCGA-COAD (n = 23) failed the events-per-variable gate and is included for completeness only. Downloadable supplementary file: [AppendixTableA5_primary_crc_univariable.csv](tables/AppendixTableA5_primary_crc_univariable.csv). Source: CRC_Statistical_Summary.csv, run CRC_DTP_20260806_1651.

## Appendix Table A6. Metastasis differential-expression and enrichment statistics

Referenced from Section 3.6 (Figures 12 and 13). These tables give the full statistics behind the metastasis strand: the per-signature GSEA normalised enrichment scores and FDR for each contrast, the liver-score validation, and the single-sample separability diagnostics. Because the metastasis analysis produces several distinct result files, each is provided as a downloadable supplementary file rather than printed here.

Normal-versus-primary GSEA (unadjusted): [AppendixTableA6a_mets_gsea_NvP.csv](tables/AppendixTableA6a_mets_gsea_NvP.csv). Primary-versus-metastasis GSEA on the liver-adjusted ranking: [AppendixTableA6b_mets_gsea_PvM_adjusted.csv](tables/AppendixTableA6b_mets_gsea_PvM_adjusted.csv). Primary-versus-metastasis normalised enrichment scores before and after liver adjustment, the values behind Figure 12B: [AppendixTableA6c_mets_gsea_PvM_adj_vs_unadj.csv](tables/AppendixTableA6c_mets_gsea_PvM_adj_vs_unadj.csv). Per-patient liver score across the matched normal, primary and metastasis samples, the values behind Figure 12A: [AppendixTableA6d_mets_liver_scores.csv](tables/AppendixTableA6d_mets_liver_scores.csv). Paired single-sample ROC areas under the curve for signature separability: [AppendixTableA6e_mets_roc_auc.csv](tables/AppendixTableA6e_mets_roc_auc.csv). Source: run CRC_DTP_20260808_1642, mets_de module.
