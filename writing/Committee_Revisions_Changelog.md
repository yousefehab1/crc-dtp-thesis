# Committee revisions — change log

This document lists every change made to the thesis in response to the committee feedback, grouped by reviewer, with the location and the exact text added or revised. It is a tracking aid; nothing here is part of the thesis itself.

Two items still require input from you: the **code repository URL** (placeholder left in the Data and Code Availability section) and nothing else outstanding.

---

## Dr. Whitfield — Statistics

**#1 Pooled-FDR sensitivity across all subgroup tests** — Section 3.3 (added) and Section 4.8 (added reference).
Added to 3.3: "when the false-discovery correction is pooled across all 159 testable subgroup tests, rather than applied within each stratum, only one association survives (stage III overall survival in GSE39582, q = 0.015), against fourteen under the per-stratum scheme. The subgroup analysis is therefore hypothesis-generating and does not, on its own, establish effect in any particular stratum."
Added to 4.8: "pooling the correction across all 159 testable subgroup tests, the most conservative option, leaves only the strongest subgroup association standing (Section 3.3)."

**#2 Justify removeBatchEffect + Wilcoxon vs stratified Cox** — Section 2.8 (added paragraph).
"a stratified Cox model with cohort as a stratum is the standard alternative... the per-cohort Cox forest (Figure 10) is already the stratified, within-cohort view... while the batch-corrected-score analysis exists only to give a single, interpretable cross-tissue summary... a deliberate presentational choice for the aggregate, sitting alongside, not in place of, the stratified per-cohort models."

**#3 Statistically significant vs clinically useful (concordance 0.58–0.65)** — Section 3.2 (added sentence).
"a concordance of 0.58 to 0.65 is reliably above the 0.50 of chance yet well short of the 0.70 or higher generally expected of a stand-alone clinical prognostic tool, so the score is a real but modest discriminator rather than a ready-made classifier."

**#4 Landmark pre-specified** — Section 2.5 (revised).
Now: "Outcomes were analysed at a pre-specified 36-month (three-year) landmark, fixed as a design parameter before any outcome analysis. Three years was chosen because the majority of colorectal cancer recurrences occur within this window, so the landmark captures most events while retaining sufficient follow-up."

---

## Dr. Osei — Biology

**#5 How the Up/Down signature was derived** — Section 2.2 (added).
"The signature was derived from HCT-116 colorectal cancer cells maintained under 5-fluorouracil for 14 days to select a surviving drug-tolerant persister population; the up and down components comprise the genes consistently up- and down-regulated in these persister cells relative to untreated parental cells, with the exact selection thresholds to be reported in the laboratory's forthcoming publication."

**#6 Coherent biological story for the Down score** — Section 4.4 (added paragraph).
Based on the GO over-representation analysis you ran (cytoplasmic translation, adjusted p = 1e-28): "the down set is best understood not as a second persister axis but as its opposite: a high down score marks biosynthetically active, proliferative, non-persister cells. That single interpretation ties together the otherwise disconnected places the down score appeared" — metastasis (translational demand of outgrowth), proliferative LUAD/CESC, null in primary CRC, and why Composite (Up − Down) sharpens the persister signal.

**#7 Why deconvolution was not attempted** — Section 4.4 (added).
"Computational deconvolution of the bulk profiles was considered as a partial remedy but not pursued, because without a reference expression profile for the persister state itself, deconvolution can estimate broad cell-type fractions but cannot isolate a rare within-epithelial persister subpopulation; it is therefore left, with single-cell profiling, for future work."

---

## Dr. Marchetti — Reproducibility

**#8 Delete word-count placeholder** — Title page (removed).
Deleted: "Placeholder for tracking; the declaration and final count will be completed on the submission form." The line now reads only: "Word count (body, Introduction to Conclusions): approximately 13,705 words."

**#9 Data and code availability statement** — new section added after Conclusions.
States the public accessions (GSE39582, GSE50760, TCGA/GDC, TCGA-CDR), a repository-URL placeholder ("[repository URL to be added]"), and that the DTP signature is withheld pending the Atlasi Laboratory's primary publication.

**#10 Software-version table** — Section 2.10 (added Table 2.5).
R 4.5.3; GSVA 2.4.9; GSEABase 1.72.0; limma 3.66.0; clusterProfiler 4.18.4; TCGAbiolinks 2.38.0; CMScaller 0.99.2; PDSclassifier commit c89a19c; survival 3.8-9.

**#11 Pinned PDSclassifier version in-text** — Section 2.7 (added).
"the classifier was pinned to a fixed version for reproducibility (PDSclassifier, GitHub commit sidmall/PDSclassifier@c89a19c)."

---

## Dr. Halloran — Readability

**#12 Three-prediction verdict table** — Section 3.1 (added Table 3.0).
A three-row table: Prediction 1 (primary CRC) — Supported; Prediction 2 (pan-cancer) — Partial; Prediction 3 (metastasis) — Not supported at the bulk level.

**#13 Tighten 3.5 "remaining significant results" + closing takeaway** — Section 3.5 (revised).
Condensed, and now closes with: "The takeaway is one of concentration rather than spread... the programme is prognostic in particular tissue contexts rather than uniformly across the atlas."

---

## Prof. Adeyemi — Chair / Translational

**#14 Comparison to existing clinical tools** — Section 4.7 (added).
"the score was not compared head-to-head against the validated clinical prognostic assays already used in colorectal cancer, such as Oncotype DX Colon and ColoPrint (Salazar et al. 2011), or against stage and microsatellite-status nomograms. It is presented here as a mechanistic, biologically motivated signal that is complementary to, rather than a replacement for, those tools."

**#15 Clinic-ready assay requirements** — Section 5.2 (added).
"A clinic-ready version of the assay would not require whole-transcriptome profiling: the signature could be locked onto a targeted panel, for example a NanoString or quantitative PCR readout of its member genes, with a fixed, pre-registered scoring threshold, which would make it deployable on routine formalin-fixed clinical material."

**#16 Minimum detectable effect in the metastasis limitation** — Section 4.4 (revised).
Now states the number rather than "underpowered": "for a paired design with about eighteen patients, the minimum standardised effect detectable at 80% power and a two-sided 5% level is approximately 0.66, a medium-to-large shift, so only a substantial, tissue-wide difference could have been detected."

---

## Reference list

Added: Salazar R, Roepman P, Capella G, Moreno V, Simon I, Dreezen C, et al. (2011) Gene expression signature to improve prognosis prediction of stage II and III colorectal cancer. Journal of Clinical Oncology 29:17–24. (For the ColoPrint comparison, #14.)
