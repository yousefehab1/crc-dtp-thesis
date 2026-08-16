## 3.1 Cohort and data overview

The Results test one hypothesis through three predictions, each on its own data, over a shared signature and scoring core (Figure 3). The drug-tolerant persister signature, its up, down and composite scores together with six reference signatures, was scored in every dataset by ssGSEA single-sample scoring, so the same programme is measured the same way throughout. The three strands then ask whether that programme tracks outcome in primary colorectal cancer (Sections 3.2 to 3.4), generalises across cancer types (Section 3.5), and is enriched in metastasis (Section 3.6), before Section 3.7 draws them together.

The colorectal strand used two independent cohorts, GSE39582 (Affymetrix microarray, about 585 patients) and TCGA-COAD (RNA-seq, 455 patients), each with overall survival and recurrence endpoints, analysed in the whole cohort, the adjuvant-treated subset, and within CMS, PDS, stage and microsatellite subgroups. The pan-cancer strand drew on the TCGA Clinical Data Resource: of 33 solid-tumour cohorts, two were removed as non-solid where the persister notion does not apply (TCGA-LAML and TCGA-DLBC), five fell below the 100-patient threshold (ACC, MESO, UVM, UCS and CHOL), and one dropped below it after the clinical join (TCGA-KICH, 66 patients), leaving 25 cohorts and about 8,900 patients. The metastasis strand used GSE50760, paired RNA-seq of matched normal mucosa, primary tumour and liver metastasis from roughly 18 patients. The datasets, their sizes and the analyses each entered are laid out in Figure 3. The overall verdict on the three predictions is summarised in Table 3.0 and developed in the sections that follow.

Table 3.0. Verdict on the three predictions of the hypothesis.

| Prediction | Test | Verdict |
|---|---|---|
| 1. Tracks outcome in primary colorectal cancer | Univariable and multivariable survival in GSE39582 and TCGA-COAD (Sections 3.2 to 3.4) | Supported |
| 2. Generalises across cancer types | Per-cohort and pooled survival across 25 TCGA cohorts (Section 3.5) | Partial |
| 3. Enriched in the metastatic compartment | Liver-adjusted differential expression and GSEA in GSE50760 (Section 3.6) | Not supported at the bulk level |

With the cohorts, scores and endpoints defined, the first strand asks the most direct question, whether the persister score is associated with outcome in primary colorectal cancer.

![Figure 3](figures/Fig_study_design.png)

**Figure 3. Study design and analysis overview.** The persister signature and six reference signatures were scored by ssGSEA in every dataset and tested against the three predictions of the hypothesis: association with outcome in primary colorectal cancer, generalisation across solid-tumour types, and enrichment in matched metastasis. The datasets, sample sizes, analyses and corresponding Results sections for each strand are shown in the panels.
