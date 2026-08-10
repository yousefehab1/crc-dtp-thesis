# Section 4 — Discussion: writing plan

Purpose: weigh what the three strands collectively say about the hypothesis, honestly and in the context of the literature. This is where the mixed result (supported in primary colorectal cancer, partial across tissues, null in bulk metastasis) is turned into a coherent argument rather than left as three separate verdicts. Target about 2,800 words. House style: prose, no em dashes, no drama, past tense for results and present tense for interpretation.

## Design principle: one sequential conclusion per subsection

Following the supervisor's note on the Results, every subsection ends with a single plain conclusion sentence, and that sentence is the premise the next subsection opens on. The chain is the argument. The conclusions, read in order, should stand alone as a one-paragraph summary of the Discussion:

1. The persister up programme is a real and specific prognostic signal in treated primary colorectal cancer.
2. Because it is real there, the question becomes whether it generalises; across tissues it does so only partly, so it is prognostic in particular tumours rather than a universal stress marker.
3. Because generalisation is context-dependent, the natural next test is the compartment where persisters are selected; the bulk metastasis result cannot resolve that prediction, because bulk tissue dilutes a fraction-of-cells state, so its null is uninformative rather than negative.
4. Taken together, the three predictions validate the programme as a genuine but context-limited prognostic state whose role in relapse and metastasis remains open.
5. That reading is only as strong as the design allows, so the strengths that support it and the limitations that soften the nulls must be set out, followed by the two methodological choices still to be settled.

## Subsection-by-subsection

### 4.1 Principal findings

Cover: restate the hypothesis and its three predictions in one short paragraph, then give the verdict on each without re-reporting numbers, so the reader has the whole shape before the interpretation begins. Supported for primary colorectal outcome; partial across cancer types; null in bulk metastasis after liver adjustment. Note once that the persister up component, not the down component, carried the colorectal signal, and that the composite behaved like the up score.

Simple conclusion (hands to 4.2): the clearest and most secure of the three results is the primary colorectal association, so the interpretation begins there.

### 4.2 Prediction 1 — a real and specific prognostic signal in colorectal cancer

Cover:
- Why the primary result is trustworthy: it replicated across two cohorts and two platforms, was clearest for recurrence and in the adjuvant-treated subset, and survived adjustment for stage, microsatellite status, CMS and PDS (Section 3.4), so it is not a restatement of the standard taxonomy.
- The specificity argument from the reference panel: the fetal intestinal programme tracked outcome alongside the up score while the stemness (CBC), revival (revSC), regenerative (RSC) and inflammatory (IBD) programmes did not, and MYC ran the opposite way. So the signal sits on the persister and fetal-regenerative axis, not on generic stemness or proliferation. This is the single most important specificity point in the thesis.
- The adjuvant-treated subset carrying the strongest effect fits the persister interpretation directly: cytotoxic pressure is where a drug-tolerant state should matter most (link to Rehman et al. 2021; Sharma et al. 2010).
- Literature link, to be developed rather than mentioned in passing: the persister and diapause state is defined by slow cycling (Sharma et al. 2010; Rehman et al. 2021), and an independently derived, patient-based classification has recently identified a slow-cycling colorectal subtype, PDS3, as carrying the worst prognosis in locally advanced disease (Malla et al. 2024). Develop the parallel in three steps. First, the existence of PDS3 shows that slow-cycling biology is genuinely and independently prognostic in colorectal cancer, which lends external plausibility to a slow-cycling persister score also being prognostic. Second, the analysis here supports an overlap directly: the persister score varied most across the PDS axis (Section 3.3) and PDS absorbed more of the overall-survival signal than CMS or clinicopathology (Section 3.4), so the score is partly capturing the pathway-level slow-cycling and regenerative biology that PDS also measures. Third, hold the honest boundary: the mapping is not one-to-one, because the DTP subgroup signal was clearest in the regenerative PDS2 stratum rather than in PDS3 itself, so the relationship is one of shared slow-cycling and regenerative biology rather than identity between the persister score and any single subtype. Flag that this argument rests on one external study (Malla et al. 2024) and should be framed as convergent evidence, not proof.

Simple conclusion (hands to 4.3): a signal that is real, specific and strongest under chemotherapy in colorectal cancer raises the question of whether it marks a state shared across cancers or one peculiar to this tissue.

### 4.3 Prediction 2 — prognostic in some tissues, not a universal stress marker

Cover:
- The honest reading of the pan-cancer result: the direction was broadly consistent but significance concentrated in a few tissues, the pooled effect was small, and discrimination was near chance when pooled.
- Where it was strong is the interesting part: low-grade glioma and pancreatic adenocarcinoma showed the up and composite scores significant at both endpoints with genuine discrimination (concordance up to 0.73), exceeding colorectal itself. Attempt a mechanistic explanation, hedged, for why these two tissues in particular:
  - Both are tumours in which non-genetic plasticity and slow-cycling drug-tolerance are recognised drivers rather than incidental features, which is the biology a persister programme should mark. In pancreatic adenocarcinoma, therapy resistance runs through epithelial-to-mesenchymal and fetal or regenerative (acinar-to-ductal) reprogramming and documented drug-tolerant persister states under chemotherapy, the same fetal-regenerative axis the up score tracked in colorectal cancer (Section 3.2). In low-grade glioma, prognosis is dominated by dedifferentiation and stem-like plasticity, so a slow-cycling stress programme plausibly captures the more plasticity-prone, aggressive tumours.
  - State the shared thread explicitly: in both tissues the persister and fetal-regenerative axis is a known route to aggressiveness, so a prognostic signal there is consistent with the programme marking a genuine stress state rather than tissue-specific noise.
  - Give the honest methodological caveat alongside the mechanism: low-grade glioma carries an unusually wide prognostic spread (grade and IDH status), so any biologically meaningful axis discriminates well there, which inflates concordance and means the high value is not, on its own, evidence of a persister-specific effect.
  - Citations for this paragraph must be gathered and verified before writing (pancreatic drug-tolerant persister and EMT/fetal reprogramming; glioma plasticity and dedifferentiation); do not draft the mechanism until those are pinned.
- The reversed kidney (KIRP) direction and the scattered down-score and composite hits show the programme is not uniformly oriented across tissues, which argues against a single universal stress axis and for context-dependent prognostic value.
- Caution: pan-cancer cohorts are untreated-heterogeneous and the endpoint is progression-free interval; the batch correction and its limits (Section 2.8) should be acknowledged briefly here rather than in Limitations.

Simple conclusion (hands to 4.4): if the programme is prognostic only in particular tissue contexts, the sharpest remaining test is not another tissue but the compartment where persister selection itself is expected, the metastasis.

### 4.4 Prediction 3 — the bulk metastasis null is uninformative, not negative

Cover:
- State the result plainly: after removing the liver-contamination confound, the up programme showed no enrichment in the metastasis in either direction, so this dataset does not support the metastasis prediction at the bulk level.
- The central interpretive move, set up in Section 1.2 and 1.3: a persister state is occupied by only a fraction of cells at any time, and a bulk score averages over all cells, so a small paired bulk cohort is exactly the setting in which a real but sparse signal would be invisible. The null is therefore weak evidence of absence, not evidence of a true absence.
- Coverage is no longer a caveat: after the identifier-harmonisation fix, the up set is scored on 74 percent of its genes in GSE50760, above the 70 percent gate, so the metastasis null is a genuine bulk-level finding rather than a coverage artefact. State this once so the null rests solely on the dilution argument, and connect to Methods 2.2 (coverage gate) and to Limitations (4.7).
- What did move is informative: the down programme, the MYC module, and the revival and regenerative stem programmes were up in the metastasis, while the fetal programme only trended. So the metastatic compartment in bulk engages proliferation and regenerative or revival stemness rather than the up-regulated persister programme, which is a coherent biological picture and a lead for single-cell follow-up.
- Credit the contamination control as a methodological strength (the liver-adjusted ranking changed the result and removed a spurious IBD primary-enrichment), foreshadowing 4.6.

Simple conclusion (hands to 4.5): with one prediction supported, one partial and one unresolved for a defensible reason, the three strands can be read together as a single, qualified verdict on the hypothesis.

### 4.5 What the three predictions say together

Cover:
- Draw the strands into one statement: the programme behaves as a genuine, specific, and context-limited prognostic state, strongest where cytotoxic selection is strongest, consistent with a persister interpretation but not shown to be a metastasis-seeding programme in bulk.
- The recurring specificity thread: across all three strands the up and fetal-regenerative axis carried the prognostic weight while generic stemness and, in metastasis, the down and revival programmes behaved differently, so the panel design did its job of separating the persister signal from its neighbours.
- What would distinguish a real state from an artefact, and where this work lands on that spectrum: it clears the primary-outcome and independence bars decisively, clears the cross-tissue bar partially, and leaves the metastasis bar untested rather than failed.

Simple conclusion (hands to 4.6): this reading is only as strong as the design that produced it, so its support and its soft spots must be stated directly.

### 4.6 Strengths

Cover: cross-platform and cross-cohort replication of the primary result; identical processing of the shared TCGA-COAD cohort across arms; a biologically motivated reference panel that enables the specificity argument; explicit separation of confounding and effect-modification analyses; and the direct handling of the metastasis liver-contamination confound with a negative-control validation.

Simple conclusion (hands to 4.7): these strengths make the positive result credible, but they do not remove the limitations that make the two weaker results provisional.

### 4.7 Limitations

Cover: the landmark analysis trades sample size for correct classification; the proportional-hazards assumption; sample-relative ssGSEA scores depend on consistent normalisation; interaction and metastasis tests are low-powered, so their nulls mark limited power rather than proven uniformity or absence; adjustment removes only measured confounders, so residual confounding remains possible; the bulk-versus-cell-state dilution problem is the key limit on the metastasis prediction; and multiple-testing correction reduces but does not eliminate false positives.

Simple conclusion (hands to 4.8): two of these limitations are not fixed facts of the data but decisions still open in the analysis, and they should be named as such.

### 4.8 Open methodological questions

Cover: the false-discovery correction family, narrow per-stratum versus pooled project-wide, which the sensitivity of the borderline colorectal results to the family makes a live question (Section 3.2, the recurrence Wilcoxon shift); and the PDS classifier input scale. Present both as considered, flagged decisions rather than oversights, and state the position taken. The position on the FDR family is now fixed in Methods 2.6 and should be restated here, not re-argued: narrow per-stratum correction is used; raw and corrected p-values are both reported; the confirmatory primary-colorectal conclusion survives conservative correction and does not depend on the family; and borderline results are treated as borderline, with the TCGA-COAD recurrence shift from 0.032 to 0.064 (no change in data, only in the family) given as the worked illustration. The Discussion adds only the interpretive point: this is why the headline claim rests on the confirmatory result, not on any single borderline test.

Simple conclusion (hands to Conclusions): with the evidence weighed, its strengths and soft spots set out, and the open choices named, what can and cannot be concluded from the work can now be stated plainly.

## References to add beyond those already in the Introduction

- A persister or minimal-residual-disease review for the relapse-seeding argument in 4.4 (candidate already cited: Marine et al. 2020; consider one focused on residual disease).
- Tissue-biology references for the mechanistic interpretation of the glioma and pancreatic signals in 4.3, now required rather than optional: a pancreatic drug-tolerant persister or EMT/fetal-reprogramming reference, and a glioma plasticity or dedifferentiation reference. Verify each by direct reading before drafting; do not assert a mechanism on an unverified citation.
- Batch-correction and removeBatchEffect caveats if a methods citation is wanted in 4.3 (Ritchie et al. 2015, already in Methods).
- Keep all interpretation anchored to references already gathered; do not introduce claims that need new unverified citations.

## Open items to confirm before writing

1. Resolved: 4.3 will attempt a hedged mechanistic explanation of the glioma and pancreatic signals (persister and fetal-regenerative plasticity in both tissues), paired with the low-grade-glioma prognostic-spread caveat. This requires verified tissue-biology citations, which must be gathered before that paragraph is drafted.
2. Resolved: the FDR-family position is now fixed in Methods 2.6 and Table 2.3 (narrow per-stratum correction; raw and corrected p-values both reported; confirmatory conclusion robust to conservative correction; borderline results flagged and family-sensitive, with the TCGA-COAD recurrence 0.032-to-0.064 shift as the illustration). Section 4.8 restates this position rather than re-arguing it.
3. Resolved: the slow-cycling PDS3 link in 4.2 will be developed further, as a three-step convergent-evidence argument (PDS3 shows slow-cycling biology is prognostic; the score overlaps the PDS axis in this analysis; but the mapping is not one-to-one), while flagging that it rests on one external study (Malla et al. 2024).
