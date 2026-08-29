# INSPIRE Arthroplasty Anesthesia Reproducibility Repository

**Version 1.01**

Study: **Neuraxial versus general anesthesia and early postoperative resource use
after lower-limb arthroplasty**, a retrospective cohort study using the INSPIRE
perioperative registry.

This package reflects the manuscript AFTER external review and the full round of
analysis corrections. It supersedes all earlier packages.

## Headline finding (verified against source by the author)
Neuraxial anesthesia was associated with shorter hospital length of stay and lower
odds of any postoperative ICU admission, robust across BMI adjustment, patient-
clustered standard errors, IPTW, a broader orthopedic cohort, and a one-operation-
per-patient sensitivity analysis. A prespecified subgroup analysis found the
association concentrated in HIP arthroplasty; the knee subgroup showed a near-null
LOS association and an underpowered/inconclusive ICU estimate. The anesthesia-by-
procedure interaction was NOT significant (p=0.086), so the joint-specific pattern
is reported as hypothesis-generating, not established effect modification.

Key numbers (arthroplasty, n=2,540; general 269 / neuraxial 2,271):
- Primary LOS ratio 0.892 (0.827–0.963); ICU OR 0.494 (0.315–0.775)
- BMI-adjusted + patient-clustered: LOS 0.896 (0.829–0.969); ICU 0.478 (0.301–0.761)
- Hip: LOS 0.874 (0.785–0.974); ICU 0.406 (0.241–0.684)  [significant]
- Knee: LOS 0.961 (0.888–1.039, ns); ICU 1.382 (0.356–5.37, underpowered, 35 events)
- IPTW: LOS 0.900 (0.842–0.961); ICU 0.440 (0.272–0.710)  [BMI-adjusted point est. unchanged]
- All-orthopedic (BMI-adj, clustered): LOS 0.917 (0.870–0.968); ICU 0.664 (0.563–0.782)

## /manuscript
- `INSPIRE_manuscript_CURRENT.docx` — THE live manuscript. Contains: Option-B ICU
  definition (any postoperative ICU admission, matching the code); BMI + patient-
  clustering robustness; hip/knee subgroup Table + interaction caveat; updated
  limitations (calendar-year not recoverable from de-identified offsets; anesthesia
  coded mutually exclusive; knee ICU underpowered); trimmed self-citations
  (refs 16–17 removed; list ends at 15). Four tables, references complete 1–15.
- `INSPIRE_tables_figures.docx` — standalone tables/figures reference doc.

## /figures
- `Figure1_forest_subgroups_300dpi.png` — MAIN figure: primary, sensitivity, and
  hip/knee subgroup estimates for LOS and ICU (300 dpi). Knee-ICU arrow to 5.37;
  interaction p=0.086 noted.
- `Figure2_loveplot.png` — covariate balance before/after IPTW.
- `Figure3_strobe_flow.png` — STROBE cohort flow.
- `FigureS1_los_distribution.png`, `FigureS2_los_violin.png`, `FigureS3_icu_rate.png`
  — supplementary raw-data figures.

## /code
- `inspire_analysis_working.py` — the author's working analysis script (definitive
  source of the numbers). ICU defined as icuin_time.notna() = any admission.

## /archive_earlier_versions
Superseded manuscript drafts kept for provenance only:
- `INSPIRE_JFMK_manuscript.docx` (earliest, JFMK-framed)
- `INSPIRE_manuscript_ICU_fixed_optionB.docx` (ICU wording fixed)
- `INSPIRE_manuscript_final.docx` (self-citations trimmed, pre-subgroup)
Use CURRENT, not these.

## Status: what is DONE
- Analysis complete and verified against source CSVs by the author.
- All external-reviewer analytical items addressed: ICU definition corrected;
  BMI adjustment; patient clustering; hip/knee subgroups + interaction;
  exposure definition confirmed (mutually exclusive categories); calendar-year
  limitation stated honestly.
- Reference list complete and sequential (1–15), no orphaned citations.
- Data Availability links present and well-formed.

## Status: REMAINING author to-dos (not analysis)
1. **Verify the 3 data-availability links are LIVE** (open in browser): Zenodo DOI
   10.5281/zenodo.22051133, GitHub ZaneSalman/inspire-arthroplasty-anesthesia,
   PhysioNet INSPIRE v1.4.2. A dead link here is a common submission failure.
2. **Co-author sign-off** on the current version — especially the hip/knee subgroup
   claim and the removal of refs 16–17 (those were co-authors' papers).
3. **Journal choice** (still open). Best-fit candidates from scope + institutional
   OA agreements: Journal of Orthopaedic Surgery (SAGE; APC ~1,600 USD with the
   Gold-discount agreement), Archives of Orthopaedic and Trauma Surgery (Springer;
   ref [11] published a near-identical study there), or BMC Anesthesiology.
   Reformatting to the chosen journal is mechanical and quick.
4. **Sharpen the "what this adds" paragraph** — the field is well-populated; lead on
   the distinct Korean cohort, longer baseline LOS, ICU utilization outcome, and the
   honest subgroup analysis.
5. Cover letter + point-by-point reviewer response (map cleanly to the 9 items).

## Honest caveats that MUST remain in the paper
- LOS and ICU admission are administrative markers of resource use, NOT direct
  functional/kinesiologic recovery measures.
- Observational; anesthesia non-random (confounding by indication). IPTW/adjustment
  handle measured confounders only; residual confounding possible. Language stays
  associational, never causal.
- Small general group (n=269), single-center, single-ethnicity.
- Knee ICU subgroup underpowered (35 events).
