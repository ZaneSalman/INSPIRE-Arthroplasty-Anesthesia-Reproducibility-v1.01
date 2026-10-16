# INSPIRE Arthroplasty Anesthesia Reproducibility Repository

**Version 1.02** (16 October 2026). Supersedes v1.01 ([10.5281/zenodo.22166599](https://doi.org/10.5281/zenodo.22166599)).

Analysis code and aggregate figures for:

> **Neuraxial versus General Anesthesia and Postoperative Adverse Events after Hip and Knee
> Arthroplasty: A Retrospective Cohort Study Using INSPIRE.**
> Salman ZG, Zeitouny R, Panwalker N, et al. Manuscript prepared for submission to *Complications* (MDPI).

The study compares neuraxial with general anesthesia in 2540 hip or knee arthroplasty
operations (269 general, 2271 neuraxial; 1993 patients) from INSPIRE v1.4.2, a de-identified
perioperative dataset from Seoul National University Hospital (2011–2020).

> **Important change from v1.01.** In v1.01, ICU admission was defined as any non-missing
> `icuin_time`. INSPIRE records the first ICU admission after an operation even when it occurs in
> a later hospitalization, so that definition counted 154 "admissions," of which only 58 occurred
> during the index (arthroplasty) hospitalization; 96 occurred after discharge (median ~2.6 years
> later). From v1.02, ICU admission is restricted to the index hospitalization. All ICU-based
> results were re-estimated. See [CHANGELOG.md](CHANGELOG.md).

## Key results (v1.02)

| Outcome | General | Neuraxial | Primary adjusted estimate (95% CI) | + Charlson comorbidity |
|---|---|---|---|---|
| ICU admission, index hospitalization | 19/269 (7.1%) | 39/2271 (1.7%) | OR 0.472 (0.239–0.930) | OR 0.518 (0.264–1.013) |
| Acute kidney injury (KDIGO creatinine) | 18/268 (6.7%) | 55/2257 (2.4%) | OR 0.516 (0.287–0.926) | OR 0.691 (0.385–1.242) |
| Composite (AKI, ICU, or in-hospital death) | 34/268 (12.7%) | 88/2257 (3.9%) | OR 0.471 (0.292–0.759) | OR 0.569 (0.353–0.917) |
| Length of stay | median 12.0 d | median 11.0 d | ratio 0.892 (0.827–0.963) | ratio 0.912 (0.847–0.982) |
| In-hospital death | 2/269 (0.7%) | 5/2271 (0.2%) | not modeled (Fisher p = 0.164) | — |

Primary models adjust for banded age, sex, ASA class (1–2 vs ≥3), emergency status, and
procedure (hip vs knee), with HC0 robust standard errors. Comorbidity-adjusted standardized
composite risk: 7.0% (general) vs 4.4% (neuraxial), risk difference −2.6 percentage points
(95% CI −5.6 to −0.2). In elective ASA 1–2 operations (n = 2222), ICU (OR 0.269), composite
(OR 0.463) and LOS (ratio 0.916) associations persisted; AKI did not (OR 0.727).

**Interpretation.** Associations are observational. The general-anesthesia group had a much
higher comorbidity burden (Charlson ≥2: 31.6% vs 11.9%), and comorbidity explained part of the
association, most clearly for AKI. AKI, mortality, the composite, and comorbidity adjustment were
added after the original analysis plan and are hypothesis-generating.

## Running the analysis

Source data are **not** included. INSPIRE v1.4.2 is available to credentialed users on PhysioNet:
<https://physionet.org/content/inspire/1.4.2/> ([10.13026/1eay-yc85](https://doi.org/10.13026/1eay-yc85)).
Download `operations.csv.gz`, `labs.csv.gz`, and `diagnosis.csv.gz` into one folder.

```bash
pip install -r requirements.txt
export INSPIRE_DATA_DIR=/path/to/inspire_v1.4.2      # folder with the *.csv.gz files
export INSPIRE_WORK_DIR=/path/to/private_workdir     # local only; holds patient-level intermediates
bash run_all.sh
```

Scripts run in numbered order and print aggregate results to the console. Intermediate files
(`*.pkl`) in `INSPIRE_WORK_DIR` contain operation-level data covered by the INSPIRE data use
agreement; keep that folder private and never commit it. Figures are written to
`$INSPIRE_WORK_DIR/figures/`.

| Script | Purpose |
|---|---|
| `01_rebuild_cohort.py` | Builds the 2540-operation arthroplasty cohort |
| `02_extract_creatinine.py` | Extracts creatinine values for cohort patients |
| `03_aki_and_mortality_models.py` | Reproduces v1.01 ICU/LOS estimates (validation); KDIGO AKI, mortality |
| `04_aki_composite_sensitivity.py` | IPTW, BMI/clustered, first-operation, hip/knee analyses for AKI and composite |
| `05_baseline_creatinine_testing_icu_timing.py` | Baseline creatinine, creatinine testing frequency, E-values, sex-disaggregated results |
| `06_icu_timing_audit.py` | Audit showing ICU timestamps outside the index hospitalization |
| `07_corrected_icu_and_composite.py` | Corrected (index-hospitalization) ICU outcome and composite; all-orthopedic cohort |
| `08_final_estimates.py` | Final ICU covariate model, LOS sensitivity analyses, balance, E-values |
| `09_extract_diagnoses.py` | Extracts ICD-10 diagnoses for cohort patients |
| `10_charlson_comorbidity.py` | Charlson Comorbidity Index from pre-operative ICD-10 codes (3-character Quan mapping) |
| `11_coded_complications.py` | Counts of diagnosis-coded postoperative complications (descriptive) |
| `12_comorbidity_adjusted_models.py` | Comorbidity-adjusted, IPTW, overlap-weighted, low-risk subgroup, and g-computation analyses |
| `13_comorbidity_descriptives.py` | Comorbidity descriptives and balance |
| `14_main_figures.py` | Figures 1–3 and graphical abstract (from aggregate estimates) |
| `15_supplementary_figures.py` | Supplementary Figures S1–S3 |

All v1.02 estimates were reproduced end to end from the raw INSPIRE files with these scripts
before release (Python 3.10, NumPy 2.2.6, pandas 2.3.3, Matplotlib 3.10.9).

## Repository contents

- `code/` — v1.02 analysis scripts (above).
- `figures/` — Figures 1–3 (600 dpi), graphical abstract, and Supplementary Figures S1–S3. Aggregate only.
- `archive_earlier_versions/` — superseded material kept for provenance, including the v1.01
  script, figures, and manuscript drafts. **These use the uncorrected ICU definition; do not
  use their ICU numbers.**

## Ethics and data use

Secondary analysis of de-identified data under the PhysioNet credentialed data use agreement.
Creation and release of INSPIRE were approved by the Seoul National University Hospital IRB
(H-2210-078-1368) and Data Review Board (BD-R-2022-11-02). No patient-level data are included
in this repository.

## Use of generative AI

Claude Opus 5.5 (Anthropic) assisted with writing and running the v1.02 analysis code and
figures, and OpenAI Codex assisted with earlier formatting. All analyses were run locally under
the data use agreement; the authors are responsible for reviewing the code and results.

## Citation

See [CITATION.cff](CITATION.cff). A Zenodo DOI for v1.02 will be issued on release.
