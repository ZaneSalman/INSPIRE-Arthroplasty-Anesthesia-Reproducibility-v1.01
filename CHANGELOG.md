# Changelog

## v1.02 — 16 October 2026

### Corrected
- **ICU admission definition.** v1.01 counted any non-missing `icuin_time`, which included ICU
  admissions during later hospitalizations (96 of 154). ICU admission is now restricted to
  admission between the index hospitalization's admission and discharge times (58 events).
  - Adjusted OR changed from 0.494 (0.315–0.775) to 0.472 (0.239–0.930).
  - IPTW OR changed from 0.440 to 0.370; all-orthopedic OR from 0.664 to 0.471.
  - ASA ≥3 OR for ICU changed from 5.29 to 10.07.
- **In-hospital death** restricted to deaths recorded during the index hospitalization (7 vs 17).

### Added (not part of the original analysis plan)
- KDIGO creatinine-defined acute kidney injury (with staging), 30- and 90-day mortality, and a
  composite adverse event (AKI, ICU admission, or in-hospital death).
- Baseline creatinine, creatinine testing frequency, ICU timing, E-values, sex-disaggregated
  results and anesthesia-by-sex interactions.
- Charlson Comorbidity Index from pre-operative ICD-10 codes; comorbidity-adjusted models,
  IPTW and overlap weighting including comorbidity, an elective ASA 1–2 subgroup, and
  g-computation standardized risks with patient-level bootstrap CIs.
- Descriptive counts of diagnosis-coded postoperative complications.
- New 600 dpi figures and regenerated supplementary figures (the v1.01 Figure S3 showed the
  uncorrected ICU rates).

### Unchanged
- Cohort (2540 operations), primary LOS estimate (ratio 0.892, 0.827–0.963), and LOS
  sensitivity analyses. The v1.01 ICU and LOS estimates are reproduced exactly by script 03/04
  as a validation step before the corrected definition is applied.

### Reorganized
- v1.01 script, figures, and manuscript moved to `archive_earlier_versions/`.

## v1.01 — 29 August 2026
- Initial archived release (Zenodo 10.5281/zenodo.22166599).
