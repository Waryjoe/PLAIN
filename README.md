# PLAIN

Core research code for Preoperative Leakage Assessment Integrating clinical data
and Noncontrast CT. This minimal archive contains only model generation, CT/ROI
preprocessing, and stomach/tumor dual-ROI radiomics extraction.

## Included workflows

| Workflow | Script in `scripts/` |
| --- | --- |
| CT/ROI alignment, cropping, resampling and training-derived intensity standardization | `preprocess_all_cases_blinded_v5_grape_aligned.py` |
| PyRadiomics extraction for stomach and tumor ROIs | `extract_pyradiomics.py` |
| Batch dual-ROI extraction and feature-count checks | `run_ct2_pyradiomics_v5_grape_batch.py` |
| Clinical model training and out-of-fold predictions | `benchmark_nine_preoperative_clinical_models_training_oof.R` |
| Clinical model fitting and held-out internal evaluation | `benchmark_nine_preoperative_clinical_models.R` |
| Radiomics model comparison and logistic clinical-radiomics fusion | `benchmark_v5_radiomics_ml_and_clinical_fusion.R` |

The model scripts cover logistic regression, BPNN, decision trees, GBDT, naive
Bayes, KNN, PLS-DA, random forest, SVM, LightGBM, XGBoost, and Boruta plus logistic
regression. Fusion uses the logits of clinical and radiomics probabilities.

`plain_paths.py` and `plain_paths.R` provide portable paths. `configs/` contains
the v5 extraction settings (1,595 features per ROI); `environment/` lists required
Python and R packages. Analytical formulas and core scripts are unchanged.

## Running

Run from the repository root. Install Python dependencies from
`environment/requirements-radiomics.txt` and the R packages listed in
`environment/r-packages.txt`. These lists are not a fully validated environment lock.

Set `PLAIN_PROJECT_ROOT` to a working copy of authorized private inputs. Defaults
use `private_workspace/project/`. Some scripts overwrite derived outputs, so do
not use the only copy of source data.

The historical scripts require their original workbook/manifest schemas. Clinical
training requires the preselected nine-variable candidate workbook; radiomics
training requires assembled feature matrices and a precomputed training-only
LASSO stability ranking. Fusion also requires clinical out-of-fold and validation
predictions. These private inputs and upstream screening outputs are not included.
The cohort-specific batch script retains its original 820-case checks. The v5
configuration must not be substituted for a separately locked v1 protocol.

## Scope and privacy

This is a code-only archive, not a one-command reproduction of all final manuscript
results or a clinical-use model. Historical feature selection/tuning is not fully
nested within all out-of-fold estimates. Previously post-hoc-selected cohorts do
not establish independent external validation; final cohort provenance and figure
values still require author reconciliation. No patient records, CT images, ROI
masks, pathology images, patient-level predictions, trained models, or private
linkage files are distributed. `.gitignore` is a safeguard, not de-identification.
No software license or final article DOI has been assigned in this archive.
