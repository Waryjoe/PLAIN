# Analysis Map

Manuscript panel numbers have changed during editing. This map uses analysis
roles rather than promising an exact Figure-number-to-file match.

| Role | Archived scripts | Inputs and scope |
| --- | --- | --- |
| CT preprocessing | `preprocess_all_cases_blinded_v5_grape_aligned.py` | Authorized CTs, fixed ROI masks, private case mapping; includes a masked exceptional-case branch |
| Feature extraction | `extract_pyradiomics.py`, `run_ct2_pyradiomics_v5_grape_batch.py` | v5-preprocessed images/masks and configuration; 1,595 features per ROI |
| Matrix assembly | `organize_ct2_v5_grape_radiomics_train_validation.py` | Private split mapping and feature workbooks |
| Patient split | `split_train_validation.py` | Historical 7:3 split; do not resplit the final validation cohort |
| Clinical selection | `screen_training_univariable_candidates.py`, `check_training_candidate_collinearity.py`, `run_training_priority_vif.R` | Training workbook and candidate/VIF outputs |
| Shared screening helper | `screen_training_features.py` | Imported by clinical screening; standalone technical-confounding analysis is historical/exploratory |
| Baseline tables | `build_table1_artifacts.py`, `generate_table1.R` | Historical descriptive-table generator; not the latest manually edited four-cohort Word tables |
| Clinical models | `benchmark_nine_preoperative_clinical_models_training_oof.R`, `benchmark_nine_preoperative_clinical_models.R` | Nine clinical variables; repeated-CV and internal-validation code |
| Clinical SHAP | `generate_lightgbm_performance_and_shap_figures.R` | Authorized clinical inputs, saved fitted model and OOF/validation predictions |
| Radiomics selection | `select_dual_roi_radiomics_lasso_stability.R` | Training-only feature ranking; pass v5 matrix and `--folds=10` explicitly |
| Feature count comparison | `compare_v5_dual_roi_fixed_feature_counts.R` | Historical fixed-feature-count comparison; do not use its legacy Rad-score results as LightGBM results |
| Radiomics and fusion | `benchmark_v5_radiomics_ml_and_clinical_fusion.R` | 7/10/15/20/30-feature ML comparisons and logistic combination of clinical/radiomics probability logits |
| Ten-feature evaluation | `create_final_10feature_radiomics_evaluation_package.R` | Saved model/prediction workbooks from the benchmark |
| Radiomics SHAP | `radiomics_shap_locked_20260911.R`, `finalize_radiomics_shap_workbook.py` | Selected features, fixed model, training/internal matrix; helper checks Excel output |
| Training heatmap | `20261021_training_heatmap_prepare.R`, `20261021_training_heatmap_plot.R` | Clinical/radiomics linkage, saved probabilities; includes explicitly labeled T4 sensitivity display |
| Liaoning subset display | `20261013_liaoning_standalone_figures.R` | Historical 355-case selected cohort; not independent external validation; requires saved source tables |
| Portable evaluation | `evaluate_predictions.R`, `evaluation_helpers.R` | Newly added release utilities, not historical training code; schema in DATA_SCHEMA.md |

## Not Included in This Core Package

- Raw patient linkage, private extraction/clinical-input workbooks and individual results.
- Outcome-driven deletion/search scripts and selected membership lists. Their existence
  and implications are disclosed in REPRODUCIBILITY.md; omission is not evidence that
  the reported external subsets were selected independently of outcomes.
- Third-party segmentation/pathology implementations, downloaded dependencies and weights.
- HE inference/region selection, CT case-display scripts with private case bindings,
  editable Adobe files, and transient desktop automation. These require separate
  provenance/privacy review. One patient per outcome supports descriptive case
  illustration, not patient-level group inference from repeated regions.
- Superseded clinical-plus-Rad-score fusion branches that are not the two-LightGBM
  probability-logit fusion described in the current manuscript.

Use `SOURCE_MANIFEST.csv` for the exact curated files and hashes, not a claim
that any similarly named historical folder is the final analysis.
