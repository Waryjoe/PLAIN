# Running the Code

## Environments

Run from the repository root. Use separate disposable environments for
historical extraction and general analysis. Dependencies are listed under
`environment/`; these are not a solver-verified lockfile. Do not install into
the only working clinical research environment just to run the demo.

The study source reports PyRadiomics 3.0.1 and SimpleITK 2.5.5 for the locked
extraction setup. Compatible Python, NumPy, compiled packages and platform
availability still need a clean-environment test. Historical R model objects
may depend on their original package versions and are not included.

## Portable Demo

See the top-level README. It requires standard-library Python and R `pROC`.
Tests deliberately operate on fictitious data, not on original patients.

## Private Path Configuration

The Python and R path helpers recognize:

| Variable | Use |
| --- | --- |
| `PLAIN_PROJECT_ROOT` | Working copy of the original private project tree |
| `PLAIN_HEATMAP_ROOT` | Working heatmap inputs and outputs |
| `PLAIN_LEGACY_WORK_ROOT` | Other authorized legacy working inputs, if needed |
| `PLAIN_WORKSPACE_ROOT` | Other authorized analysis workspace inputs, if needed |
| `PLAIN_DESKTOP_OUTPUTS_ROOT` | Non-project derived outputs |
| `PLAIN_PYTHON` | Python executable used by the R SHAP workbook helper |

Use standard `R_LIBS_USER` for your installed R libraries. Default data/output
locations are under the empty ignored `private_workspace/` directory. No script
downloads private data. Some historical scripts overwrite derived results;
use a working copy and a dedicated output location.

```bash
export PLAIN_PROJECT_ROOT="/path/to/authorized-working-copy"
export PLAIN_HEATMAP_ROOT="/path/to/private-heatmap-work"
export PLAIN_PYTHON="/path/to/analysis-python"
```

Original relative directory/sheet names inside that project remain intact to
avoid silently changing data joins. See path definitions near the beginning of
each archived script. Do not commit the private directory or environment file.

## Historical Order

1. Verify the private cohort, identifiers, fixed patient split, source masks and
   coding. The historical split script is for provenance; do not regenerate a
   different split for the current paper.
2. On the current manuscript branch, run v5 preprocessing, extraction and matrix
   organization. Resolve the anonymized exceptional-case mapping first. v5
   batch extraction expects the original 820-case blinded manifest. No GRAPE
   segmentation model is bundled.
3. Run training-only clinical candidate screening, collinearity and VIF scripts;
   run clinical model CV and internal-validation scripts using the same inputs.
4. Run radiomics LASSO/stability selection with explicit v5 inputs and the
   documented ten-fold setting. The archived default is five folds; it must
   not be mistaken for the manuscript ten-fold call:

```bash
Rscript scripts/select_dual_roi_radiomics_lasso_stability.R \
  --training-input=/path/to/private-v5-training.xlsx \
  --validation-input=/path/to/private-v5-validation.xlsx \
  --output-dir=/path/to/private-selection-output \
  --folds=10 --roi=combined
```

5. Run `benchmark_v5_radiomics_ml_and_clinical_fusion.R`, then
   `create_final_10feature_radiomics_evaluation_package.R`. Inputs include saved
   clinical OOF predictions and historical external radiomics matrices. These
   stages are not an automatic full-pipeline runner; inspect required source
   files and acknowledge selection/tuning limits before execution.
6. Run the relevant clinical/radiomics SHAP or heatmap scripts with the frozen
   models and verified matching tables. The heatmap script supports the
   historical `--t4-sensitivity` option; this is separate from main model fitting.
7. For a new transparent evaluation of existing predictions, use
   `evaluate_predictions.R` and the schema in DATA_SCHEMA.md. This utility does
   not reproduce the historical subset search or change cohort membership.

## Optional Historical Liaoning Display

`20261013_liaoning_standalone_figures.R` expects a saved selected-cohort source
folder (`03_Source_Data`) and an output directory with `SVG`, `PDF`, and
`Previews` subfolders already created. Third argument `auc-ci` adds DeLong CIs.
It uses macOS `quartz()` and Arial, so this archival display is platform-specific.
The generic evaluation utility is the cross-platform alternative, not an
identical re-render of that panel. Keep the script's post-hoc-selection captions.
