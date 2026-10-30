# PLAIN: Clinical and Dual-ROI Radiomics Code

Code supporting **Preoperative Leakage Assessment Integrating clinical data and
Noncontrast CT (PLAIN)**, a gastric-cancer anastomotic-leakage research project.

**Release status: curated code-only prerelease, 2026-10-29.** This is an archive
of selected analysis scripts, not a claim that every current manuscript number
has been reproduced. No clinical data, CT images, masks, pathology slides,
patient-level predictions, fitted models, or private linkage tables are included.
The synthetic example is fictitious. This repository is for research, not
clinical decision-making.

[中文说明](README.zh-CN.md) | [Analysis map](docs/ANALYSIS_MAP.md) |
[Run guide](docs/RUNNING.md) | [Reproducibility status](docs/REPRODUCIBILITY.md)

## Contents

- `scripts/`: 23 curated historical scripts, portable path helpers, and a new
  prediction-evaluation utility.
- `configs/`: byte-preserved v5 extraction configuration and a separately
  identified legacy v1 configuration. They are **not interchangeable**.
- `examples/`: deterministic synthetic-data generator, no study records.
- `tests/` and `tools/`: checks for syntax, path resolution, evaluation functions,
  file integrity, and common publication-privacy hazards.
- `environment/`: dependency lists, explicitly not a fully validated lockfile.
- `SOURCE_MANIFEST.csv`: original and packaged source hashes, transformation
  counts, and analysis scope. Original absolute locations stay private.

The original script filenames and analytical formulas were retained. Packaging
changes replace machine-specific paths with `plain_path(...)` and replace
literal private case identifiers with non-identifying placeholders. See
[reproducibility notes](docs/REPRODUCIBILITY.md) for the limitations of these edits.

## Quick Synthetic Demo

Run commands from the repository root. Python 3.9+ is sufficient for the demo
generator and package tests. R with `pROC` is needed for evaluation.

```bash
python3 -m unittest discover -s tests -v
python3 tools/check_release.py
python3 examples/make_synthetic_predictions.py outputs/demo_predictions.csv
Rscript scripts/evaluate_predictions.R outputs/demo_predictions.csv outputs/demo_evaluation
Rscript tests/test_evaluation.R
```

Install `pROC` in your own R environment if it is not already available:

```r
install.packages("pROC")
```

The example writes AUC/DeLong confidence intervals, Brier scores, ROC points,
grouped calibration, DCA values, and a vector PDF. It does not train PLAIN,
reproduce its reported performance, alter probabilities, smooth DCA, or select
patients. Use a new output directory for each run.

## Original Analyses

The archived scripts require the authorized private input workbooks/manifests,
their original sheet/column schemas, and compatible Python/R dependencies.
They are deliberately not run on synthetic inputs with invented study results.
Set `PLAIN_PROJECT_ROOT` to a **working copy** of the private study directory,
then follow [RUNNING.md](docs/RUNNING.md). Do not point it at the only copy of
source data: some historical scripts overwrite derived outputs.

The current manuscript branch uses the v5 preprocessed radiomics matrix. The
legacy v1 protocol remains unchanged and is provided for provenance only; do
not silently substitute it for v5 or use v5 for new v1-locked extractions.

## Important Interpretation Limits

The Ningbo 133-case and Liaoning 355-case displays involve historical post-hoc,
outcome/performance-informed subset selection. Standard CIs on these subsets
do not account for that selection and cannot establish independent external
validation. These subsets must not be described as ordinary missing-data
exclusions without source evidence. The subset plotting script retains this
qualification; the optimization scripts and private membership lists are not
part of this curated main-analysis package.

Code organization does not resolve these design issues. Current manuscript
claims, final cohort provenance, and manually updated figure values require
author reconciliation before a publication-linked release.

## Data, Rights, and Citation

No patient data are distributed. Access requires institutional authorization;
this package makes no promise of data availability. Do not commit private
workspaces, outputs, models, slide labels, or notebook outputs. `.gitignore` is
only a convenience, not de-identification. Review staged files before pushing.

The authors have not selected a software license or supplied a final citation
for this release. No open-source license or DOI is asserted here. See
[LICENSING.md](LICENSING.md) and [PUBLICATION_CHECKLIST.md](docs/PUBLICATION_CHECKLIST.md).
