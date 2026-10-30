# Data Schema and Privacy

## Portable Prediction Evaluation

UTF-8 CSV, one row per unique patient across disjoint cohorts:

| Column | Meaning |
| --- | --- |
| `patient_id` | Locally generated non-identifying ID; never a hospital/chart ID |
| `cohort` | Explicit label such as training OOF, internal, independent external, or exploratory subset |
| `outcome` | Numeric 0 = no anastomotic leakage, 1 = leakage |
| `clinical_probability` | Frozen clinical-model probability in [0,1] |
| `radiomics_probability` | Frozen radiomics-model probability in [0,1] |
| `fusion_probability` | Frozen fusion-model probability in [0,1] |

No silent deletion, duplicate aggregation, threshold optimization, clipping of
input probabilities, or refitting is performed by the evaluation entrypoint.
Missing/non-finite values, duplicates, invalid probabilities, and single-class
cohorts fail validation. Tied predictions stay together in calibration bins.

The synthetic demo is not derived from real patient records. Its scores are
analytical functions of random variables, not outputs of trained PLAIN models.

## Archived Inputs

Historical workbooks use their original English/Chinese sheet names and columns.
The scripts retain their schema checks (including historical row counts). They
are not generic tools for any arbitrary workbook.

Clinical modeling maps nine candidates to `x1`-`x9`: PT, cT stage,
histological differentiation, LDH, TT, basophil percentage, age, HDL-C, and A/G
ratio. Some original input headers still use `T1` or numerical differentiation
codes; preserve the source coding and confirm its dictionary. This release does
not guess mappings or relabel source values automatically.

Radiomics matrices use `stomach__...` and `tumor__...` feature names. Both ROIs
have 1,595 features. The private mapping links split membership to images;
author-supplied data must retain a consistent anonymized key without exposing
patient identities. Source paths emitted by extraction are private metadata.

## Before Sharing Any Data

Pseudonymization alone is not anonymization. Clinical covariates, dates, rare
outcomes, image headers, mask metadata, slide labels, and saved predictions may
still identify individuals. Institutional authorization and a separate release
review are required. This repository includes no example clinical rows taken
from the study and no consent/ethics approval claim for data sharing.
