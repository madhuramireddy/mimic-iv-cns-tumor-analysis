# MIMIC-IV CNS Tumor Odds Ratio Analysis

This repo contains R scripts to identify CNS tumor cohorts from MIMIC-IV,
flag antidepressant and statin exposures, and compute odds ratios.

## Run order
1. `Scripts/00_setup.R`
2. `Scripts/01_cns_tumor_cohort.R`
3. `Scripts/02_antidepressant_analysis.R`
4. `Scripts/03_statin_analysis.R`
5. `Scripts/04_summary_results.R`

*(Data/ and Output/ are intentionally ignored)*