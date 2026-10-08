# Cluster-Aware Cross-Validation for Prediction Models with Repeated-Measures Predictors

This repository contains R code supporting the manuscript:

"Quantifying the Optimism of Naive Cross-Validation for Binary Outcome Prediction with Repeated-Measures Predictors: A Simulation Study and Clinical Illustration." Joseph Hagan, ScD MSPH (ORCID: 0000-0002-3285-2131). Preprint: medRxiv, doi:10.64898/2026.05.27.26354222. Submitted to Diagnostic and Prognostic Research.

## Overview

This project investigates the performance of cross-validation (CV) strategies for prediction models when predictors are measured repeatedly over time but the binary outcome is defined at the subject level. In this setting, standard observation-level CV can introduce data leakage and lead to optimistic estimates of model performance. The motivating example ROP clinical dataset from Srivatsa et al., 2022 is not publicly available, but the simulation pipeline and figures can be reproduced independently without it.

## Repository Structure

```
.
├── simulation_factorial_v5.R   # Full factorial simulation study (pilot, main, sensitivity)
├── sim_figures.R               # Manuscript Figures 1-4 (reads deposited results)
├── clinical_analysis_v5.R      # ROP clinical illustration (requires non-public data)
├── results/
│   ├── sim_results_main.csv.gz         # Deposited per-dataset results, main grid (81,000 datasets, 162 conditions; gzip-compressed)
│   ├── sim_results_sensitivity.csv     # Deposited per-dataset results, sensitivity arm
│   ├── sim_condsummary_main.csv        # Condition summary (created by sim_figures.R if absent)
│   └── figures/                        # Figures (created by sim_figures.R)
├── LICENSE
└── README.md
```

## Reproducing the results

**Figures and summary tables (no simulation rerun needed).** Set the R working directory to the repository root and run `sim_figures.R`. If `results/sim_condsummary_main.csv` is not present, the script rebuilds it from `results/sim_results_main.csv.gz` using the same aggregation as `summarize_run()` in `simulation_factorial_v5.R`, and then writes Figures 1 to 4 (PDF and PNG) to `results/figures/`. Required packages: ggplot2, dplyr, patchwork, scales.

**Full simulation.** `simulation_factorial_v5.R` is controlled by `RUN_MODE` near the top of the script ("pilot", "main", or "sensitivity"). The pilot must be run first because it generates the frozen coefficient magnitude used by the other two modes. The output directory (`OUT_DIR`) is set by the user near the top of the script and should be a local, non-OneDrive folder for long runs. Required packages: glmnet, pROC. The main grid requires many hours of computation; the deposited per-dataset results allow all summaries and figures to be regenerated without rerunning it.

**Clinical illustration.** `clinical_analysis_v5.R` requires `final_daily_data.csv`, which is not publicly available. The input path (`DATA_FILE`) and output folder (`OUTPUT_DIR`) are set by the user near the top of the script.

## Notes on file names

Earlier commits of this repository contained scripts named `simulation_factorial.R` and `clinical_analysis.R`. These are now named `simulation_factorial_v5.R` and `clinical_analysis_v5.R`, which are the versions used for the results reported in the manuscript. The figure script, previously uploaded as `Figures.R`, is now named `sim_figures.R`.
