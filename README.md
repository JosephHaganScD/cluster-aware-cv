# cluster-aware-cv
R code for cluster-aware-cross-validation simulation study and motivating data analysis (Hagan, submitted)
## Requirements

- R (≥ 4.5.0 recommended)
- Required packages: `glmnet`, `pROC`, `ggplot2`

Install packages if needed:

```r
install.packages(c("glmnet", "pROC", "ggplot2"))
```

## Reproducibility Workflow

### 1. Simulation Study

```r
source("simulation_factorial.R")
```

Generates simulated datasets across 162 conditions, performs cross-validation analyses, and saves results to `results/simulation/sim_results_v3.csv`. Note: this step is computationally intensive (~20–40 hours depending on hardware).

### 2. Generate Figures and Summary Statistics

```r
source("sim_figures.R")
```

Reads simulation results and saves manuscript figures to `results/figures/`.

### 3. Clinical Illustration

Place the clinical dataset at `data/final_daily_data.csv`, then run:

```r
source("clinical_analysis.R")
```

Fits ridge logistic regression models using multiple CV strategies and saves results to `results/rop_analysis/`.

## Data Availability

The clinical dataset used in this analysis is not publicly available due to patient privacy restrictions. The simulation code is fully reproducible and does not require external data. The clinical dataset may be requested from the author.

## Methods Summary

- **Model:** Ridge logistic regression (glmnet)
- **Outcome:** Subject-level binary outcome
- **Predictors:** Repeated-measures (longitudinal) variables
- **Evaluation metric:** AUROC
- **Primary quantity:** Optimism (naive CV AUROC minus LOCO AUROC)

**Cross-validation strategies compared:**
- Naive observation-level 10-fold CV
- Subject-level (cluster-aware) 10-fold CV
- Leave-one-cluster-out (LOCO) CV (reference standard)

## Reproducibility Note

This repository contains code corresponding to the submitted manuscript version. The master RNG seed is set at the top of `simulation_factorial.R` (`set.seed(20260410)`); results are fully reproducible given the same R and package versions.

## Author

Joseph L. Hagan, ScD, MSPH  
Department of Pediatrics, Section of Neonatology  
Baylor College of Medicine, Houston, Texas

## License

MIT License. Please cite the associated manuscript if this code is used in academic work.
