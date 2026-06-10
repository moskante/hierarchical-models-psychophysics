# hierarchical-models-psychophysics

R scripts, Stan files, and data for the article:
**"Modeling Psychophysical Data in R: A Comparative Study of Four Model Frameworks"**

---

## Overview

This repository accompanies the paper: ``Modeling Psychophysical Data in R: A Comparative Study of Four Model Frameworks'' and provides scripts for analyzing psychophysical data using the statistical frameworks described:

| Model | Type | Scripts (R folder) |
|-------|------|--------|
| GLM | Generalized Linear Model (single-subject) | `simul_data_single_sub.R`, `vibro_single_sub.R` |
| GNM | Generalized Nonlinear Model (single-subject, guessing + lapsing) | `simul_data_single_sub.R`, `vibro_single_sub.R` |
| GLMM | Generalized Linear Mixed Model (hierarchical, frequentist) | `simul_data_hierarchic.R`, `vibro_hierarchic.R` |
| BH-GLM / BH-GNM | Bayesian Hierarchical GLM and GNM (guessing + lapsing) | `simul_data_hierarchic.R`, `vibro_hierarchic.R` |

The `R/` folder includes scripts for gnm utilities. 
The Stan models used in `simul_data_hierarchic.R`, `vibro_hierarchic.R` are in the `Stan/` folder. 
The `simulations/` folder includes code for reproducing analysis on simulated data included in manuscript. 


The key outcome measures for all models are the **Point of Subjective Equality (PSE)** — the stimulus intensity at which both comparison stimuli are judged equally likely — and the **Just Noticeable Difference (JND)** — an index of discrimination sensitivity derived from the psychometric function slope.

---

## Datasets

- **Simulated data** — generated with `PsySimulate()` from the [MixedPsy](https://cran.r-project.org/package=MixedPsy) R package (10 subjects, 160 trials each, includes guessing and lapsing).
- **Vibrotactile data** (`vibro_exp3`) — from Dallmann et al. (2015), also bundled in the MixedPsy package. Nine subjects discriminated speed differences under two vibration conditions (0 Hz and 32 Hz).

---

## Prerequisites

Install the required R packages before running any script:

```r
install.packages(c("tidyverse", "MixedPsy", "gnlm", "lme4", "lmerTest",
                   "DHARMa", "rstan", "shinystan", "boot", "patchwork"))
```

Stan models require a working C++ toolchain. See the [RStan Getting Started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

---

## Repository Structure

```
R/
  gnlm_functions_psejnd.R   # GNM helper functions: PSE and JND estimation
  gnlm_functions_slope.R    # GNM helper functions: slope estimation (vibrotactile)
  k_correction.R            # Utility: convert p(correct) to sensitivity (k statistic)
  simul_data_single_sub.R   # Example 1a: GLM and GNM on simulated data
  simul_data_hierarchic.R   # Example 1b: GLMM, BH-GLM, BH-GNM on simulated data
  vibro_single_sub.R        # Example 2a: GLM and GNM on vibrotactile data
  vibro_hierarchic.R        # Example 2b: GLMM, BH-GLM, BH-GNM on vibrotactile data
  
simulations/
  iter.R                    # Monte Carlo simulation: cross-method comparison

Stan/
  simul_bhglm.stan          # Bayesian Hierarchical GLM (simulated data)
  simul_bhgnm.stan          # Bayesian Hierarchical GNM (simulated data)
  vibro_bhglm.stan          # Bayesian Hierarchical GLM (vibrotactile data)
  vibro_bhgnm.stan          # Bayesian Hierarchical GNM (vibrotactile data)
```

---

## How to Use

### Single-subject analysis (simulated data)

Run `R/simul_data_single_sub.R` to fit GLM and GNM models subject-by-subject on simulated psychophysical data. This script:

1. Generates 10 simulated subjects with known population parameters (PSE, JND, guessing rate, lapse rate).
2. Fits a probit **GLM** to each subject using `PsychModels()` / `PsychParameters()` from MixedPsy.
3. Fits a **GNM** to each subject using `gnlm_functions_psejnd.R` (accounts for guessing and lapsing).
4. Tests whether estimated PSE and JND differ from the true generating values (one-sample *t*-tests).

### Hierarchical analysis (simulated data)

Run `R/simul_data_hierarchic.R` to fit hierarchical models on the same simulated data. This script:

1. Fits a probit **GLMM** with random intercepts and slopes (`lme4::glmer`).
2. Runs two **Bayesian Hierarchical** models (BH-GLM and BH-GNM) via Stan (`Stan/simul_bhglm.stan`, `Stan/simul_bhgnm.stan`).
3. Computes residuals and Sum of Squared Errors (SSE) for each model.
4. Runs DHARMa diagnostic plots to assess model fit.

> **Note:** The Stan models require `Stan/simul_bhglm.stan` and `Stan/simul_bhgnm.stan` to be compiled on first run. Compilation may take a few minutes.

### Single-subject analysis (vibrotactile data)

Run `R/vibro_single_sub.R` to apply GLM and GNM models to the real `vibro_exp3` dataset. This script:

1. Loads `vibro_exp3` from the MixedPsy package (9 subjects, 2 vibration conditions: 0 Hz and 32 Hz).
2. Fits a probit **GLM** per subject-by-condition combination, extracts PSE, JND, and slope.
3. Fits a **GNM** using `gnlm_functions_slope.R` (slope-parameterized variant with CI from 500 bootstrap replicates).
4. Tests whether vibration affects discrimination slope (paired *t*-tests: 32 Hz vs. 0 Hz).

### Hierarchical analysis (vibrotactile data)

Run `R/vibro_hierarchic.R` to fit hierarchical models to the vibrotactile data. This script:

1. Fits a probit **GLMM** with a speed × vibration interaction and random slopes.
2. Computes bootstrap CIs for PSE, JND, slope, and condition differences via `pseMer()`.
3. Runs **BH-GLM** and **BH-GNM** in Stan (`Stan/vibro_bhglm.stan`, `Stan/vibro_bhgnm.stan`).
4. Validates each model with DHARMa residual diagnostics.

### Monte Carlo simulation study

Run `simulations/iter.R` to reproduce the cross-method comparison reported in the paper. This script:

1. Repeatedly simulates new datasets and fits all five models (GLM, GNM, GLMM, BH-GLM, BH-GNM).
2. Collects exactly 150 successful fits per method (failed fits are retried on a new dataset).
3. Computes bias, RMSE, and Type I error rates relative to the true generating parameters.

> Set `run_stan <- FALSE` at the top of `iter.R` to skip the Bayesian models for a quick GLM/GNM/GLMM-only run.

---

## Helper Files

- **`R/gnlm_functions_psejnd.R`** — defines the GNM psychometric function, starting-value heuristics, parameter extraction, and bootstrap wrapper used in Examples 1 and the simulation loop. Source this file before using `process_subject()`.
- **`R/gnlm_functions_slope.R`** — slope-parameterized counterpart for the vibrotactile analysis. Source before using `process_subject_vibro()`.
- **`R/k_correction.R`** — defines `compute_k()`, which converts proportion correct to the underlying sensitivity metric *k* after correcting for guessing and lapsing.

---


## Reference

Dallmann, C. J., Ernst, M. O., & Moscatelli, A. (2015). The role of vibration in tactile speed perception. *Journal of Neurophysiology*, 114(6), 3131–3139.
