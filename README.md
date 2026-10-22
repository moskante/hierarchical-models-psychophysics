# hierarchical-models-psychophysics

R scripts and Stan programs for the article:
**"Modeling Psychophysical Data in R: A Comparative Study of Four Model Frameworks"**

---

## Overview

The scripts reproduce the analyses in the paper. They fit psychometric functions to data from 2AFC discrimination tasks with the model frameworks compared in the article (acronyms as in Table 1 of the paper):

| Acronym | Model | Level |
|---------|-------|-------|
| GLM | Generalized Linear Model | single subject (two-level approach) |
| GNM | Generalized Non-Linear Model, with guess and lapse rates | single subject (two-level approach) |
| GLMM | Generalized Linear Mixed Model | hierarchical, frequentist |
| BH-GLM | Bayesian Hierarchical GLM | hierarchical, Bayesian |
| BH-GNM | Bayesian Hierarchical GNM, with guess and lapse rates | hierarchical, Bayesian |

The main outcome measures are the **Point of Subjective Equality (PSE)**, the **Just Noticeable Difference (JND)**, and the **slope** of linear predictor $\beta_1$.

Script names follow the pattern `example<N>_<models>.R`, where `<N>` is the example in the paper and `<models>` are the acronyms of the models fitted in the script.

---

## Datasets

No data files are needed: both datasets come from the [MixedPsy](https://CRAN.R-project.org/package=MixedPsy) package.

- **Example 1: simulated data.** Generated with `MixedPsy::PsySimulate()` (10 subjects, 160 trials, non-zero guess and lapse rates). Each script regenerates the same dataset with `set.seed(123)`.
- **Example 2: experimental data.** `vibro_exp3` from Dallmann et al. (2015), bundled in MixedPsy. Nine participants compared the speed of two tactile motion stimuli under two vibration conditions (0 Hz and 32 Hz).
- **Large simulation study.** 150 independent datasets simulated with `PsySimulate()`.

---

## Repository structure

```
R/
  example1_GLM_GNM.R                 # Example 1: GLM and GNM, per subject
  example1_GLMM_BHGLM_BHGNM.R        # Example 1: GLMM, BH-GLM and BH-GNM (Stan)
  example1_BHGLM_BHGNM_brms.R        # Example 1: BH-GLM and BH-GNM with brms
  example1_BHGLM_brms_lognormal.R    # Example 1: BH-GLM with brms (lognormal sigma), compared with Stan
  example2_GLM_GNM.R                 # Example 2: GLM and GNM, per subject and condition
  example2_GLMM_BHGLM_BHGNM.R        # Example 2: GLMM, BH-GLM and BH-GNM (Stan)
  large_simulation_study.R           # Large simulation study (150 datasets, all five models)
  GNM_functions_psejnd.R             # helper functions for the GNM (PSE/JND), sourced by the scripts
  GNM_functions_slope.R              # helper functions for the GNM (slope), sourced by the scripts

Stan/
  example1_BHGLM.stan                # BH-GLM, Example 1 and large simulation study
  example1_BHGNM.stan                # BH-GNM, Example 1 and large simulation study
  example1_BHGLM_noncentred.stan     # alternative BH-GLM: non-centred sigma (Example 1)
  example1_BHGNM_halfnormal.stan     # alternative BH-GNM: half-normal sigma, hierarchical gamma/lambda
  example2_BHGLM.stan                # BH-GLM, Example 2
  example2_BHGLM_lognormal.stan      # alternative BH-GLM: lognormal prior on the slopes
  example2_BHGNM.stan                # BH-GNM, Example 2

```

The helper files only define functions. They are sourced by the scripts that need them, so you don't need to run them yourself.

---

## Requirements

### R packages

```r
install.packages(c("tidyverse", "MixedPsy", "gnlm", "boot", "lme4", "lmerTest",
                   "DHARMa", "rstan", "shinystan", "brms", "posterior", "xtable", "gridExtra"))
```

The Bayesian models need a working C++ toolchain for Stan (Rtools on Windows, Xcode command line tools on macOS, a C++ compiler on Linux). See the [RStan Getting Started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

**StanHeaders version.** rstan 2.32.7 must be used with StanHeaders 2.32.x. Installing rstan from CRAN can bring in a newer StanHeaders (e.g. 2.39), and then Stan models fail to compile, at least on Windows. After installing the packages above, install the matching StanHeaders (this works on all platforms):

```r
install.packages("https://cran.r-project.org/src/contrib/Archive/StanHeaders/StanHeaders_2.32.10.tar.gz",
                 repos = NULL, type = "source")
```

The Stan programs use the current array syntax (`array[N] real x;`), so they also work with newer versions of Stan, CmdStan and cmdstanr.

### Reproducible environment (optional)

The file `renv.lock` records the exact versions of R packages used for testing. To install them in a project library with [renv](https://rstudio.github.io/renv/):

```r
install.packages("renv")
renv::restore()
```

### Tested versions

The code was tested with R 4.5.3 on Windows 11, and on macOS and Ubuntu, in RStudio. Package versions:

| Package | Version |
|---------|---------|
| MixedPsy | 1.3.0 |
| gnlm | 1.1.2 |
| boot | 1.3.32 |
| lme4 | 2.0.6 |
| lmerTest | 3.2.1 |
| DHARMa | 0.5.0 |
| rstan | 2.32.7 |
| StanHeaders | 2.32.10 |
| shinystan | 2.7.0 |
| brms | 2.23.0 |
| posterior | 1.7.0 |
| tidyverse | 2.0.0 |
| ggplot2 | 4.0.1 (≥ 3.5.0 required) |
| dplyr | 1.1.4 |
| stringr | 1.6.0 |
| xtable | 1.8.8 |
| gridExtra | 2.3.1 |

---

## How to run

1. Clone or download the repository.
2. Set the **repository root** as the working directory, e.g. by opening the folder as an RStudio project, or with `setwd("path/to/hierarchical-models-psychophysics")`. All paths in the scripts (`R/...`, `Stan/...`, `Figs/`) are relative to the root.
3. Run the scripts in the order below, which follows the paper. Each script is self-contained: it generates or loads its own data, so any script can also be run on its own.

Figures are saved as PDF in a `Figs/` folder, which is created automatically. Figures are not tracked by git.

### Example 1: simulated dataset

**Step 1. `R/example1_GLM_GNM.R`: GLM and GNM (two-level approach)**
1. Simulates the dataset and recovers the true generating PSE and JND of each subject.
2. Fits a probit **GLM** to each subject (`MixedPsy::PsychModels()`, `PsychParameters()`).
3. Fits a **GNM** with guess and lapse rates to each subject (`gnlm::gnlr()`, functions in `R/GNM_functions_psejnd.R`), with bootstrap confidence intervals.
4. Compares the estimated PSE and JND with the true values (one-sample *t*-tests) and plots the estimates.

**Step 2. `R/example1_GLMM_BHGLM_BHGNM.R`: hierarchical models**
1. Fits a probit **GLMM** with random intercepts and slopes (`lme4::glmer()`), with bootstrap confidence intervals (`MixedPsy::pseMer()`).
2. Fits the **BH-GLM** (`Stan/example1_BHGLM.stan`) and the **BH-GNM** (`Stan/example1_BHGNM.stan`) with `rstan`.
3. Computes the Sum of Squared Errors (SSE) and runs DHARMa diagnostics for each model.
4. Plots the population-level estimates of the three models.

The script also explains, in comments, how to fit the alternative models:
- `Stan/example1_BHGLM_noncentred.stan`: the same BH-GLM with a non-centred parameterisation of the subject-level sigma (Eq. 28 in the paper). The priors and the posterior are the same, but sampling is usually more efficient, with fewer divergent transitions. The paper presents the centred model because it is easier to read.
- `Stan/example1_BHGNM_halfnormal.stan`: half-normal prior on sigma and hierarchical guess and lapse rates.

**Step 3. `R/example1_BHGLM_BHGNM_brms.R`: Bayesian hierarchical models with brms**
Fits the BH-GLM and the BH-GNM with `brms`, the high-level alternative to native Stan code described in the paper, and computes the SSE and DHARMa diagnostics.

**Step 4. `R/example1_BHGLM_brms_lognormal.R`: brms model with lognormal sigma**
Fits the BH-GLM in `brms` with the same parametrization as the Stan model (lognormal prior on sigma, population-level PSE and JND), fits the Stan model for comparison, and runs DHARMa diagnostics.

### Example 2: experimental dataset (`vibro_exp3`)

**Step 5. `R/example2_GLM_GNM.R`: GLM and GNM (two-level approach)**
1. Fits a probit **GLM** to each subject and vibration condition, and extracts PSE, JND and slope.
2. Fits a **GNM** with guess and lapse rates (functions in `R/GNM_functions_slope.R`), with bootstrap confidence intervals.
3. Tests the effect of vibration on the slope (paired *t*-tests, 32 Hz vs 0 Hz) and plots the estimates.

**Step 6. `R/example2_GLMM_BHGLM_BHGNM.R`: hierarchical models**
1. Fits a probit **GLMM** with a speed × vibration interaction, with bootstrap confidence intervals for PSE, JND, slope and the differences between conditions (`pseMer()`).
2. Fits the **BH-GLM** (`Stan/example2_BHGLM.stan`) and the **BH-GNM** (`Stan/example2_BHGNM.stan`).
3. Computes SSE and DHARMa diagnostics, and plots the differences between conditions.

The script also explains, in comments, how to fit the alternative model `Stan/example2_BHGLM_lognormal.stan`.

### Large simulation study

**Step 7. `R/large_simulation_study.R`**
1. Simulates new datasets and fits all five models (GLM, GNM, GLMM, BH-GLM, BH-GNM) to each one.
2. Keeps a dataset only if all five models converge, until 150 datasets have been collected.
3. Computes bias and RMSE of PSE, JND, guess and lapse rates with respect to the generating values, plus the mean SSE.
4. Prints the summary tables and saves the violin plots of the estimates.

> This script fits 750 models and takes a long time to run.

---

## Priors of the Bayesian models

Prior distributions of all Bayesian hierarchical models. **Paper** marks the models presented in the article; the others are alternatives available in this repository. The second argument of N, LogNormal and half-Cauchy is a standard deviation; *i* indexes participants and *h* the vibration conditions of Example 2.

### Example 1: simulated data

| Model | File | Paper | Subject location | Subject scale / slope | Population location | Population scale / slope | Standard deviations | Guess γ / lapse λ |
|---|---|:---:|---|---|---|---|---|---|
| BH-GLM | `Stan/example1_BHGLM.stan` | ✓ | pse<sub>i</sub>&nbsp;∼&nbsp;N(PSE,&nbsp;τ<sub>pse</sub>) | σ<sub>i</sub>&nbsp;∼&nbsp;LogNormal(log&nbsp;Σ,&nbsp;τ<sub>σ</sub>) | PSE&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>PSE</sub>) | log&nbsp;Σ&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>Σ</sub>) | all τ&nbsp;∼&nbsp;half‑Cauchy(0,&nbsp;2.5) | – |
| BH-GLM, non-centred | `Stan/example1_BHGLM_noncentred.stan` |  | as BH-GLM | σ<sub>i</sub>&nbsp;=&nbsp;exp(log&nbsp;Σ&nbsp;+&nbsp;τ<sub>σ</sub>&nbsp;z<sub>i</sub>), z<sub>i</sub>&nbsp;∼&nbsp;N(0,&nbsp;1); same prior as BH-GLM | as BH-GLM | as BH-GLM | as BH-GLM | – |
| BH-GLM (brms) | `R/example1_BHGLM_BHGNM_brms.R` | ✓ | intercept b<sub>0i</sub>, normal | slope b<sub>1i</sub>, normal | β<sub>0</sub>&nbsp;∼&nbsp;N(0,&nbsp;3) | β<sub>1</sub>&nbsp;∼&nbsp;N(0,&nbsp;3) | sd&nbsp;∼&nbsp;half‑Cauchy(0,&nbsp;2.5); correlation&nbsp;∼&nbsp;LKJ(2) | – |
| BH-GLM, lognormal (brms) | `R/example1_BHGLM_brms_lognormal.R` |  | as BH-GLM (Stan) | as BH-GLM (Stan) | as BH-GLM (Stan) | as BH-GLM (Stan) | as BH-GLM (Stan) | – |
| BH-GNM | `Stan/example1_BHGNM.stan` | ✓ | μ<sub>i</sub>&nbsp;∼&nbsp;N(M,&nbsp;τ<sub>μ</sub>) | σ<sub>i</sub>&nbsp;∼&nbsp;LogNormal(log&nbsp;Σ,&nbsp;τ<sub>σ</sub>) | M&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>M</sub>) | log&nbsp;Σ&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>Σ</sub>) | all τ&nbsp;∼&nbsp;half‑Cauchy(0,&nbsp;2.5) | γ<sub>i</sub>&nbsp;∼&nbsp;U(0,&nbsp;1), λ<sub>i</sub>&nbsp;∼&nbsp;U(0,&nbsp;1&nbsp;−&nbsp;γ<sub>i</sub>) |
| BH-GNM, half-normal | `Stan/example1_BHGNM_halfnormal.stan` |  | μ<sub>i</sub>&nbsp;∼&nbsp;N(M,&nbsp;τ<sub>μ</sub>) | σ<sub>i</sub>&nbsp;∼&nbsp;N(Σ,&nbsp;τ<sub>σ</sub>), σ<sub>i</sub>&nbsp;>&nbsp;0 | M&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>M</sub>) | Σ&nbsp;∼&nbsp;half‑N(0,&nbsp;τ<sub>Σ</sub>) | τ&nbsp;∼&nbsp;half‑Cauchy(0,&nbsp;2.5); τ<sub>γ</sub>, τ<sub>λ</sub>&nbsp;∼&nbsp;half‑Cauchy(0,&nbsp;0.5) | γ<sub>i</sub>&nbsp;∼&nbsp;N(Γ,&nbsp;τ<sub>γ</sub>) on [0,&nbsp;1], Γ&nbsp;∼&nbsp;N(0.05,&nbsp;0.1) on [0,&nbsp;1]; same for λ |
| BH-GNM (brms) | `R/example1_BHGLM_BHGNM_brms.R` |  | μ<sub>i</sub>, normal | σ<sub>i</sub>, normal | μ&nbsp;∼&nbsp;N(75,&nbsp;50) | σ&nbsp;∼&nbsp;half‑N(0,&nbsp;50) | sd&nbsp;∼&nbsp;half‑Cauchy(0,&nbsp;2.5) | γ<sub>i</sub>, λ<sub>i</sub>&nbsp;∼&nbsp;U(0,&nbsp;1) |

### Example 2: experimental data (`vibro_exp3`)

| Model | File | Paper | Subject intercept | Subject slope | Population intercept | Population slope | Standard deviations | Guess γ / lapse λ |
|---|---|:---:|---|---|---|---|---|---|
| BH-GLM | `Stan/example2_BHGLM.stan` | ✓ | b<sub>0ih</sub>&nbsp;∼&nbsp;N(β<sub>0h</sub>,&nbsp;τ<sub>b0</sub>) | b<sub>1ih</sub>&nbsp;∼&nbsp;N(β<sub>1h</sub>,&nbsp;τ<sub>b1</sub>) | β<sub>0h</sub>&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>β0</sub>) | β<sub>1h</sub>&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>β1</sub>) | all τ&nbsp;∼&nbsp;half‑Cauchy(0,&nbsp;2.5) | – |
| BH-GNM | `Stan/example2_BHGNM.stan` | ✓ | as BH-GLM | as BH-GLM | as BH-GLM | as BH-GLM | as BH-GLM | γ<sub>ih</sub>&nbsp;∼&nbsp;U(0,&nbsp;1), λ<sub>ih</sub>&nbsp;∼&nbsp;U(0,&nbsp;1&nbsp;−&nbsp;γ<sub>ih</sub>) |
| BH-GLM, lognormal | `Stan/example2_BHGLM_lognormal.stan` |  | b<sub>0ih</sub>&nbsp;∼&nbsp;N(β<sub>0h</sub>,&nbsp;τ<sub>b0</sub>) | b<sub>1ih</sub>&nbsp;∼&nbsp;LogNormal(log&nbsp;β<sub>1h</sub>,&nbsp;τ<sub>b1</sub>), non-centred | β<sub>0h</sub>&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>β0</sub>) | log&nbsp;β<sub>1h</sub>&nbsp;∼&nbsp;N(0,&nbsp;τ<sub>β1</sub>) | τ<sub>b0</sub>, τ<sub>β0</sub>&nbsp;∼&nbsp;half‑N(0,&nbsp;2.5); τ<sub>b1</sub>, τ<sub>β1</sub>&nbsp;∼&nbsp;half‑N(0,&nbsp;1) | – |

---

## Reference

Dallmann, C. J., Ernst, M. O., & Moscatelli, A. (2015). The role of vibration in tactile speed perception. *Journal of Neurophysiology*, 114(6), 3131–3139.
