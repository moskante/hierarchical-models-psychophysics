# simul_data_single_sub.R
#
# Example 1a: single-subject GLM and GNM analysis on simulated psychophysical data.
#
# Workflow:
#   1. Simulate 10 subjects with known population parameters.
#   2. Extract the true PSE and JND for each subject (used as reference values).
#   3. Fit a probit GLM per subject (standard approach, no guessing/lapsing).
#   4. Fit a GNM per subject (extends GLM with free guessing and lapsing parameters).
#   5. Test whether estimated PSE and JND deviate from the true values.
#
# Required packages: tidyverse, MixedPsy, gnlm
# Source file:       R/gnlm_functions_psejnd.R

library(tidyverse)
library(MixedPsy)
library(gnlm)

# ---------------------------------------------------------------------------
# 1. Simulate dataset
#
# PsySimulate() generates binomial psychophysical data for nsubjects subjects.
# Each subject completes ntrials trials; the true psychometric function
# includes guessing (gamma) and lapsing (lambda) rates drawn from the
# population distributions specified in the MixedPsy package documentation.
# The Subject column is factored so that downstream functions treat it as a
# grouping variable.
# ---------------------------------------------------------------------------
set.seed(123)
simul_data <- PsySimulate(ntrials = 160, nsubjects = 10, guess = TRUE, lapse = TRUE)
simul_data$Subject <- factor(simul_data$Subject)

# ---------------------------------------------------------------------------
# 2. Extract true population parameters
#
# PsySimulate() stores the generating intercept, slope, gamma, and lambda for
# each subject as repeated values within each subject's rows. We take the
# first row per subject to recover these true values, then derive:
#   PSE = -Intercept / Slope   (stimulus level where P(response) = 0.5)
#   JND = qnorm(0.75) / Slope  (spread at the 75th quantile of the probit)
# The mean across subjects becomes the reference value for inference.
# ---------------------------------------------------------------------------
parameters_simul <- simul_data %>%
  group_by(Subject) %>%
  summarise(across(everything(), first),
            sigma = 1/Slope,
            gamma  = Gamma,
            lambda = Lambda,
            k = qnorm((0.5 - gamma) / (1 - gamma - lambda)) * sigma,
            pse    = -Intercept / Slope + k,
            jnd    = qnorm((0.75 - gamma) / (1 - gamma - lambda)) * sigma - k
            ) %>%
  select(Subject, pse, jnd, gamma, lambda)

sample_mean_pse <- mean(parameters_simul$pse)   # sample mean PSE (true reference)
sample_mean_jnd <- mean(parameters_simul$jnd)   # sample mean JND (true reference)

# ---------------------------------------------------------------------------
# 3. GLM: probit model fit for each subject
#
# PsychModels() fits separate probit GLMs to each subject's data. The
# response is a two-column matrix (successes, failures).
# PsychParameters() extracts PSE and JND from the list of fitted models.
# ---------------------------------------------------------------------------

# GLM ----
glm_list <- PsychModels(formula        = cbind(Longer, Total - Longer) ~ X,
                        data           = simul_data,
                        group_factors  = "Subject")

parameters_glm <- PsychParameters(glm_list)

# Inference: does the GLM recover the true population PSE and JND?
# One-sample t-tests compare estimated subject-level values against the
# known generating mean.
## inference tests -----
t_glm_pse <- t.test(parameters_glm$pse, mu = sample_mean_pse)
t_glm_jnd <- t.test(parameters_glm$jnd, mu = sample_mean_jnd)


# ---------------------------------------------------------------------------
# 4. GNM: nonlinear model with guessing and lapsing
#
# gnlm_functions_psejnd.R defines process_subject(), which fits a GNM to one
# subject's data and returns PSE, JND, gamma, lambda, and their bootstrap SEs.
# group_split() splits the data frame by Subject; map_dfr() applies
# process_subject() to each split and stacks the results into one tibble.
# ---------------------------------------------------------------------------

# GNL ----
source("R/gnlm_functions_psejnd.R")
parameters_gnm <- simul_data %>%
  group_split(Subject) %>%
  map_dfr(process_subject)

# Inference: does the GNM recover the true PSE and JND?
## inference tests -----
t_gnm_pse <- t.test(parameters_gnm$pse, mu = sample_mean_pse)
t_gnm_jnd <- t.test(parameters_gnm$jnd, mu = sample_mean_jnd)

