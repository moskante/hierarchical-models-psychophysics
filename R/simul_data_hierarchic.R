# simul_data_hierarchic.R
#
# Example 1b: hierarchical model analysis on simulated psychophysical data.
#
# Workflow:
#   1. Simulate 10 subjects (same setup as simul_data_single_sub.R).
#   2. Fit a probit GLMM with random intercepts and slopes (frequentist).
#   3. Validate the GLMM with DHARMa residual diagnostics.
#   4. Fit a Bayesian Hierarchical GLM (BH-GLM) via Stan.
#   5. Fit a Bayesian Hierarchical GNM (BH-GNM) via Stan (adds guessing/lapsing).
#   6. Validate each Bayesian model with DHARMa diagnostics using posterior
#      predictive simulations.
#
# Required packages: tidyverse, MixedPsy, lme4, lmerTest, DHARMa, rstan, shinystan
# Stan models:       Stan/simul_bhglm.stan, Stan/simul_bhgnm.stan

library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(DHARMa)
library(rstan)
library(shinystan)

# Simulate the same dataset structure as in Example 1a
set.seed(123)
simul_data <- PsySimulate(ntrials = 160, nsubjects = 10, guess = TRUE, lapse = TRUE)

# ---------------------------------------------------------------------------
# GLMM: probit mixed model with random intercepts and slopes
#
# The fixed effect of X estimates the population-level psychometric slope.
# (1 + X | Subject) allows both the intercept (related to PSE) and the slope
# to vary randomly across subjects, capturing individual differences.
# ---------------------------------------------------------------------------

# GLMM ----
glmm <- glmer(cbind(Longer, Total - Longer) ~ X + (1 + X | Subject),
              family = binomial(link = "probit"),
              data   = simul_data)

# Residuals and SSE: response-scale residuals (observed proportion - fitted probability)
## residuals ---
residuals_glmm <- residuals(glmm, type = "response")
SSE_glmm       <- sum(residuals_glmm^2)

# Bootstrap CI for PSE and JND from the GLMM fixed effects (500 resamples)
## bootstrap -----
glmm_CI <- pseMer(glmm, B = 500)

# ---------------------------------------------------------------------------
# DHARMa diagnostics for the GLMM
#
# DHARMa uses simulation-based scaled residuals to assess model assumptions.
# testUniformity checks whether residuals are uniformly distributed (as
# expected for a correctly specified model). testDispersion checks for
# over- or under-dispersion. testQuantiles tests whether quantile predictions
# match observations. testOutliers flags extreme residuals.
# ---------------------------------------------------------------------------
## DHARMa diagnostic -----
testDispersion(glmm)
simulationOutput <- simulateResiduals(fittedModel = glmm)
plot(simulationOutput)
testUniformity(simulationOutput)
testDispersion(simulationOutput)
testQuantiles(simulationOutput)
testOutliers(simulationOutput)

# ---------------------------------------------------------------------------
# Prepare data list for Stan
#
# Stan requires a named list. Each element corresponds to a variable declared
# in the data{} block of the .stan file.
#   y       -- number of "longer" responses per trial row (binomial successes)
#   n       -- total trials per row (binomial denominator)
#   x       -- stimulus level (continuous predictor)
#   subject -- integer-coded subject ID (1 to nsubj)
#   nsubj   -- number of unique subjects
# ---------------------------------------------------------------------------
# BHMs data ----
datistan        <- list()
datistan$y      <- simul_data$Longer
datistan$n      <- simul_data$Total
datistan$nobs   <- length(datistan$y)
datistan$x      <- simul_data$X
datistan$subject <- simul_data$Subject
datistan$nsubj  <- length(unique(datistan$subject))
datistan$ncond  <- length(unique(datistan$condition))

# ---------------------------------------------------------------------------
# BH-GLM: Bayesian Hierarchical GLM
#
# The Stan model (simul_bhglm.stan) places hierarchical priors on per-subject
# PSE and sigma (slope), with Cauchy(0, 2.5) hyperpriors on the population
# standard deviations.
#
# Initialisation: starting the chains near plausible values (PSE ~ 90, sigma ~ 15)
# improves mixing and reduces the warmup needed.
# Chains: 3 chains, 3000 warmup + 2000 sampling iterations each.
# ---------------------------------------------------------------------------

# BH-GLM ----
init_fun <- function(init_pse = 90, init_sigma = 15) {
  list(
    pse   = rep(init_pse,   datistan$nsubj),
    sigma = rep(init_sigma, datistan$nsubj)
  )
}

fit_bhglm <- stan(
  file   = "Stan/simul_bhglm.stan",
  data   = datistan,
  chains = 3,
  warmup = 3000,
  iter   = 5000,   # total iterations per chain (including warmup)
  cores  = 3,
  refresh = 10,
  init   = list(init_fun(), init_fun(), init_fun())
)

# Extract posterior predictive distributions for residual diagnostics
## Residuals ----
samples_bhglm      <- extract(fit_bhglm)
posteriorPredDistr <- samples_bhglm$predProb   # posterior predictive probabilities
posteriorPredSim   <- samples_bhglm$y_sim      # posterior predictive counts (for DHARMa)

# SSE: compare median posterior prediction to observed proportions
fitted_probs  <- apply(posteriorPredDistr, 2, median)
y_obs  <- datistan$y
trials <- datistan$n
SSE_bhglm <- sum((y_obs / trials - fitted_probs)^2)

# DHARMa diagnostics using posterior predictive simulations
# posteriorPredSim is a (draws x observations) matrix; DHARMa expects
# (observations x draws), hence the transpose.
## DHARMa diagnostic ------
sim_bhglm <- createDHARMa(
  simulatedResponse        = t(posteriorPredSim),
  observedResponse         = datistan$y,
  fittedPredictedResponse  = apply(posteriorPredDistr, 2, median),
  integerResponse          = TRUE
)

plot(sim_bhglm)
testUniformity(sim_bhglm)
testDispersion(sim_bhglm)
testOutliers(sim_bhglm)
testQuantiles(sim_bhglm)


# ---------------------------------------------------------------------------
# BH-GNM: Bayesian Hierarchical GNM
#
# Extends BH-GLM by estimating per-subject guessing (gamma) and lapsing
# (lambda) rates with uniform priors. All other hyperparameters are shared
# with the BH-GLM specification.
# ---------------------------------------------------------------------------

# BH-GNM ----
# init_fun <- function(init_pse = 90, init_sigma = 15,
#                      init_gamma = 0.01, init_lambda = 0.01) {
#   list(
#     pse    = rep(init_pse,    datistan$nsubj),
#     sigma  = rep(init_sigma,  datistan$nsubj),
#     gamma  = rep(init_gamma,  datistan$nsubj),
#     lambda = rep(init_lambda, datistan$nsubj)
#   )
# }

init_fun <- function(init_mu = 90, init_sigma = 15,
                     init_gamma = 0.01, init_lambda_raw = 0.01) {
  list(
    # Subject-level parameters
    mu         = rep(init_mu, datistan$nsubj),  # mu = pse when gamma == lambda
    sigma      = rep(init_sigma, datistan$nsubj),
    gamma      = rep(init_gamma, datistan$nsubj),
    lambda_raw = rep(init_lambda_raw, datistan$nsubj),
    
    # Population-level hyper-parameters (Good practice to initialize as well)
    MU           = init_mu,
    SIGMA        = init_sigma,
    GAMMA        = init_gamma,
    LAMBDA_tilde = init_lambda_raw
  )
}


fit_bhgnm <- stan(
  file   = "Stan/simul_bhgnm_GL.stan",
  data   = datistan,
  chains = 3,
  warmup = 3000,
  iter   = 5000,
  cores  = 3,
  refresh = 10,
  init   = list(init_fun(), init_fun(), init_fun())
  # control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15)
)

## Residuals ----
samples_bhgnm      <- extract(fit_bhgnm)
posteriorPredDistr <- samples_bhgnm$predProb
posteriorPredSim   <- samples_bhgnm$y_sim

fitted_probs  <- apply(posteriorPredDistr, 2, median)
y_obs  <- datistan$y
trials <- datistan$n
SSE_bhgnm <- sum((y_obs / trials - fitted_probs)^2)

## DHARMa diagnostic ------
sim_bhgnm <- createDHARMa(
  simulatedResponse        = t(posteriorPredSim),
  observedResponse         = datistan$y,
  fittedPredictedResponse  = apply(posteriorPredDistr, 2, median),
  integerResponse          = TRUE
)

testUniformity(sim_bhgnm)
testDispersion(sim_bhgnm)
testOutliers(sim_bhgnm)
testQuantiles(sim_bhgnm)
