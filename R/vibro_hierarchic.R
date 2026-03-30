# vibro_hierarchic.R
#
# Example 2b: hierarchical model analysis on vibrotactile data.
#
# Workflow:
#   1. Fit a probit GLMM with a speed x vibration interaction and random slopes.
#   2. Bootstrap CIs for PSE, JND, slope, and condition differences.
#   3. Fit BH-GLM and BH-GNM via Stan.
#   4. Validate each model with DHARMa diagnostics.
#
# Required packages: tidyverse, MixedPsy, lme4, lmerTest, DHARMa, rstan, shinystan
# Stan models:       Stan/vibro_bhglm.stan, Stan/vibro_bhgnm.stan
# Dataset:           vibro_exp3 (from MixedPsy)

library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(DHARMa)
library(rstan)
library(shinystan)

# vibro_exp3 is bundled with the MixedPsy package.

# ---------------------------------------------------------------------------
# GLMM: probit mixed model with speed x vibration interaction
#
# Fixed effects:
#   speed            -- main effect of speed (common slope, 0 Hz condition)
#   vibration        -- shift in intercept at 32 Hz
#   speed:vibration  -- shift in slope at 32 Hz (key test of the vibration effect)
# Random effects:
#   (1 + speed | subject) -- subjects differ in their baseline PSE and sensitivity
# ---------------------------------------------------------------------------

# GLMM -----
glmm_vibro <- glmer(cbind(faster, slower) ~ speed * vibration + (1 + speed | subject),
                    family = binomial(link = "probit"),
                    data   = vibro_exp3)

# Residuals and SSE for model comparison
## Residuals -----
residuals_glmm_vibro <- residuals(glmm_vibro, type = "response")
SSE_glmm_vibro       <- sum(residuals_glmm_vibro^2)

# DHARMa diagnostics for GLMM
simulatedRes_glmm_vibro <- simulateResiduals(fittedModel = glmm_vibro)
testUniformity(simulatedRes_glmm_vibro)
testDispersion(simulatedRes_glmm_vibro)
testOutliers(simulatedRes_glmm_vibro)
testQuantiles(simulatedRes_glmm_vibro)

# ---------------------------------------------------------------------------
# Bootstrap CIs from the GLMM (500 replicates)
#
# fun2mod() is a custom extraction function passed to pseMer(). It extracts:
#   PSE and JND for both vibration conditions (0 Hz and 32 Hz)
#   Differences in PSE, JND, and slope across conditions
#
# The GLMM is parameterised as:
#   Intercept       = beta_0 (baseline intercept, 0 Hz condition)
#   speed           = beta_1 (baseline slope, 0 Hz)
#   vibration32     = beta_2 (intercept shift at 32 Hz)
#   speed:vibration = beta_4 (slope shift at 32 Hz)
#
#   PSE_0  = -beta_0 / beta_1
#   PSE_32 = -(beta_0 + beta_2) / (beta_1 + beta_4)
#   JND    = qnorm(0.75) / slope
# ---------------------------------------------------------------------------
## bootstrap -----
fun2mod <- function(mer.obj) {
  jndpse        <- vector(mode = "numeric", length = 7)
  names(jndpse) <- c("PSE 0", "PSE 32", "JND 0", "JND 32",
                     "delta PSE", "delta JND", "delta slope")
  jndpse[1] <- -fixef(mer.obj)[1] / fixef(mer.obj)[2]                                  # PSE (0 Hz)
  jndpse[2] <- -(fixef(mer.obj)[1] + fixef(mer.obj)[3]) /
               (fixef(mer.obj)[2] + fixef(mer.obj)[4])                                  # PSE (32 Hz)
  jndpse[3] <- qnorm(0.75) / fixef(mer.obj)[2]                                         # JND (0 Hz)
  jndpse[4] <- qnorm(0.75) / (fixef(mer.obj)[2] + fixef(mer.obj)[4])                  # JND (32 Hz)
  jndpse[5] <- jndpse[2] - jndpse[1]   # Difference in PSE (32 Hz - 0 Hz)
  jndpse[6] <- jndpse[4] - jndpse[3]   # Difference in JND (32 Hz - 0 Hz)
  jndpse[7] <- fixef(mer.obj)[4]        # Raw slope difference (interaction coefficient)
  return(jndpse)
}

boot_vibro <- pseMer(glmm_vibro, B = 500, FUN = fun2mod)

# ---------------------------------------------------------------------------
# Prepare data list for Stan
#
# vibration is re-coded as 0/1 integer (vibr). Subject labels are mapped
# to consecutive integers (soggetto) because Stan requires integer indices.
# The structure datistan$x, datistan$subject, and datistan$vibration must
# align row-by-row with datistan$y (number of "faster" responses).
# ---------------------------------------------------------------------------
# BHMs data -----
attach(vibro_exp3)
# Encode vibration condition as 0/1
vibr <- 0 * speed
vibr[vibration == "32"] <- 1
# Encode subject labels as consecutive integers
soggetto <- 0 * speed
soggetto[subject == "AK"] <- 1
soggetto[subject == "AR"] <- 2
soggetto[subject == "DN"] <- 3
soggetto[subject == "FA"] <- 4
soggetto[subject == "MA"] <- 5
soggetto[subject == "MI"] <- 6
soggetto[subject == "NI"] <- 7
soggetto[subject == "NN"] <- 8
soggetto[subject == "RV"] <- 9

datistan          <- list()
datistan$y        <- faster
datistan$n        <- faster + slower
datistan$nobs     <- length(datistan$y)
datistan$x        <- speed
datistan$vibration <- vibr
datistan$subject  <- soggetto
datistan$nsubj    <- length(unique(datistan$subject))
datistan$ncond    <- length(unique(datistan$vibration))
detach(vibro_exp3)

# ---------------------------------------------------------------------------
# BH-GLM: Bayesian Hierarchical GLM (two-condition version)
#
# The Stan model (vibro_bhglm.stan) uses a matrix parameterisation:
#   b0[s, v] and b1[s, v] are per-subject intercepts and slopes for
#   condition v (1 = 0 Hz, 2 = 32 Hz).
# Initialisation near zero / small positive values for the linear predictor
# keeps the probit function in a reasonable range at the start of sampling.
# ---------------------------------------------------------------------------

# BH-GLM ----
init_fun <- function(init_b0 = -2, init_b1 = 0.1) {
  list(
    b0 = matrix(init_b0, datistan$nsubj, 2),
    b1 = matrix(init_b1, datistan$nsubj, 2)
  )
}

fit_bhglm_vibro <- stan(
  file   = "Stan/vibro_bhglm.stan",
  data   = datistan,
  chains = 3,
  warmup = 3000,
  iter   = 5000,
  cores  = 3,
  refresh = 10,
  init   = list(init_fun(), init_fun(), init_fun())
)

## Residuals -----
samples_bhglm_vibro <- extract(fit_bhglm_vibro)
posteriorPredDistr  <- samples_bhglm_vibro$predProb
posteriorPredSim    <- samples_bhglm_vibro$y_sim

fitted_probs     <- apply(posteriorPredDistr, 2, median)
y_obs   <- datistan$y
trials  <- datistan$n
SSE_bhglm_vibro  <- sum((y_obs / trials - fitted_probs)^2)

## DHARMa diagnostic ------
sim_bhglm_vibro <- createDHARMa(
  simulatedResponse        = t(posteriorPredSim),
  observedResponse         = datistan$y,
  fittedPredictedResponse  = apply(posteriorPredDistr, 2, median),
  integerResponse          = TRUE
)

testUniformity(sim_bhglm_vibro)
testDispersion(sim_bhglm_vibro)
testOutliers(sim_bhglm_vibro)
testQuantiles(sim_bhglm_vibro)

# ---------------------------------------------------------------------------
# BH-GNM: Bayesian Hierarchical GNM (two-condition version)
#
# Extends BH-GLM with condition-specific guessing and lapsing matrices:
#   gamma[s, v] and lambda[s, v].
# ---------------------------------------------------------------------------

# BH-GNM ------
init_fun <- function(init_b0 = -2, init_b1 = 0.1,
                     init_gamma = 0.01, init_lambda = 0.01) {
  list(
    b0     = matrix(init_b0,     datistan$nsubj, 2),
    b1     = matrix(init_b1,     datistan$nsubj, 2),
    gamma  = matrix(init_gamma,  datistan$nsubj, 2),
    lambda = matrix(init_lambda, datistan$nsubj, 2)
  )
}

fit_bhgnm_vibro <- stan(
  file   = "Stan/vibro_bhgnm.stan",
  data   = datistan,
  chains = 3,
  warmup = 3000,
  iter   = 5000,
  cores  = 3,
  refresh = 10,
  init   = list(init_fun(), init_fun(), init_fun())
)

## Residuals -----
samples_bhgnm_vibro <- extract(fit_bhgnm_vibro)
posteriorPredDistr  <- samples_bhgnm_vibro$predProb
posteriorPredSim    <- samples_bhgnm_vibro$y_sim

fitted_probs     <- apply(posteriorPredDistr, 2, median)
y_obs   <- datistan$y
trials  <- datistan$n
SSE_bhgnm_vibro  <- sum((y_obs / trials - fitted_probs)^2)

## DHARMa diagnostic ------
sim_bhgnm_vibro <- createDHARMa(
  simulatedResponse        = t(posteriorPredSim),
  observedResponse         = datistan$y,
  fittedPredictedResponse  = apply(posteriorPredDistr, 2, median),
  integerResponse          = TRUE
)

testUniformity(sim_bhgnm_vibro)
testDispersion(sim_bhgnm_vibro)
testOutliers(sim_bhgnm_vibro)
testQuantiles(sim_bhgnm_vibro)
