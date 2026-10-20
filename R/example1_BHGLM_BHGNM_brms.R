# ==============================================================================
# Bayesian hierarchical models with brms (alternative to native Stan)
# ==============================================================================

library(tidyverse)
library(MixedPsy)
library(brms)
library(DHARMa)

# 1. Data preparation ----
set.seed(123)
simul_data <- PsySimulate(ntrials = 160, nsubjects = 10, guess = TRUE, lapse = TRUE)
simul_data$Subject <- factor(simul_data$Subject)

# 2. Model 1: BH-GLM (no guess/lapse rates) ----
# brms formula syntax is the same as lme4

# 3. Fit the model with tuned priors and initial values
fit_bhglm <- brm(
  formula = Longer | trials(Total) ~ X + (1 + X | Subject),
  data = simul_data,
  family = binomial(link = "probit"), 
  prior = c(
    prior(normal(0, 3), class = "Intercept"), # weakly informative priors on the probit scale
    prior(normal(0, 3), class = "b"),
    prior(cauchy(0, 2.5), class = "sd"),
    prior(lkj(2), class = "cor")              # prior on the random-effects correlation
  ),
  chains = 3, 
  iter = 5000, 
  warmup = 3000,
  refresh = 10, 
  init = 0,               # initialise all parameters at 0 (unconstrained scale) for stable sampling
  backend = "rstan"       # or "cmdstanr" if installed
)

# 4. Check the fitted object
if (is(fit_bhglm, "brmsfit")) {
  summary(fit_bhglm)
} else {
  print("Il modello non è stato creato correttamente.")
}


## ---  Predictive Distributions with brms ---
# expected counts (n * pi)
pred_counts <- posterior_epred(fit_bhglm) 

# divide by number of trials to get probabilities (pi)
trials_vec <- simul_data$Total
posteriorPredDistr_brms <- sweep(pred_counts, 2, trials_vec, "/")

posteriorPredSim_brms <- posterior_predict(fit_bhglm)

# fitted probabilities and SSE
fitted_probs <- apply(posteriorPredDistr_brms, 2, median)

# observed proportions from the data frame
y_obs <- simul_data$Longer
trials <- simul_data$Total

# Sum of Squared Errors (SSE)
SSE_bhglm <- sum( (y_obs/trials - fitted_probs)^2 )

print(paste("SSE bhglm:", SSE_bhglm))


# create DHARMa residuals object
sim_bhglm <- createDHARMa(
  simulatedResponse = t(posteriorPredSim_brms), 
  observedResponse = y_obs,
  fittedPredictedResponse = fitted_probs,
  integerResponse = TRUE
)

plot(sim_bhglm)
testUniformity(sim_bhglm)
testDispersion(sim_bhglm)
testOutliers(sim_bhglm)

extract_psychophysic_params <- function(model) {
  # 1. Extract fixed effects
  fe <- fixef(model)
  b0 <- fe["Intercept", "Estimate"]
  b1 <- fe["X", "Estimate"]
  
  # 2. Extract subject-level random effects
  re <- ranef(model)$Subject
  
  # 3. Population-level parameters
  # PSE = -b0 / b1
  # Sigma = 1 / b1
  pse_group <- -b0 / b1
  sigma_group <- 1 / b1
  jnd_group <- 0.6745 * sigma_group
  
  # 4. Subject-level parameters
  # fixed effect + subject-specific random effect
  subjects <- rownames(re)
  subj_results <- data.frame(
    Subject = subjects,
    PSE = numeric(length(subjects)),
    Sigma = numeric(length(subjects)),
    JND = numeric(length(subjects))
  )
  
  for (i in seq_along(subjects)) {
    # subject-specific intercept and slope
    s_b0 <- b0 + re[i, "Estimate", "Intercept"]
    s_b1 <- b1 + re[i, "Estimate", "X"]
    
    subj_results$PSE[i]   <- -s_b0 / s_b1
    subj_results$Sigma[i] <- 1 / s_b1
    subj_results$JND[i]   <- 0.6745 * (1 / s_b1)
  }
  
  return(list(
    Population = data.frame(PSE = pse_group, Sigma = sigma_group, JND = jnd_group),
    Subjects = subj_results
  ))
}


params <- extract_psychophysic_params(fit_bhglm)
print(params$Population)
print(params$Subjects)

# 5. Model 2: BH-GNM (with guess/lapse rates) ----
# Nonlinear model including the lower and upper asymptotes (approximates Stan/example1_BHGNM.stan).

# Psychometric function: gamma + (1 - gamma - lambda) * Phi((X - mu) / sigma)

formula_bhgnm <- brms::bf(
  Longer | trials(Total) ~ gammaPar + (1 - gammaPar - lambdaPar) * Phi((X - muPar) / sigmaPar), # mu and sigma are Reserved Names in brms
  muPar ~ 1 + (1 | Subject),
  sigmaPar ~ 1 + (1 | Subject),
  
  gammaPar ~ 0 + Subject,
  lambdaPar ~ 0 + Subject,
  
  nl = TRUE
)

# Alternative closer to the Stan model (log-scale sigma):
#
# formula_bhgnm <- brms::bf(
#   Longer | trials(Total) ~ gammaPar + (1 - gammaPar - lambdaPar) * Phi((X - muPar) / exp(logSigmaPar)),
#   muPar ~ 1 + (1 | Subject),
#   logSigmaPar ~ 1 + (1 | Subject), # exp(logSigmaPar) guarantees positivity
#   gammaPar ~ 0 + Subject,
#   lambdaPar ~ 0 + Subject,
#   nl = TRUE
# )

# Priors
# lb/ub bounds keep sigma positive and gamma, lambda in [0, 1]
priors_bhgnm <- c(
  # population-level location and scale, on the stimulus scale (X from 40 to 120)
  prior(normal(75, 50), class = "b", nlpar = "muPar"),
  prior(normal(0, 50), class = "b", nlpar = "sigmaPar", lb = 0),
  
  # group-level SDs (equivalent to tau_pse and tau_sigma)
  prior(cauchy(0, 2.5), class = "sd", nlpar = "muPar"),
  prior(cauchy(0, 2.5), class = "sd", nlpar = "sigmaPar"),
  
  # priors on gamma and lambda (bounded in [0, 1])
  prior(uniform(0, 1), class = "b", nlpar = "gammaPar", lb = 0, ub = 1),
  prior(uniform(0, 1), class = "b", nlpar = "lambdaPar", lb = 0, ub = 1)
)

# Fit the model
fit_bhgnm_brms <- brm(
  formula = formula_bhgnm,
  data = simul_data, 
  family = binomial(link = "identity"), 
  prior = priors_bhgnm,
  chains = 3,
  iter =5000,
  warmup = 3000,
  cores = 3,
  init = 0 
)



## ---  Predictive Distributions with brms ---
# expected counts (n * pi)
pred_counts <- posterior_epred(fit_bhgnm_brms) 

# divide by number of trials to get probabilities (pi)
trials_vec <- simul_data$Total
posteriorPredDistr_brms_bhgnm <- sweep(pred_counts, 2, trials_vec, "/")

posteriorPredSim_brms_bhgnm <- posterior_predict(fit_bhgnm_brms)

# fitted probabilities and SSE
fitted_probs <- apply(posteriorPredDistr_brms_bhgnm, 2, median)

# observed proportions from the data frame
y_obs <- simul_data$Longer
trials <- simul_data$Total

# Sum of Squared Errors (SSE)
SSE_bhgnm <- sum( (y_obs/trials - fitted_probs)^2 )

print(paste("SSE bhgnm:", SSE_bhgnm))



# create DHARMa residuals object
sim_bhgnm <- createDHARMa(
  simulatedResponse = t(posteriorPredSim_brms_bhgnm), 
  observedResponse = y_obs,
  fittedPredictedResponse = fitted_probs,
  integerResponse = TRUE
)

plot(sim_bhgnm)
testUniformity(sim_bhgnm)
testDispersion(sim_bhgnm)
testOutliers(sim_bhgnm)
