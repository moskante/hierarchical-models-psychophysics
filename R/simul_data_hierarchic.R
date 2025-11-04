library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(DHARMa)
library(rstan)
library(shinystan)

set.seed(123)
simul_data <- PsySimulate(ntrials = 160, nsubjects = 10, guess = TRUE, lapse = TRUE)

# GLMM ----
glmm <- glmer(cbind(Longer, Total-Longer) ~ X + (1+X|Subject), family = binomial(link = "probit"), data = simul_data)

## residuals ---
residuals_glmm <- residuals(glmm, type = "response")
SSE_glmm <- sum(residuals_glmm^2)
#loglik_glmm <- as.numeric(logLik(glmm))


## DHARMa diagnostic -----
testDispersion(glmm)
simulationOutput <- simulateResiduals(fittedModel = glmm)
plot(simulationOutput)
testUniformity(simulationOutput)
testDispersion(simulationOutput)
testQuantiles(simulationOutput)
testOutliers(simulationOutput)

## bootstrap -----
glmm_CI <- pseMer(glmm, B = 500)

# BHMs data ----
datistan=list()
datistan$y=simul_data$Longer
datistan$n=simul_data$Total
datistan$nobs=length(datistan$y)
datistan$x=simul_data$X
datistan$subject=simul_data$Subject
datistan$nsubj=length(unique(datistan$subject))
datistan$ncond=length(unique(datistan$condition))


# BH-GLM ----
init_fun <- function(init_pse = 90, init_beta = 15) {
  list(
    pse = rep(init_pse, datistan$nsubj),
    beta  = rep(init_beta, datistan$nsubj))}

fit_bhglm <- stan(
  file = "Stan/simul_bhglm.stan",  
  data = datistan,    
  chains = 3,             
  warmup = 3000,          
  iter = 5000,         
  cores = 3,              
  refresh = 10,             
  init=list(init_fun(),init_fun(),init_fun()))  

## Residuals ----
samples_bhglm <- extract(fit_bhglm)
posteriorPredDistr <- samples_bhglm$predProb
posteriorPredSim   <- samples_bhglm$y_sim  

fitted_probs <- apply(posteriorPredDistr, 2, median)
y_obs <- datistan$y
trials <- datistan$n
SSE_bhglm <- sum( (y_obs/trials - fitted_probs)^2 )

## DHARMa diagnostic ------
sim_bhglm <- createDHARMa(
  simulatedResponse = t(posteriorPredSim),
  observedResponse = datistan$y,
  fittedPredictedResponse = apply(posteriorPredDistr, 2, median),
  integerResponse = TRUE
)

plot(sim_bhglm)
testUniformity(sim_bhglm)
testDispersion(sim_bhglm)
testOutliers(sim_bhglm)
testQuantiles(sim_bhglm)


# BH-GNM ----
init_fun <- function(init_pse = 90, init_beta = 15,
                     init_gamma = 0.01, init_lambda = 0.01) {
  list(
    pse = rep(init_pse, datistan$nsubj),
    beta  = rep(init_beta, datistan$nsubj),
    gamma  = rep(init_gamma, datistan$nsubj),
    lambda  = rep(init_lambda, datistan$nsubj)
  )
}
  
  
fit_bhgnm <- stan(
  file = "Stan/simul_bhgnm.stan",  
  data = datistan,    
  chains = 3,             
  warmup = 3000,          
  iter = 5000,          
  cores = 3,              
  refresh = 10,             
  init=list(init_fun(),init_fun(),init_fun())
  # control = list(adapt_delta = 0.999, stepsize = 0.01, max_treedepth = 15)
)  

## Residuals ----
samples_bhgnm <- extract(fit_bhgnm)
posteriorPredDistr <- samples_bhgnm$predProb
posteriorPredSim   <- samples_bhgnm$y_sim  

fitted_probs <- apply(posteriorPredDistr, 2, median)
y_obs <- datistan$y
trials <- datistan$n
SSE_bhgnm <- sum( (y_obs/trials - fitted_probs)^2 )

## DHARMa diagnostic ------
sim_bhgnm <- createDHARMa(
  simulatedResponse = t(posteriorPredSim),
  observedResponse = datistan$y,
  fittedPredictedResponse = apply(posteriorPredDistr, 2, median),
  integerResponse = TRUE
)

testUniformity(sim_bhgnm)
testDispersion(sim_bhgnm)
testOutliers(sim_bhgnm)
testQuantiles(sim_bhgnm)

