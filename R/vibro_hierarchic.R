library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(DHARMa)
library(rstan)
library(shinystan)

# data are already in MixedPsy package. Alternatively can be loaded with:
# set path befor
# load("vibro_exp3.RData")

# GLMM -----
glmm_vibro <- glmer(cbind(faster, slower) ~ speed * vibration + (1+speed | subject),
                    family = binomial(link = "probit"),
                    data = vibro_exp3)

## Residuals -----
residuals_glmm_vibro <- residuals(glmm_vibro, type = "response")
SSE_glmm_vibro <- sum(residuals_glmm_vibro^2)
loglik_glmm_vibro <- as.numeric(logLik(glmm_vibro))

simulatedRes_glmm_vibro <- simulateResiduals(fittedModel = glmm_vibro)

testUniformity(simulatedRes_glmm_vibro)
testDispersion(simulatedRes_glmm_vibro)
testOutliers(simulatedRes_glmm_vibro)
testQuantiles(simulatedRes_glmm_vibro)

## bootstrap -----
fun2mod <- function(mer.obj){
  jndpse = vector(mode = "numeric", length = 7)
  names(jndpse) = c("PSE 0", "PSE 32","JND 0", "JND 32", "delta PSE", "delta JND", "delta slope")
  jndpse[1] = -fixef(mer.obj)[1]/fixef(mer.obj)[2] #pse_0
  jndpse[2] = -(fixef(mer.obj)[1]+fixef(mer.obj)[3])/(fixef(mer.obj)[2]+ fixef(mer.obj)[4]) #pse_32
  jndpse[3] = qnorm(0.75)/fixef(mer.obj)[2] #jnd_0
  jndpse[4] = qnorm(0.75)/(fixef(mer.obj)[2]+ fixef(mer.obj)[4]) #jnd_32
  jndpse[5] = jndpse[2] - jndpse[1] # Difference in PSE (32 Hz - 0 Hz)
  jndpse[6] = jndpse[4] - jndpse[3] # Difference in JND (32 Hz - 0 Hz)
  jndpse[7] = fixef(mer.obj)[4]
  return(jndpse)
}

boot_vibro <- pseMer(glmm_vibro, B = 500, FUN = fun2mod)

# BHMs data -----
attach(vibro_exp3)
vibr=0*speed
vibr[vibration=="32"]=1
soggetto=0*speed
soggetto[subject=="AK"]=1
soggetto[subject=="AR"]=2
soggetto[subject=="DN"]=3
soggetto[subject=="FA"]=4
soggetto[subject=="MA"]=5
soggetto[subject=="MI"]=6
soggetto[subject=="NI"]=7
soggetto[subject=="NN"]=8
soggetto[subject=="RV"]=9

datistan=list()
datistan$y=faster
datistan$n=faster+slower
datistan$nobs=length(datistan$y)
datistan$x=speed
datistan$vibration = vibr
datistan$subject = soggetto
datistan$nsubj=length(unique(datistan$subject))
datistan$ncond=length(unique(datistan$vibration))
detach(vibro_exp3)

# BH-GLM ----
init_fun <- function(init_alpha = -2, init_beta = 0.1) {
  list(
    alpha = matrix(init_alpha, datistan$nsubj, 2),
    beta  = matrix(init_beta, datistan$nsubj, 2)
  )
}
  
fit_bhglm_vibro <- stan(
  file = "Stan/vibro_bhglm.stan", 
  data = datistan,    
  chains = 3,             
  warmup = 3000,          
  iter = 5000,          
  cores = 3,              
  refresh = 10,             
  init=list(init_fun(),init_fun(),init_fun())
)  

## Residuals -----
samples_bhglm_vibro <- extract(fit_bhglm_vibro)
posteriorPredDistr <- samples_bhglm_vibro$predProb
posteriorPredSim   <- samples_bhglm_vibro$y_sim  

fitted_probs <- apply(posteriorPredDistr, 2, median)
y_obs <- datistan$y
trials <- datistan$n
SSE_bhglm_vibro <- sum( (y_obs/trials - fitted_probs)^2 )
loglik_vec <- dbinom(y_obs, size = trials, prob = fitted_probs, log = TRUE)
loglik_bhglm_vibro <- sum(loglik_vec)

## DHARMa diagnostic ------
sim_bhglm_vibro <- createDHARMa(
  simulatedResponse = t(posteriorPredSim),
  observedResponse = datistan$y,
  fittedPredictedResponse = apply(posteriorPredDistr, 2, median),
  integerResponse = TRUE
)

testUniformity(sim_bhglm_vibro)
testDispersion(sim_bhglm_vibro)
testOutliers(sim_bhglm_vibro)
testQuantiles(sim_bhglm_vibro)

# BH-GNM ------
init_fun <- function(init_alpha = -2, init_beta = 0.1, 
                     init_gamma = 0.01, init_lambda = 0.01) {
  list(
    alpha = matrix(init_alpha, datistan$nsubj, 2),
    beta  = matrix(init_beta, datistan$nsubj, 2),
    gamma  = matrix(init_gamma, datistan$nsubj, 2),
    lambda  = matrix(init_lambda, datistan$nsubj, 2)
  )
}

fit_bhgnm_vibro <- stan(
  file = "Stan/vibro_bhgnm.stan",  
  data = datistan,   
  chains = 3,            
  warmup = 3000,         
  iter = 5000,          
  cores = 3,             
  refresh = 10,
  init=list(init_fun(),init_fun(),init_fun()))  

## Residuals -----
samples_bhgnm_vibro <- extract(fit_bhgnm_vibro)
posteriorPredDistr <- samples_bhgnm_vibro$predProb
posteriorPredSim   <- samples_bhgnm_vibro$y_sim  

fitted_probs <- apply(posteriorPredDistr, 2, median)
y_obs <- datistan$y
trials <- datistan$n
SSE_bhgnm_vibro <- sum( (y_obs/trials - fitted_probs)^2 )
loglik_vec <- dbinom(y_obs, size = trials, prob = fitted_probs, log = TRUE)
loglik_bhgnm_vibro <- sum(loglik_vec)

## DHARMa diagnostic ------
sim_bhgnm_vibro <- createDHARMa(
  simulatedResponse = t(posteriorPredSim),
  observedResponse = datistan$y,
  fittedPredictedResponse = apply(posteriorPredDistr, 2, median),
  integerResponse = TRUE
)

testUniformity(sim_bhgnm_vibro)
testDispersion(sim_bhgnm_vibro)
testOutliers(sim_bhgnm_vibro)
testQuantiles(sim_bhgnm_vibro)
