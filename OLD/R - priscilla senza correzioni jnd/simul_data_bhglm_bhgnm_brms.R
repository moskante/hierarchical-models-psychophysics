# ==============================================================================
# utilizzo di brms (Sostituisce Stan puro)
# ==============================================================================

library(tidyverse)
library(MixedPsy)
library(brms)
library(DHARMa)

# 1. Preparazione Dati ----
set.seed(123)
simul_data <- PsySimulate(ntrials = 160, nsubjects = 10, guess = TRUE, lapse = TRUE)
simul_data$Subject <- factor(simul_data$Subject)

# 2. Modello 1: Bayes Hierarchical GLM (Senza Lapses) ----
# # Usiamo la sintassi brms che è identica a lme4
# 1. Elimina il file precedente se corrotto
if (file.exists("fits/bhglm_fit.rds")) file.remove("fits/bhglm_fit.rds")

# 2. Esegui il modello con prior e init ottimizzati
fit_bhglm <- brm(
  formula = Longer | trials(Total) ~ X + (1 + X | Subject),
  data = simul_data,
  family = binomial(link = "probit"), 
  prior = c(
    prior(normal(0, 3), class = "Intercept"), # Prior più ragionevoli per probit
    prior(normal(0, 3), class = "b"),
    prior(cauchy(0, 2.5), class = "sd"),
    prior(lkj(2), class = "cor")              # Prior per la correlazione random effects
  ),
  chains = 3, 
  iter = 5000, 
  warmup = 3000,
  refresh = 10, 
  init = 0,               # Fondamentale: forza l'inizio da zero per evitare log(0)
  backend = "rstan",      # O "cmdstanr" se lo hai configurato
  file = "fits/bhglm_fit"
)

# 3. Verifica l'oggetto prima di usarlo
if (is(fit_bhglm, "brmsfit")) {
  summary(fit_bhglm)
} else {
  print("Il modello non è stato creato correttamente.")
}


## ---  Predictive Distributions with brms ---
# Ottieni i conteggi attesi (n * pi)
pred_counts <- posterior_epred(fit_bhglm) 

# Dividi per il vettore dei trials per ottenere le probabilità (pi)
trials_vec <- simul_data$Total
posteriorPredDistr_brms <- sweep(pred_counts, 2, trials_vec, "/")

posteriorPredSim_brms <- posterior_predict(fit_bhglm)

# calcolo di fitted_probs e il SSE
fitted_probs <- apply(posteriorPredDistr_brms, 2, median)

# Estrai i dati osservati direttamente dal dataframe
y_obs <- simul_data$Longer
trials <- simul_data$Total

# Calcola il SSE (Sum of Squared Errors)
SSE_bhglm <- sum( (y_obs/trials - fitted_probs)^2 )

print(paste("SSE bhglm:", SSE_bhglm))


# Crea l'oggetto residui di DHARMa
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
  # 1. Estraiamo i coefficienti fissi (Fixed Effects)
  fe <- fixef(model)
  b0 <- fe["Intercept", "Estimate"]
  b1 <- fe["X", "Estimate"]
  
  # 2. Estraiamo gli effetti casuali (Random Effects per soggetto)
  re <- ranef(model)$Subject
  
  # 3. Calcolo dei parametri a livello di popolazione (Group Level)
  # PSE = -b0 / b1
  # Sigma = 1 / b1
  pse_group <- -b0 / b1
  sigma_group <- 1 / b1
  jnd_group <- 0.6745 * sigma_group
  
  # 4. Calcolo dei parametri per ogni singolo soggetto
  # Sommiamo l'effetto fisso all'effetto casuale specifico
  subjects <- rownames(re)
  subj_results <- data.frame(
    Subject = subjects,
    PSE = numeric(length(subjects)),
    Sigma = numeric(length(subjects)),
    JND = numeric(length(subjects))
  )
  
  for (i in seq_along(subjects)) {
    # Intercetta e pendenza specifica del soggetto
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


# 3. Modello 2: Bayes Hierarchical GNM (Con Lapses) ----
# Dovrebbe parzialmente corrispondere al vecchio simul_bhgnm.stan.
# Questo è il modello "Non-Lineare" che include gli asintoti.

# Definiamo la formula psicometrica: gamma + (1 - gamma - lambda) * Phi(beta0 + beta1*x)

formula_bhgnm <- bf(
  Longer | trials(Total) ~ gamma + (1 - gamma - lambda) * Phi((X - pse) / sigma),
  pse ~ 1 + (1 | Subject),
  sigma ~ 1 + (1 | Subject),
  
  gamma ~ 0 + Subject,
  lambda ~ 0 + Subject,
  
  nl = TRUE
)


# 2. Definiamo le Priors
# Nota: in brms dobbiamo specificare i limiti inferiori/superiori (lb/ub)
# per garantire che sigma sia positivo e che i rate siano tra 0 e 1.
priors_bhgnm <- c(
  # Intercette 
  prior(normal(0, 10), class = "b", nlpar = "pse"), 
  prior(normal(0, 10), class = "b", nlpar = "sigma", lb = 0),
  
  # Deviazioni standard di gruppo (equivalenti a tau_pse e tau_sigma)
  prior(cauchy(0, 2.5), class = "sd", nlpar = "pse"),
  prior(cauchy(0, 2.5), class = "sd", nlpar = "sigma"),
  
  # Priors per gamma e lambda (limitate tra 0 e 1)
  prior(uniform(0, 1), class = "b", nlpar = "gamma", lb = 0, ub = 1),
  prior(uniform(0, 1), class = "b", nlpar = "lambda", lb = 0, ub = 1)
)

# 3. Fit del Modello
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
# Ottieni i conteggi attesi (n * pi)
pred_counts <- posterior_epred(fit_bhgnm_brms) 

# Dividi per il vettore dei trials per ottenere le probabilità (pi)
trials_vec <- simul_data$Total
posteriorPredDistr_brms_bhgnm <- sweep(pred_counts, 2, trials_vec, "/")

posteriorPredSim_brms_bhgnm <- posterior_predict(fit_bhgnm_brms)

# calcolo di fitted_probs e il SSE
fitted_probs <- apply(posteriorPredDistr_brms_bhgnm, 2, median)

# Estrai i dati osservati direttamente dal dataframe
y_obs <- simul_data$Longer
trials <- simul_data$Total

# Calcola il SSE (Sum of Squared Errors)
SSE_bhgnm <- sum( (y_obs/trials - fitted_probs)^2 )

print(paste("SSE bhgnm:", SSE_bhgnm))



# Crea l'oggetto residui di DHARMa
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

extract_psychophysic_params_brms <- function(model) {
  # 1. Estraiamo gli effetti fissi (Population Level)
  fe <- fixef(model)
  
   pse_group   <- fe["pse_Intercept", "Estimate"]
  sigma_group <- fe["sigma_Intercept", "Estimate"]
  jnd_group   <- 0.6745 * sigma_group
  
  # 2. Estraiamo i coefficienti specifici per soggetto (Subject Level)
  subj_coefs <- coef(model)$Subject
  
  # 3. Creiamo il dataframe per i soggetti estraendo i valori
  subjects <- rownames(subj_coefs)
  
  subj_results <- data.frame(
    Subject = subjects,
    PSE     = subj_coefs[, "Estimate", "pse_Intercept"],
    Sigma   = subj_coefs[, "Estimate", "sigma_Intercept"]
  )
  
  # Calcoliamo la JND per i soggetti
  subj_results$JND <- 0.6745 * subj_results$Sigma
  
  # Pulizia degli indici di riga
  rownames(subj_results) <- NULL
  
  # 4. Restituiamo i risultati
  return(list(
    Population = data.frame(PSE = pse_group, Sigma = sigma_group, JND = jnd_group),
    Subjects   = subj_results
  ))
}


params <- extract_psychophysic_params_brms(fit_bhgnm_brms)
print(params$Population)
print(params$Subjects)

### calcolo dei valori Veri per verifica

true_params <- unique(simul_data[, c("Subject", "Intercept", "Slope")])
true_params$True_PSE   <- -true_params$Intercept / true_params$Slope
true_params$True_Sigma <- 1 / true_params$Slope
true_params$True_JND   <- 0.6745 * true_params$True_Sigma

print(true_params)



