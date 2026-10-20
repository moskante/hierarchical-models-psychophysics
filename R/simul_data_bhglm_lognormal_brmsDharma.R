library(brms)
library(posterior)
library(dplyr)
library(DHARMa)
library(tidyr)
library(stringr)
library(MixedPsy)

set.seed(123)

simul_data <- PsySimulate(
  ntrials   = 160,
  nsubjects = 10,
  guess     = TRUE,
  lapse     = TRUE
)

# ============================================================
#  PREPARE DATA
# ============================================================

# BHMs data ----
datistan        <- list()
datistan$y      <- simul_data$Longer
datistan$n      <- simul_data$Total
datistan$nobs   <- length(datistan$y)
datistan$x      <- simul_data$X
datistan$subject <- (simul_data$Subject)
datistan$nsubj  <- length(unique(datistan$subject))
datistan$ncond  <- length(unique(datistan$condition))

# Model: PI = Phi((x - pse) / sigma)
# pse_i   ~ Normal(PSE, tau_pse)
# lsd_i   = log(sigma_i) ~ Normal(log_SIGMA, tau_sigma)  
formula <- bf(
  y | trials(n) ~ Phi((x - pse) / sigmai),
  nlf(sigmai ~ exp(lsd)),
  pse ~ 1 + (1 | subject),
  lsd ~ 1 + (1 | subject),
  nl = TRUE
)

# Famiglia binomiale con link identity (la probit è già dentro Phi)
fam <- binomial(link = "identity")

init_fun_brms <- function() {
  list(
    b_pse     = as.array(90),               # PSE 
    b_lsd     = as.array(log(15))     # sigma iniziale 
   # sd_1      = as.array(c(15 / 10)),      # sd random effect pse (piccola ma non zero)
  #  sd_2      = as.array(0.1),              # sd random effect lsd
  #  tauPSE   = 15,
  #  tauSIGMA = 1
  )
}
gq_code <- "
  real PSE = b_pse[1];
  real log_SIGMA = b_lsd[1];
  real SIGMA = exp(b_lsd[1] + 0.5 * square(sd_2[1]));  // E[sigma_i]
  real JND = 0.6745 * SIGMA;
  vector[N_1] pse;
  vector[N_2] sigma;
  vector[N_2] jnd;
  vector[N] predProb;
  array[N] int y_sim;

  for (j in 1:N_1) pse[j] = b_pse[1] + r_1_pse_1[j];
  for (j in 1:N_2) {
    sigma[j] = exp(b_lsd[1] + r_2_lsd_1[j]);
    jnd[j]   = 0.6745 * sigma[j];
  }
  for (n in 1:N) {
    predProb[n] = Phi((C_1[n] - pse[J_1[n]]) / sigma[J_2[n]]);
    y_sim[n]    = binomial_rng(trials[n], predProb[n]);
  }
"



stanvars <- 
  stanvar(scode = "  real<lower=0> tauPSE;",   block = "parameters") +
  stanvar(scode = "  real<lower=0> tauSIGMA;", block = "parameters") +
  stanvar(scode = "  tauPSE ~ cauchy(0, 2.5);
  tauSIGMA ~ cauchy(0, 2.5);", block = "model")+
  stanvar(scode = gq_code, block = "genquant")

priors <- c(
  prior(normal(0, tauPSE),   class = "b", coef = "Intercept", nlpar = "pse"),
  prior(normal(0, tauSIGMA), class = "b", coef = "Intercept", nlpar = "lsd"),
  prior(cauchy(0, 2.5), class = "sd", nlpar = "pse"),   # tau_pse
  prior(cauchy(0, 2.5), class = "sd", nlpar = "lsd")    # tau_sigma
)

fit2 <- brm(
  formula, data = datistan, family = fam,
  prior = priors, stanvars = stanvars,
  chains = 3, iter = 5000, warmup = 3000, cores = 3,
  init = init_fun_brms,
  backend = "rstan"
)



draws <- as_draws_df(fit2)
# equivalent to hyperparameters_bhglm
hyper_brms <- summarise_draws(
  subset_draws(draws, variable = c("PSE", "SIGMA", "JND", "tauPSE", "tauSIGMA")),
  mean, sd,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  rhat, ess_bulk, ess_tail
)

###STAN implementation (comparison scope)
# BH-GLM ----
init_fun <- function(init_pse = 90, init_sigma = 15) {
  list(
    pse   = rep(init_pse,   datistan$nsubj),
    sigma = rep(init_sigma, datistan$nsubj)
  )
}

fit_bhglm <- stan(
  file   = "Stan/simul_bhglm_lognormal.stan",
  data   = datistan,
  chains = 3,
  warmup = 3000,
  iter   = 5000,   # total iterations per chain (including warmup)
  cores  = 3,
  refresh = 10,
  init   = list(init_fun(), init_fun(), init_fun())
)



hyperparameters_bhglm <- summary(fit_bhglm)$summary %>%
  as_tibble(rownames = "params") %>%
  filter(str_detect(params, "PSE|JND"), params != "tau_PSE")

hyper_brms %>% filter(variable %in% c("PSE", "JND")) %>% select(variable, mean, q2.5, q97.5)
hyperparameters_bhglm%>% filter(params %in% c("PSE", "JND")) %>% select(params, mean, `2.5%`, `97.5%`)



draws <- as_draws_matrix(fit2)
variables(draws)   # controllo

# Iperparametri (equivalente a hyperparameters_bhglm)
hyper_brms <- summarise_draws(
  subset_draws(draws, variable = c("PSE", "SIGMA", "JND")),
  mean, sd, ~quantile2(.x, probs = c(0.025, 0.975)), rhat, ess_bulk
)

# Parametri per soggetto (equivalente a parameters_bhglm)
parameters_brms <- summarise_draws(
  subset_draws(draws, variable = c("pse", "jnd")),
  mean, ~quantile2(.x, probs = c(0.025, 0.975)), rhat, ess_bulk
) %>%
  separate(variable, into = c("params", "Subject"), sep = "\\[|\\]", extra = "drop")

# SSE e DHARMa
y_sim    <- as.matrix(subset_draws(draws, variable = "y_sim"))
predProb <- as.matrix(subset_draws(draws, variable = "predProb"))
fitted_probs   <- apply(predProb, 2, median)
SSE_bhglm_brms <- sum((datistan$y / datistan$n - fitted_probs)^2)


y_mat <- matrix(as.numeric(y_sim), nrow = nrow(y_sim), ncol = ncol(y_sim))  

    
sim_bhglm_brms <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit2)),
  observedResponse = datistan$y,
  fittedPredictedResponse = apply(t(posterior_epred(fit2)), 1, mean),
  integerResponse = TRUE)


plot(sim_bhglm_brms); testUniformity(sim_bhglm_brms); testDispersion(sim_bhglm_brms)
testOutliers(sim_bhglm_brms); testQuantiles(sim_bhglm_brms)
# Per il grafico di confronto
df_bhglm <- hyperparameters_bhglm %>%
  rename(Parameter = params, Estimate = mean,
         Lower = `2.5%`, Upper = `97.5%`) %>%
  select(Parameter, Estimate, Lower, Upper) %>%
  mutate(Model = "BH-GLM (Stan)")

# brms
df_bhglm_brms <- hyper_brms %>%
  filter(variable %in% c("PSE", "JND")) %>%
  transmute(Parameter = variable, Estimate = mean,
            Lower = q2.5, Upper = q97.5, Model = "BH-GLM (brms)")

# Combinazione
df_combined <- bind_rows(df_bhglm, df_bhglm_brms) %>%
  mutate(Model = factor(Model,
                        levels = c("BH-GLM (Stan)", "BH-GLM (brms)")))

