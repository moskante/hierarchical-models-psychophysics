# ============================================================================
# Monte Carlo simulation study: cross-method comparison of GLM, GNM, GLMM,
# BH-GLM, and BH-GNM on independently simulated psychophysical datasets.
# ============================================================================
# NOMENCLATURE & PARAMETER DEFINITIONS
# ============================================================================
# The generative population fixed effects represent the sole target for 
# calculating recovery biases and RMSE.
# ============================================================================

library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(rstan)
library(gnlm)

source("R/gnlm_functions_psejnd.R")   # contains PsychParametersGNM()

# ============================================================
# Monte Carlo loop configuration
# ============================================================
n_target  <- 150    # number of successful iterations required per method
ntrials   <- 160
nsubjects <- 10
run_stan  <- TRUE   # set FALSE for quick GLM/GNM/GLMM-only runs

# --- True Population Values (Generative Ground Truth) ---
true_fixeff_mu      <- -(-7) / 0.0875        

# Theoretical means of generative CRAN uniform ranges c(0, 0.05)
true_gamma          <- 0.025                 
true_lambda         <- 0.025                 
true_sigma          <- 1 / 0.0875
true_k_val          <- qnorm((0.5 - true_gamma) / (1 - true_gamma - true_lambda)) * true_sigma

true_pse            <- true_fixeff_mu + true_k_val
true_jnd            <- qnorm((0.75 - true_gamma) / (1 - true_gamma - true_lambda)) * true_sigma - true_k_val

if (run_stan) {
  stan_bhglm <- stan_model(file = "Stan/simul_bhglm.stan")
  stan_bhgnm <- stan_model(file = "Stan/simul_bhgnm_GL.stan") 
}

# ============================================================
# Clean Bayesian Helpers
# ============================================================
stan_sse <- function(stan_fit, y_obs, trials) {
  fitted_probs <- apply(extract(stan_fit)$predProb, 2, median)
  sum((y_obs / trials - fitted_probs)^2)
}

# Super-simplified parameters extraction (no scale factors needed!)
extract_stan_all_params <- function(stan_fit) {
  s <- extract(stan_fit)
  
  list(
    PSE    = median(s$PSE, na.rm = TRUE), 
    JND    = median(s$JND, na.rm = TRUE), 
    GAMMA  = if ("GAMMA" %in% names(s)) median(s$GAMMA, na.rm = TRUE) else 0,
    LAMBDA = if ("LAMBDA" %in% names(s)) median(s$LAMBDA, na.rm = TRUE) else 0
  )
}

simulate_dataset <- function() {
  simul_data         <- PsySimulate(ntrials   = ntrials,
                                    nsubjects = nsubjects,
                                    guess     = TRUE,
                                    lapse     = TRUE)
  simul_data$Subject <- factor(simul_data$Subject)
  
  # Estraiamo la stima campionaria reale dei soggetti generati in questa iterazione[cite: 5]
  params_true <- simul_data %>%
    group_by(Subject) %>%
    summarise(
      Intercept_i = Intercept[1],
      Slope_i     = Slope[1],
      .groups     = "drop"
    ) %>%
    mutate(
      pse_sim = -Intercept_i / Slope_i
    )
  
  # Ritorniamo sia il dataset che la stima campionaria per l'inizializzazione
  list(
    data            = simul_data,
    true_sample_pse = mean(params_true$pse_sim, na.rm = TRUE)
  )
}

# ============================================================
# Fit functions
# ============================================================
fit_GLM <- function(iter, sim) {
  tryCatch({
    simul_data <- sim$data
    glm_list   <- PsychModels(cbind(Longer, Total - Longer) ~ X, data = simul_data, group_factors = "Subject")
    params_glm <- PsychParameters(glm_list) %>% 
      mutate(Subject = factor(Subject, levels = levels(simul_data$Subject)))
    
    sse_glm <- simul_data %>%
      left_join(params_glm, by = "Subject") %>%
      mutate(eta = qnorm(0.75) * (X - pse) / jnd, pred_prob = pnorm(eta), obs_prop = Longer / Total) %>%
      summarise(sse = sum((obs_prop - pred_prob)^2, na.rm = TRUE)) %>% pull(sse)
    
    # Test against the absolute fixed-effect latent targets
    t_pse <- t.test(params_glm$pse, mu = true_pse)
    t_jnd <- t.test(params_glm$jnd, mu = true_jnd)
    
    tibble(iter       = iter, 
           method     = "GLM",
           fit_pse    = mean(params_glm$pse, na.rm = TRUE), 
           fit_jnd    = mean(params_glm$jnd, na.rm = TRUE),
           fit_gamma  = 0, 
           fit_lambda = 0, 
           t_pse_stat = as.numeric(t_pse$statistic), t_pse_p = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic), t_jnd_p = t_jnd$p.value, 
           SSE        = sse_glm)
  }, error = function(e) NULL)
}

fit_GNM <- function(iter, sim) {
  tryCatch({
    simul_data    <- sim$data
    subjects_data <- simul_data %>% group_split(Subject)
    sub_names     <- map_chr(subjects_data, ~ as.character(.x$Subject[1]))
    
    params_gnm <- map_dfr(subjects_data, process_subject) %>% 
      mutate(Subject = factor(sub_names, levels = levels(simul_data$Subject)))
    
    sse_gnm <- simul_data %>%
      left_join(params_gnm, by = "Subject") %>%
      mutate(
        v_pse = (0.5 - gamma) / (1 - gamma - lambda),
        v_jnd = (0.75 - gamma) / (1 - gamma - lambda),
        denom = qnorm(pmax(1e-5, pmin(1-1e-5, v_jnd))) - qnorm(pmax(1e-5, pmin(1-1e-5, v_pse))),
        subj_sigma = jnd / denom,
        subj_mu    = pse - qnorm(pmax(1e-5, pmin(1-1e-5, v_pse))) * subj_sigma,
        eta = (X - subj_mu) / subj_sigma,
        pred_prob = gamma + (1 - gamma - lambda) * pnorm(eta),
        obs_prop = Longer / Total
      ) %>%
      summarise(sse = sum((obs_prop - pred_prob)^2, na.rm = TRUE)) %>% pull(sse)
    
    # GNM recovers adjusted metrics, so test against the adjusted targets
    t_pse <- t.test(params_gnm$pse, mu = true_pse)
    t_jnd <- t.test(params_gnm$jnd, mu = true_jnd)
    
    tibble(iter       = iter, 
           method     = "GNM",
           fit_pse    = mean(params_gnm$pse, na.rm = TRUE), 
           fit_jnd    = mean(params_gnm$jnd, na.rm = TRUE),
           fit_gamma  = mean(params_gnm$gamma, na.rm = TRUE), 
           fit_lambda = mean(params_gnm$lambda, na.rm = TRUE),
           t_pse_stat = as.numeric(t_pse$statistic), t_pse_p = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic), t_jnd_p = t_jnd$p.value, 
           SSE        = sse_gnm)
  }, error = function(e) NULL)
}

fit_GLMM <- function(iter, sim) {
  tryCatch({
    simul_data <- sim$data
    glmm_fit   <- glmer(cbind(Longer, Total - Longer) ~ X + (1 + X | Subject), family = binomial(link = "probit"), data = simul_data)
    
    fe       <- fixef(glmm_fit)
    fit_pse  <- as.numeric(-fe["(Intercept)"] / fe["X"])
    fit_jnd  <- as.numeric(qnorm(0.75) / fe["X"])
    
    tibble(iter       = iter, 
           method     = "GLMM",
           fit_pse    = fit_pse, 
           fit_jnd    = fit_jnd,
           fit_gamma  = 0, 
           fit_lambda = 0, 
           t_pse_stat = NA_real_, t_pse_p = NA_real_, t_jnd_stat = NA_real_, t_jnd_p = NA_real_, 
           SSE        = sum(residuals(glmm_fit, type = "response")^2))
  }, error = function(e) NULL)
}

fit_BHGLM <- function(iter, sim) {
  tryCatch({
    simul_data      <- sim$data
    true_sample_pse <- sim$true_sample_pse 
    
    datistan <- list(y = simul_data$Longer, n = simul_data$Total, 
                     nobs = nrow(simul_data), x = simul_data$X, 
                     subject = as.integer(simul_data$Subject), nsubj = nlevels(simul_data$Subject))
    
    # UNIFIED MINIMAL INITIALIZATION: dynamically seeded by sample-level PSE[cite: 5]
    init_fn  <- function() list(
      pse   = rep(true_sample_pse, datistan$nsubj),
      sigma = rep(15,              datistan$nsubj)
    )
    
    fit <- sampling(stan_bhglm, data = datistan, chains = 3, warmup = 3000,
                    iter = 5000, cores = 3, refresh = 0, init = list(init_fn(), init_fn(), init_fn()))
    p   <- extract_stan_all_params(fit)
    
    tibble(iter       = iter, 
           method     = "BH-GLM",
           fit_pse    = p$PSE, 
           fit_jnd    = p$JND,
           fit_gamma  = 0, 
           fit_lambda = 0, 
           t_pse_stat = NA_real_, t_pse_p = NA_real_, t_jnd_stat = NA_real_, t_jnd_p = NA_real_, 
           SSE        = stan_sse(fit, datistan$y, datistan$n))
  }, error = function(e) NULL)
}

fit_BHGNM <- function(iter, sim) {
  tryCatch({
    simul_data      <- sim$data
    true_sample_pse <- sim$true_sample_pse # dynamically unpacked[cite: 5]
    
    datistan <- list(y = simul_data$Longer, n = simul_data$Total, nobs = nrow(simul_data),
                     x = simul_data$X, subject = as.integer(simul_data$Subject), nsubj = nlevels(simul_data$Subject))
    
    # UNIFIED MINIMAL INITIALIZATION: identically seeded by sample-level PSE[cite: 5]
    init_fn  <- function() list(
      mu    = rep(true_sample_pse, datistan$nsubj),
      sigma = rep(15,              datistan$nsubj)
    )
    
    fit <- sampling(stan_bhgnm, data = datistan, chains = 3, 
                    warmup = 3000, iter = 5000, cores = 3, refresh = 0, init = list(init_fn(), init_fn(), init_fn()))
    p   <- extract_stan_all_params(fit)
    
    tibble(iter       = iter, 
           method     = "BH-GNM",
           fit_pse    = p$PSE, 
           fit_jnd    = p$JND,
           fit_gamma  = p$GAMMA, 
           fit_lambda = p$LAMBDA,
           t_pse_stat = NA_real_, t_pse_p = NA_real_, t_jnd_stat = NA_real_, t_jnd_p = NA_real_, 
           SSE        = stan_sse(fit, datistan$y, datistan$n))
  }, error = function(e) NULL)
}

# ============================================================
# Main Loop Execution
# ============================================================
all_methods  <- if (run_stan) c("GLM", "GNM", "GLMM", "BH-GLM", "BH-GNM") else c("GLM", "GNM", "GLMM")
results_list <- setNames(vector("list", length(all_methods)), all_methods)
n_done       <- setNames(integer(length(all_methods)), all_methods)
iter_counter <- 0L

while (any(n_done < n_target)) {
  sim          <- simulate_dataset() # Sim output contains dataset + true_sample_pse[cite: 5]
  iter_counter <- iter_counter + 1L
  for (m in all_methods) {
    if (n_done[m] >= n_target) next
    row <- switch(m,
                  "GLM"    = fit_GLM(  iter_counter, sim),
                  "GNM"    = fit_GNM(  iter_counter, sim),
                  "GLMM"   = fit_GLMM( iter_counter, sim),
                  "BH-GLM" = fit_BHGLM(iter_counter, sim),
                  "BH-GNM" = fit_BHGNM(iter_counter, sim))
    if (!is.null(row)) {
      n_done[m] <- n_done[m] + 1L
      results_list[[m]][[n_done[m]]] <- row
    }
  }
  message(sprintf("Attempt %d | successes: %s", iter_counter, paste(names(n_done), n_done, sep = "=", collapse = " | ")))
}

# ============================================================
# Post-Processing & Tables
# ============================================================
results <- bind_rows(lapply(results_list, bind_rows)) %>%
  mutate(
    # 1. Behavioral Biases (vs Theoretical Adjusted Truth)
    bias_pse       = fit_pse - true_pse,
    bias_jnd       = fit_jnd - true_jnd,
    
    # 2. Asymptotes Biases
    bias_gamma     = fit_gamma - true_gamma,
    bias_lambda    = fit_lambda - true_lambda
  )

# ============================================================
# Summary Tables
# ============================================================
build_summary <- function(df) {
  df %>%
    group_by(method) %>%
    summarise(
      # --- Behavioral Performance (Adjusted Truth) ---
      bias_pse       = mean(bias_pse, na.rm = TRUE),
      rmse_pse       = sqrt(mean( (fit_pse - true_pse)^2, na.rm = TRUE)),
      bias_jnd       = mean(bias_jnd, na.rm = TRUE),
      rmse_jnd       = sqrt(mean( (fit_jnd - true_jnd)^2, na.rm = TRUE)),
      
      # --- Asymptotes ---
      bias_gamma     = mean(bias_gamma, na.rm = TRUE),
      rmse_gamma     = sqrt(mean( (fit_gamma - true_gamma)^2, na.rm = TRUE)),
      bias_lambda    = mean(bias_lambda, na.rm = TRUE),
      rmse_lambda    = sqrt(mean( (fit_lambda - true_lambda)^2, na.rm = TRUE)),
      
      mean_SSE       = mean(SSE, na.rm = TRUE),
      .groups = "drop"
    )
}

print("--- GLOBAL SUMMARY TABLE ---")
print(build_summary(results))

print("--- POSITIVE PSE SUMMARY TABLE ---")
print(build_summary(results %>% dplyr::filter(fit_pse > 0)))