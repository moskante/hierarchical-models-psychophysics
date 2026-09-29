# ============================================================================
# Monte Carlo simulation study: cross-method comparison of GLM, GNM, GLMM,
# BH-GLM, and BH-GNM on independently simulated psychophysical datasets.
#
# ============================================================================
# NOMENCLATURE & PARAMETER DEFINITIONS
# ============================================================================
# 1. "true_" Prefixes (Generative Ground Truth):
#    - true_fixeff_pse / true_fixeff_jnd: The fixed, absolute parameters of the 
#      underlying population distribution used to seed the simulation.
#    - true_sample_pse / true_sample_jnd: The empirical arithmetic mean of the
#      10 specific subjects drawn in a given iteration (includes sample noise).
#
# 2. "fit_" Prefixes (Estimated Model Parameters):
#    - fit_sample_mean_pse / fit_sample_mean_jnd: The mean of the subject-specific
#      parameters estimated by the model for that iteration's sample (simplified).
#    - fit_sample_mean_pse_adj / fit_sample_mean_jnd_adj: The mean of the subject-specific
#      parameters mathematically adjusted to account for estimated guess/lapse rates.
#    - fit_pop_pse / fit_pop_jnd: The global population hyper-parameters directly
#      estimated by hierarchical models (fixed effects in GLMM, hyper-medians in Stan),
#      or the sample mean proxy for non-hierarchical models (GLM, GNM).
# ============================================================================

library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(rstan)
library(gnlm)

source("../R/gnlm_functions_psejnd.R")   # defines process_subject() for GNM fitting

# ============================================================
# Monte Carlo loop configuration
# ============================================================
n_target  <- 150    # number of successful iterations required per method
ntrials   <- 160
nsubjects <- 10
run_stan  <- TRUE   # set FALSE for quick GLM/GNM/GLMM-only runs

# --- True Population Values (Generative Ground Truth) ---
true_fixeff_pse     <- -(-7) / 0.0875        
true_fixeff_jnd     <- qnorm(0.75) / 0.0875  

# Mathematically Adjusted Population Parameters (using mean of uniform ranges)
pop_mean_gamma      <- 0.075                 # mean of c(0.05, 0.10)
pop_mean_lambda     <- 0.025                 # mean of c(0, 0.05)
pop_sigma           <- 1 / 0.0875
k_pop_val           <- qnorm((0.5 - pop_mean_gamma) / (1 - pop_mean_gamma - pop_mean_lambda)) * pop_sigma

true_fixeff_pse_adj <- true_fixeff_pse + k_pop_val
true_fixeff_jnd_adj <- qnorm((0.75 - pop_mean_gamma) / (1 - pop_mean_gamma - pop_mean_lambda)) * pop_sigma - k_pop_val

if (run_stan) {
  stan_bhglm <- stan_model(file = "../Stan/simul_bhglm.stan")
  stan_bhgnm <- stan_model(file = "../Stan/simul_bhgnm.stan")
}

# ============================================================
# Helpers
# ============================================================
stan_sse <- function(stan_fit, y_obs, trials) {
  fitted_probs <- apply(extract(stan_fit)$predProb, 2, median)
  sum((y_obs / trials - fitted_probs)^2)
}

# Sample-level helper (Mean of subject posterior medians)
stan_params <- function(stan_fit, param_pse = "pse", param_jnd = "jnd") {
  s <- extract(stan_fit)
  list(pse = mean(apply(s[[param_pse]], 2, median), na.rm = TRUE),
       jnd = mean(apply(s[[param_jnd]], 2, median), na.rm = TRUE))
}

# Population-level helper (Direct hyper-parameter estimation)
stan_pop_params <- function(stan_fit, param_pse = "PSE", param_jnd = "JND") {
  s <- extract(stan_fit)
  list(pse = median(s[[param_pse]], na.rm = TRUE),
       jnd = median(s[[param_jnd]], na.rm = TRUE))
}

# Bayesian Adjusted parameters helper (incorporating gamma and lambda)
stan_params_adjusted <- function(stan_fit) {
  s <- extract(stan_fit)
  nsubj <- dim(s$pse)[2]
  
  adj_pse_subj <- numeric(nsubj)
  adj_jnd_subj <- numeric(nsubj)
  
  for (j in 1:nsubj) {
    pse_draws    <- s$pse[, j]
    sigma_draws  <- s$sigma[, j]
    gamma_draws  <- s$gamma[, j]
    lambda_draws <- s$lambda[, j]
    
    # Calculate Shift k_i
    val_pse <- (0.5 - gamma_draws) / (1 - gamma_draws - lambda_draws)
    val_pse <- pmax(1e-5, pmin(1 - 1e-5, val_pse)) 
    k_draws <- qnorm(val_pse) * sigma_draws
    
    # Calculate Adjusted JND(0.75)
    val_jnd <- (0.75 - gamma_draws) / (1 - gamma_draws - lambda_draws)
    val_jnd <- pmax(1e-5, pmin(1 - 1e-5, val_jnd))
    jnd_draws <- qnorm(val_jnd) * sigma_draws - k_draws
    
    adj_pse_subj[j] <- median(pse_draws + k_draws, na.rm = TRUE)
    adj_jnd_subj[j] <- median(jnd_draws, na.rm = TRUE)
  }
  
  list(pse = mean(adj_pse_subj, na.rm = TRUE),
       jnd = mean(adj_jnd_subj, na.rm = TRUE))
}

simulate_dataset <- function() {
  simul_data         <- PsySimulate(ntrials   = ntrials,
                                    nsubjects = nsubjects,
                                    guess     = TRUE,
                                    lapse     = TRUE)
  simul_data$Subject <- factor(simul_data$Subject)
  params_true <- simul_data %>%
    group_by(Subject) %>%
    summarise(across(everything(), first),
              pse = -Intercept / Slope,
              jnd = qnorm(0.75) / Slope,
              .groups = "drop")
  list(data            = simul_data,
       true_sample_pse = mean(params_true$pse),
       true_sample_jnd = mean(params_true$jnd))
}

# ============================================================
# Fit functions
# ============================================================
fit_GLM <- function(iter, sim) {
  tryCatch({
    simul_data      <- sim$data
    true_sample_pse <- sim$true_sample_pse
    true_sample_jnd <- sim$true_sample_jnd
    
    glm_list   <- PsychModels(cbind(Longer, Total - Longer) ~ X,
                              data = simul_data,
                              group_factors = "Subject")
    params_glm <- PsychParameters(glm_list)
    
    params_glm_clean <- params_glm %>%
      mutate(Subject = factor(Subject, levels = levels(simul_data$Subject)))
    
    sse_glm <- simul_data %>%
      left_join(params_glm_clean, by = "Subject") %>%
      mutate(
        eta       = qnorm(0.75) * (X - pse) / jnd,
        pred_prob = pnorm(eta),
        obs_prop  = Longer / Total
      ) %>%
      summarise(sse = sum((obs_prop - pred_prob)^2, na.rm = TRUE)) %>%
      pull(sse)
    
    t_pse      <- t.test(params_glm$pse,  mu = true_sample_pse)
    t_jnd      <- t.test(params_glm$jnd,  mu = true_sample_jnd)
    
    m_pse <- mean(params_glm$pse, na.rm = TRUE)
    m_jnd <- mean(params_glm$jnd, na.rm = TRUE)
    
    tibble(iter = iter, method = "GLM",
           fit_sample_mean_pse     = m_pse, fit_sample_mean_jnd     = m_jnd,
           fit_sample_mean_pse_adj = m_pse, fit_sample_mean_jnd_adj = m_jnd, # Same (gamma=lambda=0)
           fit_pop_pse             = m_pse, fit_pop_jnd             = m_jnd, 
           t_pse_stat = as.numeric(t_pse$statistic), t_pse_p = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic), t_jnd_p = t_jnd$p.value,
           SSE = sse_glm, true_sample_pse = true_sample_pse, true_sample_jnd = true_sample_jnd,
           true_fixeff_pse = true_fixeff_pse_val, true_fixeff_jnd = true_fixeff_jnd_val)
  }, error = function(e) NULL)
}

fit_GNM <- function(iter, sim) {
  tryCatch({
    simul_data      <- sim$data
    true_sample_pse <- sim$true_sample_pse
    true_sample_jnd <- sim$true_sample_jnd
    
    subjects_data <- simul_data %>% group_split(Subject)
    sub_names     <- map_chr(subjects_data, ~ as.character(.x$Subject[1]))
    params_gnm    <- map_dfr(subjects_data, process_subject)
    
    g_col <- intersect(c("gamma", "guess"), names(params_gnm))
    l_col <- intersect(c("lambda", "lapse"), names(params_gnm))
    
    params_gnm_clean <- params_gnm %>%
      mutate(
        Subject  = factor(sub_names, levels = levels(simul_data$Subject)),
        p_gamma  = if (length(g_col) > 0) .data[[g_col[1]]] else 0,
        p_lambda = if (length(l_col) > 0) .data[[l_col[1]]] else 0
      )
    
    sse_gnm <- simul_data %>%
      left_join(params_gnm_clean, by = "Subject") %>%
      mutate(
        eta       = qnorm(0.75) * (X - pse) / jnd,
        pred_prob = p_gamma + (1 - p_gamma - p_lambda) * pnorm(eta),
        obs_prop  = Longer / Total
      ) %>%
      summarise(sse = sum((obs_prop - pred_prob)^2, na.rm = TRUE)) %>%
      pull(sse)
    
    t_pse <- t.test(params_gnm$pse, mu = true_sample_pse)
    t_jnd <- t.test(params_gnm$jnd, mu = true_sample_jnd)
    
    m_pse <- mean(params_gnm$pse, na.rm = TRUE)
    m_jnd <- mean(params_gnm$jnd, na.rm = TRUE)
    
    # --- GNM Point-estimate Adjustment ---
    sigmas  <- params_gnm_clean$jnd / qnorm(0.75)
    gammas  <- params_gnm_clean$p_gamma
    lambdas <- params_gnm_clean$p_lambda
    
    val_pse <- (0.5 - gammas) / (1 - gammas - lambdas)
    val_pse <- pmax(1e-5, pmin(1 - 1e-5, val_pse))
    k_vals  <- qnorm(val_pse) * sigmas
    
    val_jnd <- (0.75 - gammas) / (1 - gammas - lambdas)
    val_jnd <- pmax(1e-5, pmin(1 - 1e-5, val_jnd))
    
    jnd_adj_vals <- qnorm(val_jnd) * sigmas - k_vals
    pse_adj_vals <- params_gnm_clean$pse + k_vals
    
    m_pse_adj <- mean(pse_adj_vals, na.rm = TRUE)
    m_jnd_adj <- mean(jnd_adj_vals, na.rm = TRUE)
    
    tibble(iter = iter, method = "GNM",
           fit_sample_mean_pse     = m_pse,     fit_sample_mean_jnd     = m_jnd,
           fit_sample_mean_pse_adj = m_pse_adj, fit_sample_mean_jnd_adj = m_jnd_adj,
           fit_pop_pse             = m_pse,     fit_pop_jnd             = m_jnd,
           t_pse_stat = as.numeric(t_pse$statistic), t_pse_p = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic), t_jnd_p = t_jnd$p.value,
           SSE = sse_gnm, true_sample_pse = true_sample_pse, true_sample_jnd = true_sample_jnd,
           true_fixeff_pse = true_fixeff_pse_val, true_fixeff_jnd = true_fixeff_jnd_val)
  }, error = function(e) NULL)
}

fit_GLMM <- function(iter, sim) {
  tryCatch({
    simul_data      <- sim$data
    true_sample_pse <- sim$true_sample_pse
    true_sample_jnd <- sim$true_sample_jnd
    
    glmm_fit <- glmer(
      cbind(Longer, Total - Longer) ~ X + (1 + X | Subject),
      family = binomial(link = "probit"), data = simul_data
    )
    
    cc       <- coef(glmm_fit)$Subject
    subj_pse <- -cc[,"(Intercept)"] / cc[,"X"]
    subj_jnd <- qnorm(0.75) / cc[,"X"]
    
    fe       <- fixef(glmm_fit)
    m_pse    <- mean(subj_pse, na.rm = TRUE)
    m_jnd    <- mean(subj_jnd, na.rm = TRUE)
    
    tibble(iter = iter, method = "GLMM",
           fit_sample_mean_pse     = m_pse, fit_sample_mean_jnd     = m_jnd,
           fit_sample_mean_pse_adj = m_pse, fit_sample_mean_jnd_adj = m_jnd, # Same (gamma=lambda=0)
           fit_pop_pse         = as.numeric(-fe["(Intercept)"] / fe["X"]),
           fit_pop_jnd         = as.numeric(qnorm(0.75) / fe["X"]),
           t_pse_stat = NA_real_, t_pse_p = NA_real_,
           t_jnd_stat = NA_real_, t_jnd_p = NA_real_,
           SSE        = sum(residuals(glmm_fit, type = "response")^2),
           true_sample_pse = true_sample_pse, true_sample_jnd = true_sample_jnd,
           true_fixeff_pse = true_fixeff_pse_val, true_fixeff_jnd = true_fixeff_jnd_val)
  }, error = function(e) NULL)
}

fit_BHGLM <- function(iter, sim) {
  tryCatch({
    simul_data      <- sim$data
    true_sample_pse <- sim$true_sample_pse
    true_sample_jnd <- sim$true_sample_jnd
    
    datistan <- list(
      y = simul_data$Longer, n = simul_data$Total, nobs = nrow(simul_data),
      x = simul_data$X, subject = as.integer(simul_data$Subject), nsubj = nlevels(simul_data$Subject)
    )
    
    init_fn <- function() list(pse   = rep(true_sample_pse, datistan$nsubj),
                               sigma = rep(15,       datistan$nsubj))
    fit <- sampling(stan_bhglm, data = datistan,
                    chains = 3, warmup = 3000, iter = 5000,
                    cores = 3, refresh = 0,
                    init = list(init_fn(), init_fn(), init_fn()))
    
    p_samp <- stan_params(fit, param_pse = "pse", param_jnd = "jnd")
    p_pop  <- stan_pop_params(fit, param_pse = "PSE", param_jnd = "JND")
    
    tibble(iter = iter, method = "BH-GLM",
           fit_sample_mean_pse     = p_samp$pse, fit_sample_mean_jnd     = p_samp$jnd,
           fit_sample_mean_pse_adj = p_samp$pse, fit_sample_mean_jnd_adj = p_samp$jnd, # Same (gamma=lambda=0)
           fit_pop_pse             = p_pop$pse,  fit_pop_jnd             = p_pop$jnd,
           t_pse_stat = NA_real_, t_pse_p = NA_real_,
           t_jnd_stat = NA_real_, t_jnd_p = NA_real_,
           SSE = stan_sse(fit, datistan$y, datistan$n),
           true_sample_pse = true_sample_pse, true_sample_jnd = true_sample_jnd,
           true_fixeff_pse = true_fixeff_pse_val, true_fixeff_jnd = true_fixeff_jnd_val)
  }, error = function(e) NULL)
}

fit_BHGNM <- function(iter, sim) {
  tryCatch({
    simul_data      <- sim$data
    true_sample_pse <- sim$true_sample_pse
    true_sample_jnd <- sim$true_sample_jnd
    
    datistan <- list(
      y = simul_data$Longer, n = simul_data$Total, nobs = nrow(simul_data),
      x = simul_data$X, subject = as.integer(simul_data$Subject), nsubj = nlevels(simul_data$Subject)
    )
    
    init_fn <- function() list(pse    = rep(true_sample_pse, datistan$nsubj),
                               sigma  = rep(15,       datistan$nsubj),
                               gamma  = rep(0.01,     datistan$nsubj),
                               lambda = rep(0.01,     datistan$nsubj))
    fit <- sampling(stan_bhgnm, data = datistan,
                    chains = 3, warmup = 3000, iter = 5000,
                    cores = 3, refresh = 0,
                    init = list(init_fn(), init_fn(), init_fn()))
    
    p_samp     <- stan_params(fit, param_pse = "pse", param_jnd = "jnd")
    p_pop      <- stan_pop_params(fit, param_pse = "PSE", param_jnd = "JND")
    
    # Bayesian Adjusted PSE and JND 
    p_samp_adj <- stan_params_adjusted(fit)
    
    tibble(iter = iter, method = "BH-GNM",
           fit_sample_mean_pse     = p_samp$pse,     fit_sample_mean_jnd     = p_samp$jnd,
           fit_sample_mean_pse_adj = p_samp_adj$pse, fit_sample_mean_jnd_adj = p_samp_adj$jnd,
           fit_pop_pse             = p_pop$pse,      fit_pop_jnd             = p_pop$jnd,
           t_pse_stat = NA_real_, t_pse_p = NA_real_,
           t_jnd_stat = NA_real_, t_jnd_p = NA_real_,
           SSE = stan_sse(fit, datistan$y, datistan$n),
           true_sample_pse = true_sample_pse, true_sample_jnd = true_sample_jnd,
           true_fixeff_pse = true_fixeff_pse_val, true_fixeff_jnd = true_fixeff_jnd_val)
  }, error = function(e) NULL)
}

# ============================================================
# Main Loop Execution
# ============================================================
all_methods <- if (run_stan) c("GLM", "GNM", "GLMM", "BH-GLM", "BH-GNM") else c("GLM", "GNM", "GLMM")

results_list <- setNames(vector("list", length(all_methods)), all_methods)
n_done       <- setNames(integer(length(all_methods)), all_methods)
iter_counter <- 0L

while (any(n_done < n_target)) {
  
  sim <- simulate_dataset()
  iter_counter <- iter_counter + 1L
  
  for (m in all_methods) {
    if (n_done[m] >= n_target) next
    row <- switch(m,
                  "GLM"    = fit_GLM( iter_counter, sim),
                  "GNM"    = fit_GNM( iter_counter, sim),
                  "GLMM"   = fit_GLMM(iter_counter, sim),
                  "BH-GLM" = fit_BHGLM(iter_counter, sim),
                  "BH-GNM" = fit_BHGNM(iter_counter, sim)
    )
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
    true_fixeff_pse_adj = true_fixeff_pse_adj,
    true_fixeff_jnd_adj = true_fixeff_jnd_adj,
    bias_samp_pse     = fit_sample_mean_pse - true_sample_pse,
    bias_samp_jnd     = fit_sample_mean_jnd - true_sample_jnd,
    bias_samp_pse_adj = fit_sample_mean_pse_adj - true_fixeff_pse_adj, # Added
    bias_samp_jnd_adj = fit_sample_mean_jnd_adj - true_fixeff_jnd_adj, # Added
    bias_pop_pse      = fit_pop_pse - true_fixeff_pse,
    bias_pop_jnd      = fit_pop_jnd - true_fixeff_jnd
  )

summary_table <- results %>%
  group_by(method) %>%
  summarise(
    # --- Sample Level Metrics (Simplified) ---
    bias_samp_pse       = mean(fit_sample_mean_pse - true_fixeff_pse_adj, na.rm = TRUE),
    rmse_samp_pse       = sqrt(mean((fit_sample_mean_pse - true_fixeff_pse_adj)^2, na.rm = TRUE)),
    bias_samp_jnd       = mean(fit_sample_mean_jnd - true_fixeff_jnd_adj, na.rm = TRUE),
    rmse_samp_jnd       = sqrt(mean((fit_sample_mean_jnd - true_fixeff_jnd_adj)^2, na.rm = TRUE)),
    
    # --- Sample Level Metrics (Adjusted for Lapse/Guess) ---
    bias_samp_pse_adj   = mean(fit_sample_mean_pse_adj - true_fixeff_pse_adj, na.rm = TRUE),
    rmse_samp_pse_adj   = sqrt(mean((fit_sample_mean_pse_adj - true_fixeff_pse_adj)^2, na.rm = TRUE)),
    bias_samp_jnd_adj   = mean(fit_sample_mean_jnd_adj - true_fixeff_jnd_adj, na.rm = TRUE),
    rmse_samp_jnd_adj   = sqrt(mean((fit_sample_mean_jnd_adj - true_fixeff_jnd_adj)^2, na.rm = TRUE)),
    
    # --- Population Level Metrics ---
    bias_pop_pse        = mean(fit_pop_pse - true_fixeff_pse_adj, na.rm = TRUE),
    rmse_pop_pse        = sqrt(mean((fit_pop_pse - true_fixeff_pse_adj)^2, na.rm = TRUE)),
    bias_pop_jnd        = mean(fit_pop_jnd - true_fixeff_jnd_adj, na.rm = TRUE),
    rmse_pop_jnd        = sqrt(mean((fit_pop_jnd - true_fixeff_jnd_adj)^2, na.rm = TRUE)),
    
    mean_SSE            = mean(SSE, na.rm = TRUE),
    reject_rate_pse     = mean(t_pse_p < 0.05, na.rm = TRUE),
    reject_rate_jnd     = mean(t_jnd_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

print("--- GLOBAL SUMMARY TABLE ---")
print(summary_table)

summary_table_pos <- results %>%
  dplyr::filter(fit_sample_mean_pse > 0) %>%
  group_by(method) %>%
  summarise(
    # --- Sample Level Metrics (Simplified) ---
    bias_samp_pse       = mean(fit_sample_mean_pse - true_sample_pse, na.rm = TRUE),
    rmse_samp_pse       = sqrt(mean((fit_sample_mean_pse - true_sample_pse)^2, na.rm = TRUE)),
    bias_samp_jnd       = mean(fit_sample_mean_jnd - true_sample_jnd, na.rm = TRUE),
    rmse_samp_jnd       = sqrt(mean((fit_sample_mean_jnd - true_sample_jnd)^2, na.rm = TRUE)),
    
    # --- Sample Level Metrics (Adjusted for Lapse/Guess) ---
    bias_samp_pse_adj   = mean(fit_sample_mean_pse_adj - true_sample_pse, na.rm = TRUE),
    rmse_samp_pse_adj   = sqrt(mean((fit_sample_mean_pse_adj - true_sample_pse)^2, na.rm = TRUE)),
    bias_samp_jnd_adj   = mean(fit_sample_mean_jnd_adj - true_sample_jnd, na.rm = TRUE),
    rmse_samp_jnd_adj   = sqrt(mean((fit_sample_mean_jnd_adj - true_sample_jnd)^2, na.rm = TRUE)),
    
    # --- Population Level Metrics ---
    bias_pop_pse        = mean(fit_pop_pse - true_fixeff_pse, na.rm = TRUE),
    rmse_pop_pse        = sqrt(mean((fit_pop_pse - true_fixeff_pse)^2, na.rm = TRUE)),
    bias_pop_jnd        = mean(fit_pop_jnd - true_fixeff_jnd, na.rm = TRUE),
    rmse_pop_jnd        = sqrt(mean((fit_pop_jnd - true_fixeff_jnd)^2, na.rm = TRUE)),
    
    mean_SSE            = mean(SSE, na.rm = TRUE),
    reject_rate_pse     = mean(t_pse_p < 0.05, na.rm = TRUE),
    reject_rate_jnd     = mean(t_jnd_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

print("--- POSITIVE PSE SUMMARY TABLE ---")
print(summary_table_pos)

# ============================================================
# Plotting
# ============================================================
violin_n150 <- list()
violin_n150_filename <- list("violin_150_pse.pdf", "violin_150_jnd.pdf")

violin_n150[["pse"]] <- ggplot(data = results, mapping = aes(y = bias_samp_pse, x = method)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "PSE sample bias", x = NULL) +
  coord_cartesian(ylim = c(-10, 10))+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")

violin_n150[["jnd"]] <- ggplot(data = results, mapping = aes(y = bias_samp_jnd, x = method)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "JND sample bias", x = NULL) +
  coord_cartesian(ylim = c(-3, 10)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")

map2(.x = violin_n150_filename, .y = violin_n150, .f = ggsave)

library(patchwork)
combined_plot <- violin_n150[["jnd"]] / violin_n150[["pse"]]
combined_plot <- combined_plot + plot_annotation(tag_levels = 'A')

ggsave(filename = "combined_violins_vertical.pdf", 
       plot = combined_plot, device = "pdf", width = 6, height = 9)