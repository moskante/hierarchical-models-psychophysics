# Monte Carlo simulation study: cross-method comparison of GLM, GNM, GLMM,
# BH-GLM, and BH-GNM on independently simulated psychophysical datasets.
#
#   Each iteration draws a fresh dataset, fits all five models, and records
#   estimated PSE and JND for each. The loop continues until exactly n_target
#   *successful* fits have been collected per method. If a method fails on a
#   given dataset (e.g. optimisation divergence, Stan non-convergence), the
#   outer while loop draws a new dataset and retries only the failing methods.
#   This ensures balanced comparison without biasing results by discarding
#   entire iterations when only one method fails.


library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(rstan)
library(gnlm)

source("../R/gnlm_functions_psejnd.R")   # defines process_subject() for GNM fitting

# ============================================================
# Monte Carlo loop: compare GLM, GNM, GLMM, BH-GLM, BH-GNM
# across n_iter independent simulated datasets.
# ============================================================

n_target  <- 150    # number of successful iterations required per method
ntrials   <- 160
nsubjects <- 10
run_stan  <- TRUE   # set FALSE for quick GLM/GNM/GLMM-only runs

# Compile Stan models once at the start to avoid repeated (slow) compilation.

if (run_stan) {
  stan_bhglm <- stan_model(file = "../Stan/simul_bhglm.stan")
  stan_bhgnm <- stan_model(file = "../Stan/simul_bhgnm.stan")
}

# ============================================================
# Helpers
# ============================================================

# stan_sse(): compute Sum of Squared Errors for a Stan model.
# Compares the median posterior predictive probability to the observed
# proportion correct, summed across all observations.
stan_sse <- function(stan_fit, y_obs, trials) {
  fitted_probs <- apply(extract(stan_fit)$predProb, 2, median)
  sum((y_obs / trials - fitted_probs)^2)
}

# stan_params(): extract population-mean PSE and JND from a Stan fit.
# Takes the per-subject posterior medians, then averages across subjects.
stan_params <- function(stan_fit, param_pse = "pse", param_jnd = "jnd") {
  s <- extract(stan_fit)
  list(pse = mean(apply(s[[param_pse]], 2, median), na.rm = TRUE),
       jnd = mean(apply(s[[param_jnd]], 2, median), na.rm = TRUE))
}

# ============================================================
# Simulate one dataset and extract true parameters
# ============================================================
simulate_dataset <- function() {
  simul_data         <- PsySimulate(ntrials   = ntrials,
                                    nsubjects = nsubjects,
                                    guess     = TRUE,
                                    lapse     = TRUE)
  simul_data$Subject <- factor(simul_data$Subject)
  # Recover the true per-subject PSE and JND from the generating parameters.
  params_true <- simul_data %>%
    group_by(Subject) %>%
    summarise(across(everything(), first),
              pse = -Intercept / Slope,
              jnd = qnorm(0.75) / Slope,
              .groups = "drop")
  list(data     = simul_data,
       true_PSE = mean(params_true$pse),
       true_JND = mean(params_true$jnd))
}

# ============================================================
# Fit functions: one per method, each returns a result row
# or NULL on failure. Each accepts a pre-simulated dataset.
# ============================================================

fit_GLM <- function(iter, simul_data, true_PSE, true_JND) {
  tryCatch({
    glm_list   <- PsychModels(cbind(Longer, Total - Longer) ~ X,
                              data = simul_data,
                              group_factors = "Subject")
    params_glm <- PsychParameters(glm_list)
    # One-sample t-tests against the true generating mean 
    t_pse      <- t.test(params_glm$pse,  mu = true_PSE)
    t_jnd      <- t.test(params_glm$jnd,  mu = true_JND)
    tibble(iter = iter, method = "GLM",
           mean_pse   = mean(params_glm$pse, na.rm = TRUE),
           mean_jnd   = mean(params_glm$jnd, na.rm = TRUE),
           t_pse_stat = as.numeric(t_pse$statistic), t_pse_p = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic), t_jnd_p = t_jnd$p.value,
           SSE = NA_real_, true_PSE = true_PSE, true_JND = true_JND)
  }, error = function(e) NULL)
}

fit_GNM <- function(iter, simul_data, true_PSE, true_JND) {
    tryCatch({
    params_gnm <- simul_data %>%
      group_split(Subject) %>%
      map_dfr(process_subject)
    
    t_pse <- t.test(params_gnm$pse, mu = true_PSE)
    t_jnd <- t.test(params_gnm$jnd, mu = true_JND)
    tibble(iter = iter, method = "GNM",
           mean_pse   = mean(params_gnm$pse, na.rm = TRUE),
           mean_jnd   = mean(params_gnm$jnd, na.rm = TRUE),
           t_pse_stat = as.numeric(t_pse$statistic), t_pse_p = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic), t_jnd_p = t_jnd$p.value,
           SSE = NA_real_, true_PSE = true_PSE, true_JND = true_JND)
  }, error = function(e) NULL)
}

fit_GLMM <- function(iter, simul_data, true_PSE, true_JND) {
  tryCatch({
    glmm_fit <- glmer(
      cbind(Longer, Total - Longer) ~ X + (1 + X | Subject),
      family = binomial(link = "probit"), data = simul_data
    )
    fe      <- fixef(glmm_fit)
    # Population-level PSE and JND derived from fixed-effect intercept and slope
    tibble(iter = iter, method = "GLMM",
           mean_pse   = as.numeric(-fe["(Intercept)"] / fe["X"]),
           mean_jnd   = as.numeric(qnorm(0.75) / fe["X"]),
           t_pse_stat = NA_real_, t_pse_p = NA_real_,
           t_jnd_stat = NA_real_, t_jnd_p = NA_real_,
           SSE        = sum(residuals(glmm_fit, type = "response")^2),
           true_PSE   = true_PSE, true_JND = true_JND)
  }, error = function(e) NULL)
}

fit_BHGLM <- function(iter, simul_data, true_PSE, true_JND, datistan) {
  tryCatch({
    # Initialise chains near the true PSE to improve mixing
    init_fn <- function() list(pse   = rep(true_PSE, datistan$nsubj),
                               sigma = rep(15,       datistan$nsubj))
    fit <- sampling(stan_bhglm, data = datistan,
                    chains = 3, warmup = 3000, iter = 5000,
                    cores = 3, refresh = 0,
                    init = list(init_fn(), init_fn(), init_fn()))
    p <- stan_params(fit, param_pse = "pse", param_jnd = "jnd")
    tibble(iter = iter, method = "BH-GLM",
           mean_pse = p$pse, mean_jnd = p$jnd,
           t_pse_stat = NA_real_, t_pse_p = NA_real_,
           t_jnd_stat = NA_real_, t_jnd_p = NA_real_,
           SSE = stan_sse(fit, datistan$y, datistan$n),
           true_PSE = true_PSE, true_JND = true_JND)
  }, error = function(e) NULL)
}

fit_BHGNM <- function(iter, simul_data, true_PSE, true_JND, datistan) {
  tryCatch({
    init_fn <- function() list(pse    = rep(true_PSE, datistan$nsubj),
                               sigma  = rep(15,       datistan$nsubj),
                               gamma  = rep(0.01,     datistan$nsubj),
                               lambda = rep(0.01,     datistan$nsubj))
    fit <- sampling(stan_bhgnm, data = datistan,
                    chains = 3, warmup = 3000, iter = 5000,
                    cores = 3, refresh = 0,
                    init = list(init_fn(), init_fn(), init_fn()))
    p <- stan_params(fit, param_pse = "pse", param_jnd = "jnd")
    tibble(iter = iter, method = "BH-GNM",
           mean_pse = p$pse, mean_jnd = p$jnd,
           t_pse_stat = NA_real_, t_pse_p = NA_real_,
           t_jnd_stat = NA_real_, t_jnd_p = NA_real_,
           SSE = stan_sse(fit, datistan$y, datistan$n),
           true_PSE = true_PSE, true_JND = true_JND)
  }, error = function(e) NULL)
}

# ============================================================
# Main loop: collect exactly n_target successes per method.
# On any method failure, resimulate and retry only that method.
# ============================================================
all_methods <- if (run_stan) c("GLM", "GNM", "GLMM", "BH-GLM", "BH-GNM") else c("GLM", "GNM", "GLMM")

# Initialise per-method result lists and counters
results_list <- setNames(vector("list", length(all_methods)), all_methods)
n_done       <- setNames(integer(length(all_methods)), all_methods)
iter_counter <- 0L   # global attempt counter for labelling

while (any(n_done < n_target)) {
  
  # Draw a fresh dataset 
  sim        <- simulate_dataset()
  simul_data <- sim$data
  true_PSE   <- sim$true_PSE
  true_JND   <- sim$true_JND
  
  # Stan data list (built once per dataset, reused by both Stan models)
  if (run_stan) {
    datistan <- list(
      y       = simul_data$Longer,
      n       = simul_data$Total,
      nobs    = nrow(simul_data),
      x       = simul_data$X,
      subject = as.integer(simul_data$Subject),
      nsubj   = nlevels(simul_data$Subject)
    )
  }
  
  iter_counter <- iter_counter + 1L
  
  # Try each method that still needs more successes
  for (m in all_methods) {
    if (n_done[m] >= n_target) next   # already complete, skip
    
    row <- switch(m,
                  "GLM"    = fit_GLM( iter_counter, simul_data, true_PSE, true_JND),
                  "GNM"    = fit_GNM( iter_counter, simul_data, true_PSE, true_JND),
                  "GLMM"   = fit_GLMM(iter_counter, simul_data, true_PSE, true_JND),
                  "BH-GLM" = fit_BHGLM(iter_counter, simul_data, true_PSE, true_JND, datistan),
                  "BH-GNM" = fit_BHGNM(iter_counter, simul_data, true_PSE, true_JND, datistan)
    )
    
    if (!is.null(row)) {
      n_done[m] <- n_done[m] + 1L
      results_list[[m]][[n_done[m]]] <- row
    }
    # If NULL: method failed on this dataset; the outer while loop
    # will draw a new dataset on the next iteration and retry.
  }
  
  message(sprintf("Attempt %d | successes: %s",
                  iter_counter,
                  paste(names(n_done), n_done, sep = "=", collapse = " | ")))
}

results <- bind_rows(lapply(results_list, bind_rows)) %>%
  mutate(
    bias_pse = mean_pse - true_PSE,
    bias_jnd = mean_jnd - true_JND
  )

# ============================================================
# Summary table: objective cross-method comparison
# ============================================================
summary_table <- results %>%
  group_by(method) %>%
  summarise(
    # Bias and RMSE relative to the true generating values
    bias_pse        = mean(mean_pse - true_PSE, na.rm = TRUE),
    rmse_pse        = sqrt(mean((mean_pse - true_PSE)^2, na.rm = TRUE)),
    bias_jnd        = mean(mean_jnd - true_JND, na.rm = TRUE),
    rmse_jnd        = sqrt(mean((mean_jnd - true_JND)^2, na.rm = TRUE)),
    # Mean SSE (available for GLMM, BH-GLM, BH-GNM; NA for GLM/GNM)
    mean_SSE        = mean(SSE, na.rm = TRUE),
    # Empirical Type I error rate: proportion of iterations where H0 was
    # rejected at alpha = 0.05 (should be ~0.05 under a correct model)
    reject_rate_pse = mean(t_pse_p < 0.05, na.rm = TRUE),
    reject_rate_jnd = mean(t_jnd_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

# Plotting results 

violin_n150 <- list()
violin_n150_filename <- list("violin_150_pse.pdf", "violin_150_jnd.pdf")

violin_n150[["pse"]] <- ggplot(data = results, mapping = aes(y = bias_pse, x = method)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "PSE bias", x = NULL) +
  coord_cartesian(ylim = c(-10, 10))+
  geom_hline( yintercept = 0,color = "red",linetype = "dashed")

violin_n150[["jnd"]] <- ggplot(data = results, mapping = aes(y = bias_jnd, x = method)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "JND bias", x = NULL) +
  coord_cartesian(ylim = c(-3, 10)) +
  geom_hline( yintercept = 0,color = "red",linetype = "dashed")

map2(.x = violin_n150_filename, .y = violin_n150, .f = ggsave)
