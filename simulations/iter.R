# iter.R
#
# Monte Carlo simulation study: cross-method comparison of GLM, GNM, GLMM,
# BH-GLM, and BH-GNM on independently simulated psychophysical datasets.
#
# Design:
#   Each iteration draws a fresh dataset, fits all five models, and records
#   estimated PSE and JND for each. The loop continues until exactly n_target
#   *successful* fits have been collected per method. If a method fails on a
#   given dataset (e.g. optimisation divergence, Stan non-convergence), the
#   outer while loop draws a new dataset and retries only the failing methods.
#   This ensures balanced comparison without biasing results by discarding
#   entire iterations when only one method fails.
#
# Outputs:
#   simulation_loop_results.csv  -- per-iteration estimates (one row per method)
#   simulation_loop_summary.csv  -- bias, RMSE, SSE, Type I error per method
#   figure_violin_pse_jnd.pdf    -- parameter recovery distributions
#   figure_density_error.pdf     -- estimation error densities
#   figure_sse_boxplot.pdf       -- SSE distributions (hierarchical models only)
#
# Required packages: tidyverse, MixedPsy, lme4, lmerTest, rstan, patchwork
# Source file:       R/gnlm_functions_psejnd.R

library(tidyverse)
library(MixedPsy)
library(lme4)
library(lmerTest)
library(rstan)

source("R/gnlm_functions_psejnd.R")   # defines process_subject() for GNM fitting

# ============================================================
# Monte Carlo loop: compare GLM, GNM, GLMM, BH-GLM, BH-GNM
# across n_iter independent simulated datasets.
# ============================================================

n_target  <- 3    # exact number of successful iterations required per method
ntrials   <- 160
nsubjects <- 10
run_stan  <- TRUE   # set FALSE for quick GLM/GNM/GLMM-only runs

# Compile Stan models once at the start; the compiled objects are reused
# across all iterations to avoid repeated (slow) compilation.
if (run_stan) {
  stan_bhglm <- stan_model(file = "Stan/simul_bhglm.stan")
  stan_bhgnm <- stan_model(file = "Stan/simul_bhgnm.stan")
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
  # Recover the true per-subject PSE and JND from the generating parameters
  # stored by PsySimulate() in each row (first row per subject is sufficient).
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
    # One-sample t-tests against the true generating mean (Type I error assessment)
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
  
  # Draw a fresh dataset for this attempt
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

results <- bind_rows(lapply(results_list, bind_rows))

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

# ============================================================
# Save outputs
# ============================================================
#write_csv(results,       "simulation_loop_results.csv")
#write_csv(summary_table, "simulation_loop_summary.csv")

# ============================================================
# Plots
# ============================================================

# Violin plots of estimated PSE and JND per model across simulation iterations.
# Horizontal reference line = population parameters used for data generation.

ref_PSE <- 80
ref_JND <- 7.7

# ---- Method factor order (simple -> hierarchical) -------------
method_order <- c("GLM", "GNM", "GLMM", "BH-GLM", "BH-GNM")

results <- results %>%
  mutate(method = factor(method, levels = method_order))

# ---- Palette ------------------------------------------------------
method_colours <- c(
  "GLM"    = "#4E79A7",
  "GNM"    = "#F28E2B",
  "GLMM"   = "#59A14F",
  "BH-GLM" = "#B07AA1",
  "BH-GNM" = "#E15759"
)

# ---- Shared theme -------------------------------------------------
violin_theme <- theme_classic(base_size = 12) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(size = 13, face = "bold"),
    axis.title.x      = element_blank(),
    legend.position   = "none",
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4)
  )

# ================================================================
# Panel A – PSE
# ================================================================
p_pse <- ggplot(results, aes(x = method, y = mean_pse, fill = method)) +
  geom_hline(yintercept = ref_PSE,
             linetype = "dashed", linewidth = 0.7, colour = "grey30") +
  geom_violin(trim = FALSE, alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               colour = "grey20", fill = "white", linewidth = 0.5) +
  annotate("text",
           x = 0.6, y = ref_PSE,
           label = sprintf("True mean = %.1f", ref_PSE),
           hjust = 0, vjust = -0.5, size = 3.3, colour = "grey30") +
  scale_fill_manual(values = method_colours) +
  scale_y_continuous(name = "Estimated PSE") +
  violin_theme

# ================================================================
# Panel B – JND
# ================================================================
p_jnd <- ggplot(results, aes(x = method, y = mean_jnd, fill = method)) +
  geom_hline(yintercept = ref_JND,
             linetype = "dashed", linewidth = 0.7, colour = "grey30") +
  geom_violin(trim = FALSE, alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               colour = "grey20", fill = "white", linewidth = 0.5) +
  annotate("text",
           x = 0.6, y = ref_JND,
           label = sprintf("True mean = %.2f", ref_JND),
           hjust = 0, vjust = -0.5, size = 3.3, colour = "grey30") +
  scale_fill_manual(values = method_colours) +
  scale_y_continuous(name = "Estimated JND") +
  violin_theme

# ================================================================
# Combine with patchwork
# ================================================================
library(patchwork)

fig_violin <- (p_pse | p_jnd) +
  plot_annotation(
    title   = "Parameter recovery across simulated datasets",
    subtitle = sprintf("n = %d iterations  ·  dashed line = true generating mean",
                       floor(nrow(results) / length(unique(results$method)))),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )

# ================================================================
# Save violin figure
# ================================================================
ggsave("figure_violin_pse_jnd.pdf", fig_violin, width = 10, height = 5)

message("Saved: figure_violin_pse_jnd.pdf / .png")

# ================================================================
# Estimation error density figure
#
# error = estimated parameter - true generating parameter.
# A well-calibrated model should have errors centred on zero.
# Wider distributions indicate lower precision.
# ================================================================
error_data <- results %>%
  mutate(
    PSE = mean_pse - true_PSE,
    JND = mean_jnd - true_JND
  ) %>%
  pivot_longer(c(PSE, JND),
               names_to  = "parameter",
               values_to = "error") %>%
  mutate(parameter = factor(parameter, levels = c("PSE", "JND")))

fig_density <- ggplot(error_data,
                      aes(x = error, fill = method, colour = method)) +
  geom_vline(xintercept = 0,
             linetype = "dashed", linewidth = 0.7, colour = "grey30") +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  facet_wrap(~ parameter, scales = "free", nrow = 1) +
  scale_fill_manual(values   = method_colours,
                    breaks   = method_order,
                    name     = "Model") +
  scale_colour_manual(values = method_colours,
                      breaks = method_order,
                      name   = "Model") +
  labs(
    title    = "Estimation error across simulated datasets",
    subtitle = sprintf("n = %d iterations  ·  dashed line = zero error",
                       floor(nrow(results) / length(unique(results$method)))),
    x = "Estimated \u2212 True",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(size = 13, face = "bold"),
    legend.position    = "right",
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
    plot.title         = element_text(size = 14, face = "bold"),
    plot.subtitle      = element_text(size = 10, colour = "grey40")
  )

ggsave("figure_density_error.pdf", fig_density, width = 10, height = 5)

message("Saved: figure_density_error.pdf / .png")


# ================================================================
# SSE distribution figure (hierarchical models only)
#
# SSE is defined only for GLMM, BH-GLM, and BH-GNM; GLM and GNM are
# excluded because they operate at the single-subject level and do not
# produce a single model-level SSE.
# ================================================================
fig_sse <- results %>%
  filter(!is.na(SSE)) %>%
  ggplot(aes(x = method, y = SSE, fill = method)) +
  geom_boxplot(alpha = 0.75, outlier.shape = 21,
               outlier.fill = NA, colour = "grey20", linewidth = 0.5) +
  scale_fill_manual(values = method_colours) +
  scale_x_discrete(limits = intersect(method_order,
                                      unique(results$method[!is.na(results$SSE)]))) +
  labs(
    title    = "SSE distribution across simulated datasets",
    subtitle = sprintf("n = %d iterations",
                       floor(nrow(results) / length(unique(results$method)))),
    x = NULL,
    y = "SSE"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
    plot.title         = element_text(size = 14, face = "bold"),
    plot.subtitle      = element_text(size = 10, colour = "grey40")
  )

#ggsave("figure_sse_boxplot.pdf", fig_sse, width = 6, height = 5)

#message("Saved: figure_sse_boxplot.pdf / .png")
