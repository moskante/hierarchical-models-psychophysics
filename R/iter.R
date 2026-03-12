library(tidyverse)
library(MixedPsy)
library(gnlm)
library(lme4)
library(lmerTest)
library(rstan)

source("R/gnlm_functions_psejnd.R")  # provides mu, pstart, PsychParametersGNM, process_subject

# ============================================================
# Monte Carlo loop: compare GLM, GNM, GLMM, BH-GLM, BH-GNM
# across n independent simulated datasets.
#
# Motivation: R2 and the Editor note that cross-method
# comparisons based on a single simulated dataset are
# insufficient. This script replicates the analysis n_iter
# times with independent draws (no set.seed) and collects
# objective per-iteration summaries for each method.
# ============================================================

n_iter    <- 150    # Monte Carlo replications
# Stan models are slow: consider n_iter = 20-50
# and/or set run_stan = FALSE during development.
ntrials   <- 160
nsubjects <- 10

# Set run_stan = FALSE to skip BH-GLM and BH-GNM during
# quick checks.
run_stan  <- TRUE

# ---- Stan: compile once outside the loop -------------------
# sampling() reuses the compiled object across all iterations,
# avoiding expensive recompilation at each step.
if (run_stan) {
  stan_bhglm <- stan_model(file = "Stan/simul_bhglm.stan")
  stan_bhgnm <- stan_model(file = "Stan/simul_bhgnm.stan")
}

# ============================================================
# Helper: SSE from a Stan posterior predictive distribution
# ============================================================
stan_sse <- function(stan_fit, y_obs, trials) {
  samples      <- extract(stan_fit)
  fitted_probs <- apply(samples$predProb, 2, median)
  sum((y_obs / trials - fitted_probs)^2)
}

# ============================================================
# Helper: population-average PSE / JND from a Stan fit.
# Computes the median of each subject's posterior median,
# then averages across subjects.
# Adapt param_pse / param_jnd to match your .stan files.
# ============================================================
stan_params <- function(stan_fit, param_pse = "pse",
                        param_jnd  = "jnd") {
  samples  <- extract(stan_fit)
  pse_subj <- apply(samples[[param_pse]], 2, median)
  jnd_subj <- apply(samples[[param_jnd]], 2, median)
  list(pse = mean(pse_subj, na.rm = TRUE),
       jnd = mean(jnd_subj, na.rm = TRUE))
}

# ============================================================
# Core: run every method on one simulated dataset
# ============================================================
run_one_iteration <- function(iter) {
  
  message(sprintf("--- Iteration %d / %d ---", iter, n_iter))
  
  # ---- 1. Simulate (no set.seed: each call is independent) --
  simul_data         <- PsySimulate(ntrials   = ntrials,
                                    nsubjects = nsubjects,
                                    guess     = TRUE,
                                    lapse     = TRUE)
  simul_data$Subject <- factor(simul_data$Subject)
  
  # ---- 2. True generating parameters -----------------------
  parameters_simul <- simul_data %>%
    group_by(Subject) %>%
    summarise(across(everything(), first),
              pse    = -Intercept / Slope,
              jnd    = qnorm(0.75) / Slope,
              gamma  = Gamma,
              lambda = Lambda,
              .groups = "drop") %>%
    select(Subject, pse, jnd, gamma, lambda)
  
  true_PSE <- mean(parameters_simul$pse)
  true_JND <- mean(parameters_simul$jnd)
  
  rows <- list()
  
  # ---- 3. GLM -----------------------------------------------
  rows[["GLM"]] <- tryCatch({
    glm_list       <- PsychModels(formula = cbind(Longer, Total - Longer) ~ X,
                                  data    = simul_data,
                                  group_factors = "Subject")
    params_glm     <- PsychParameters(glm_list)
    t_pse          <- t.test(params_glm$pse, mu = true_PSE)
    t_jnd          <- t.test(params_glm$jnd, mu = true_JND)
    tibble(iter = iter, method = "GLM",
           mean_pse   = mean(params_glm$pse, na.rm = TRUE),
           mean_jnd   = mean(params_glm$jnd, na.rm = TRUE),
           t_pse_stat = as.numeric(t_pse$statistic),
           t_pse_p    = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic),
           t_jnd_p    = t_jnd$p.value,
           SSE        = NA_real_,
           true_PSE   = true_PSE, true_JND = true_JND,
           converged  = TRUE)
  }, error = function(e)
    tibble(iter = iter, method = "GLM",
           mean_pse = NA, mean_jnd = NA,
           t_pse_stat = NA, t_pse_p = NA,
           t_jnd_stat = NA, t_jnd_p = NA,
           SSE = NA_real_,
           true_PSE = true_PSE, true_JND = true_JND,
           converged = FALSE))
  
  # ---- 4. GNM -----------------------------------------------
  rows[["GNM"]] <- tryCatch({
    gnm_params <- simul_data %>%
      group_split(Subject) %>%
      map_dfr(process_subject)
    t_pse <- t.test(gnm_params$pse, mu = true_PSE)
    t_jnd <- t.test(gnm_params$jnd, mu = true_JND)
    tibble(iter = iter, method = "GNM",
           mean_pse   = mean(gnm_params$pse, na.rm = TRUE),
           mean_jnd   = mean(gnm_params$jnd, na.rm = TRUE),
           t_pse_stat = as.numeric(t_pse$statistic),
           t_pse_p    = t_pse$p.value,
           t_jnd_stat = as.numeric(t_jnd$statistic),
           t_jnd_p    = t_jnd$p.value,
           SSE        = NA_real_,
           true_PSE   = true_PSE, true_JND = true_JND,
           converged  = TRUE)
  }, error = function(e)
    tibble(iter = iter, method = "GNM",
           mean_pse = NA, mean_jnd = NA,
           t_pse_stat = NA, t_pse_p = NA,
           t_jnd_stat = NA, t_jnd_p = NA,
           SSE = NA_real_,
           true_PSE = true_PSE, true_JND = true_JND,
           converged = FALSE))
  
  # ---- 5. GLMM ----------------------------------------------
  # The GLMM yields a single population-level estimate (fixef),
  # so per-subject t-tests are not applicable here. Bias and
  # RMSE are computed relative to the mean true parameter.
  rows[["GLMM"]] <- tryCatch({
    glmm_fit   <- glmer(
      cbind(Longer, Total - Longer) ~ X + (1 + X | Subject),
      family = binomial(link = "probit"),
      data   = simul_data
    )
    resid_glmm <- residuals(glmm_fit, type = "response")
    sse_glmm   <- sum(resid_glmm^2)
    fe         <- fixef(glmm_fit)
    pop_pse    <- -fe["(Intercept)"] / fe["X"]
    pop_jnd    <- qnorm(0.75) / fe["X"]
    tibble(iter = iter, method = "GLMM",
           mean_pse   = as.numeric(pop_pse),
           mean_jnd   = as.numeric(pop_jnd),
           t_pse_stat = NA_real_,   # single estimate; no cross-subject t-test
           t_pse_p    = NA_real_,
           t_jnd_stat = NA_real_,
           t_jnd_p    = NA_real_,
           SSE        = sse_glmm,
           true_PSE   = true_PSE, true_JND = true_JND,
           converged  = TRUE)
  }, error = function(e)
    tibble(iter = iter, method = "GLMM",
           mean_pse = NA, mean_jnd = NA,
           t_pse_stat = NA, t_pse_p = NA,
           t_jnd_stat = NA, t_jnd_p = NA,
           SSE = NA_real_,
           true_PSE = true_PSE, true_JND = true_JND,
           converged = FALSE))
  
  # ---- 6. Stan data list (shared by BH-GLM and BH-GNM) -----
  if (run_stan) {
    
    datistan <- list(
      y       = simul_data$Longer,
      n       = simul_data$Total,
      nobs    = nrow(simul_data),
      x       = simul_data$X,
      subject = as.integer(simul_data$Subject),
      nsubj   = nlevels(simul_data$Subject)
    )
    
    # ---- 7. BH-GLM ------------------------------------------
    rows[["BH-GLM"]] <- tryCatch({
      init_bhglm <- function() list(
        pse   = rep(true_PSE, datistan$nsubj),
        sigma = rep(15,        datistan$nsubj)
      )
      fit_bhglm <- sampling(
        stan_bhglm,  
        data    = datistan,
        chains  = 3, warmup = 3000, iter = 5000, cores = 3,
        refresh = 0,
        init    = list(init_bhglm(), init_bhglm(), init_bhglm())
      )
      sse_bhglm <- stan_sse(fit_bhglm, datistan$y, datistan$n)
      params    <- stan_params(fit_bhglm, param_pse = "pse",
                               param_jnd = "jnd")
      tibble(iter = iter, method = "BH-GLM",
             mean_pse   = params$pse,
             mean_jnd   = params$jnd,
             t_pse_stat = NA_real_,
             t_pse_p    = NA_real_,
             t_jnd_stat = NA_real_,
             t_jnd_p    = NA_real_,
             SSE        = sse_bhglm,
             true_PSE   = true_PSE, true_JND = true_JND,
             converged  = TRUE)
    }, error = function(e)
      tibble(iter = iter, method = "BH-GLM",
             mean_pse = NA, mean_jnd = NA,
             t_pse_stat = NA, t_pse_p = NA,
             t_jnd_stat = NA, t_jnd_p = NA,
             SSE = NA_real_,
             true_PSE = true_PSE, true_JND = true_JND,
             converged = FALSE))
    
    # ---- 8. BH-GNM ------------------------------------------
    rows[["BH-GNM"]] <- tryCatch({
      init_bhgnm <- function() list(
        pse    = rep(true_PSE, datistan$nsubj),
        sigma  = rep(15,        datistan$nsubj),
        gamma  = rep(0.01,      datistan$nsubj),
        lambda = rep(0.01,      datistan$nsubj)
      )
      fit_bhgnm <- sampling(
        stan_bhgnm,
        data    = datistan,
        chains  = 3, warmup = 3000, iter = 5000, cores = 3,
        refresh = 0,
        init    = list(init_bhgnm(), init_bhgnm(), init_bhgnm())
      )
      sse_bhgnm <- stan_sse(fit_bhgnm, datistan$y, datistan$n)
      params    <- stan_params(fit_bhgnm, param_pse = "pse",
                               param_jnd = "jnd")
      tibble(iter = iter, method = "BH-GNM",
             mean_pse   = params$pse,
             mean_jnd   = params$jnd,
             t_pse_stat = NA_real_,
             t_pse_p    = NA_real_,
             t_jnd_stat = NA_real_,
             t_jnd_p    = NA_real_,
             SSE        = sse_bhgnm,
             true_PSE   = true_PSE, true_JND = true_JND,
             converged  = TRUE)
    }, error = function(e)
      tibble(iter = iter, method = "BH-GNM",
             mean_pse = NA, mean_jnd = NA,
             t_pse_stat = NA, t_pse_p = NA,
             t_jnd_stat = NA, t_jnd_p = NA,
             SSE = NA_real_,
             true_PSE = true_PSE, true_JND = true_JND,
             converged = FALSE))
  }
  
  bind_rows(rows)
}

# ============================================================
# Run the loop
# ============================================================
safe_run <- possibly(run_one_iteration, otherwise = NULL)

results <- map(seq_len(n_iter), safe_run) %>%
  compact() %>%
  bind_rows()

# ============================================================
# Summary table: objective cross-method comparison
# ============================================================
summary_table <- results %>%
  group_by(method) %>%
  summarise(
    n_converged     = sum(converged, na.rm = TRUE),
    # Bias and RMSE relative to the true generating values
    bias_pse        = mean(mean_pse - true_PSE, na.rm = TRUE),
    rmse_pse        = sqrt(mean((mean_pse - true_PSE)^2, na.rm = TRUE)),
    bias_jnd        = mean(mean_jnd - true_JND, na.rm = TRUE),
    rmse_jnd        = sqrt(mean((mean_jnd - true_JND)^2, na.rm = TRUE)),
    # Mean SSE (available for GLMM, BH-GLM, BH-GNM)
    mean_SSE        = mean(SSE, na.rm = TRUE),
    # For GLM / GNM: empirical Type I error rate (should be ~0.05
    # when the model is correctly specified)
    reject_rate_pse = mean(t_pse_p < 0.05, na.rm = TRUE),
    reject_rate_jnd = mean(t_jnd_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

# ============================================================
# Save outputs
# ============================================================
write_csv(results,       "simulation_loop_results.csv")
write_csv(summary_table, "simulation_loop_summary.csv")

# ============================================================
# Plots
# ============================================================

# Violin plots of estimated PSE and JND per model across simulation iterations.
# Horizontal reference line = population parameters

ref_PSE <- 80
ref_JND <- 7.7

# ---- Method factor order (simple -> hierarchical) -------------
method_order <- c("GLM", "GNM", "GLMM", "BH-GLM", "BH-GNM")

results <- results %>%
  filter(converged) %>%
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

fig <- (p_pse | p_jnd) +
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
ggsave("figure_violin_pse_jnd.pdf", fig, width = 10, height = 5)
#ggsave("figure_violin_pse_jnd.png", fig, width = 10, height = 5, dpi = 300)

message("Saved: figure_violin_pse_jnd.pdf / .png")

# ================================================================
# Estimation error density figure
# ================================================================
# Panel order: PSE first, JND second — matching violin figure.
# Uses the same palette, theme_classic base, and strip style.

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
#ggsave("figure_density_error.png", fig_density, width = 10, height = 5, dpi = 300)

message("Saved: figure_density_error.pdf / .png")






# ================================================================
# SSE distribution figure (hierarchical models only)
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

ggsave("figure_sse_boxplot.pdf", fig_sse, width = 6, height = 5)
#ggsave("figure_sse_boxplot.png", fig_sse, width = 6, height = 5, dpi = 300)

message("Saved: figure_sse_boxplot.pdf / .png")
