# vibro_single_sub.R
#
# Example 2a: single-subject GLM and GNM analysis on vibrotactile data.
#
# Experiment: Nine subjects judged which of two moving tactile stimuli was
# faster, under two vibration conditions:
#   vibration == 0  -- no added vibration (control)
#   vibration == 32 -- 32 Hz vibration superimposed on the standard stimulus
#
# Research question: Does 32 Hz vibration shift the psychometric function
# slope (discrimination sensitivity)?
#
# Workflow:
#   1. Fit a probit GLM per subject x vibration condition.
#   2. Extract PSE, JND, and slope; test the vibration effect on slope.
#   3. Fit a GNM per subject x condition (adds guessing and lapsing).
#   4. Test the vibration effect on slope from the GNM.
#
# Required packages: tidyverse, MixedPsy, gnlm
# Source file:       R/gnlm_functions_slope.R
# Dataset:           vibro_exp3 (loaded from MixedPsy)

library(tidyverse)
library(MixedPsy)
library(gnlm)

# vibro_exp3 is bundled with the MixedPsy package.
# Columns: subject, vibration (factor: "0" or "32"), speed, faster, slower.
# faster + slower = total number of trials per row.

# ---------------------------------------------------------------------------
# GLM: probit model per subject x vibration condition
#
# PsychModels() fits a probit GLM for each combination of subject and
# vibration level. The response is (faster, slower) -- "faster" responses
# are the successes. PsychParameters() extracts PSE and JND from each fit.
# The slope coefficient is extracted manually from each model object.
# ---------------------------------------------------------------------------

# GLM ----
model_list_vibro <- PsychModels(vibro_exp3,
                                group_factors = c("subject", "vibration"),
                                formula       = cbind(faster, slower) ~ speed)
parameters_glm <- PsychParameters(model_list_vibro)

# Append slope and its SE to the parameters table
parameters_glm <- parameters_glm %>%
  mutate(
    slope    = map_dbl(model_list_vibro, ~ .x$model$coefficients["speed"]),
    slope_se = map_dbl(model_list_vibro, ~ summary(.x$model)$coefficients["speed", "Std. Error"])
  )

# Reshape to wide format for paired comparison:
# each row = one subject, columns for condition 0 and condition 32
## inference tests -----
params_glm_wider <- parameters_glm %>%
  select(-pse_se, -jnd_se, -slope_se) %>%
  pivot_wider(names_from = vibration, values_from = c(pse, jnd, slope))%>%
  mutate(slope_diff = slope_32 - slope_0)

# Paired t-test: does the slope differ between 32 Hz and 0 Hz vibration?
t_glm_slope <- t.test(params_glm_wider$slope_32, params_glm_wider$slope_0, paired = TRUE)


# ---------------------------------------------------------------------------
# GNM: nonlinear model with guessing and lapsing (slope parameterisation)
#
# gnlm_functions_slope.R defines process_subject_vibro(), which fits a GNM
# to one subject x condition split. The output includes PSE, slope, gamma,
# lambda, their bootstrap SEs, and 95% CIs.
# JND is derived post-hoc from the slope as qnorm(0.75) / slope.
# ---------------------------------------------------------------------------

# GNM ----
source("R/gnlm_functions_slope.R")
parameters_gnm <- vibro_exp3 %>%
  group_split(subject, vibration) %>%
  map_dfr(process_subject_vibro)

# Derive JND from the GNM slope estimate
parameters_gnm$jnd <- qnorm(0.75) / parameters_gnm$slope

# Wide format for paired comparison
params_gnm_wider <- parameters_gnm %>%
  select(pse, jnd, slope, vibration, subject) %>%
  pivot_wider(names_from = vibration, values_from = c(pse, jnd, slope))%>%
  mutate(slope_diff = slope_32 - slope_0)

# Paired t-test on GNM slopes
t_gnm_slope <- t.test(params_gnm_wider$slope_32, params_gnm_wider$slope_0, paired = TRUE)

# model plot JND and PSE -------------------------------------------------------

# Helper function to compute mean and 95% CI from a numeric vector
calc_vec_ci <- function(x, param_name, model_name) {
  tt <- t.test(x)
  tibble(
    Parameter = param_name,
    Estimate  = mean(x, na.rm = TRUE),
    Lower     = tt$conf.int[1],
    Upper     = tt$conf.int[2],
    Model     = model_name
  )
}

# 1. Compute summary statistics and CIs for GLM and GNM vectors
df_combined <- bind_rows(
  calc_vec_ci(params_glm_wider$slope_diff, "Slope Diff", "GLM"),
  calc_vec_ci(params_gnm_wider$slope_diff, "Slope Diff", "GNM")
) %>%
  mutate(Model = factor(Model, levels = c("GLM", "GNM")))

# 2. Define facet-specific sample mean reference values
df_ref <- tibble(
  Parameter  = c("Slope Diff"),
  yintercept = 0
)

# 3. Faceted Comparison Plot
ggplot(df_combined, aes(x = Model, y = Estimate, color = Model)) +
  # Sample mean reference lines per facet
  geom_hline(
    data = df_ref,
    aes(yintercept = yintercept),
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.7
  ) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper), size = 0.8, linewidth = 1) +
  scale_color_manual(values = c(
    "GLM" = "#E41A1C",  # Red
    "GNM" = "#377EB8"   # Blue
  )) +
  labs(
    # title = "GLM vs. GNM Vector Estimates",
    # subtitle = "Points represent sample means with 95% t-distribution CIs; dashed lines show true sample means",
    x = NULL,
    y = "Estimate (32 Hz - 0 Hz)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray92"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold")
  )

ggsave("example_2__twolevels_estimates_diffs.pdf", path = "Figs",width = 6, height = 4)
