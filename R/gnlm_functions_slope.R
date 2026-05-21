# gnlm_functions_slope.R
#
# GNM helper functions for the vibrotactile analysis (Example 2).
#
# This file is the slope-parameterised counterpart of gnlm_functions_psejnd.R.
# The key difference is that here the psychometric function is written in the
# linear predictor form (beta_0 + beta_1 * speed) rather than the location-scale
# form (mu, sigma), because the primary outcome of interest is the *slope*
# (beta_1), which directly indexes discrimination sensitivity.
#
# Relationship between parameterisations:
#   mu    (PSE)   = -beta_0 / beta_1
#   sigma         = 1 / beta_1
#   JND           = qnorm(0.75) / beta_1
#
# As in gnlm_functions_psejnd.R, gamma (guessing) and lambda (lapsing) are
# estimated as free parameters on an unconstrained scale and back-transformed
# via atan2.

# ---------------------------------------------------------------------------
# mu() -- psychometric function in linear-predictor form
#
# P(speed) = gamma + (1 - gamma - lambda) * Phi(beta_0 + beta_1 * speed)
# ---------------------------------------------------------------------------
mu <- function(p){
  beta_0  <- p[1]   # intercept of the linear predictor
  beta_1  <- p[2]   # slope: positive values -> increasing psychometric function
  # Back-transform guessing and lapse parameters from unconstrained scale
  gamma_i  <- atan2(p[3], 1) / pi + 0.5
  lambda_i <- (1 - gamma_i) * (atan2(p[4], 1) / pi + 0.5)
  gamma_i + (1 - gamma_i - lambda_i) * pnorm(beta_0 + beta_1 * x)
}

# ---------------------------------------------------------------------------
# pstart() -- generate starting values for gnlr() from a probit GLM
# ---------------------------------------------------------------------------
pstart <- function(x, y){
  glm_start  <- glm(y ~ x, family = binomial(link = "probit"))
  coef_start <- coef(glm_start)
  prop_correct   <- y[, 1] / rowSums(y)
  gamma_i_start  <- min(min(prop_correct, na.rm = TRUE), 0.99)
  lambda_i_start <- max(1 - max(prop_correct, na.rm = TRUE), 0.01)
  beta_0_start   <- coef_start[1]   # GLM intercept -> starting beta_0
  beta_1_start   <- coef_start[2]   # GLM slope     -> starting beta_1
  p3_start <- tan(pi * (gamma_i_start - 0.5))
  p4_start <- tan(pi * ((lambda_i_start / (1 - gamma_i_start)) - 0.5))
  return(c(beta_0_start, beta_1_start, p3_start, p4_start))
}

# ---------------------------------------------------------------------------
# PsychParametersGNM() -- extract PSE, slope, gamma, lambda from a fitted
# gnlr object (slope parameterisation)
#
# PSE is recovered as -beta_0 / beta_1 (the speed at which P = 0.5 in the
# absence of guessing/lapsing).
# ---------------------------------------------------------------------------
PsychParametersGNM <- function(gnm, p = 0.75){
  pout    <- gnm$coefficients
  pse     <- -pout[1] / pout[2]                       # threshold (50% point)
  slope   <- pout[2]                                   # sensitivity index
  gamma_i  <- atan2(pout[3], 1) / pi + 0.5
  lambda_i <- (1 - gamma_i) * (atan2(pout[4], 1) / pi + 0.5)
  return(c(pse = pse, slope = slope, gamma = gamma_i, lambda = lambda_i))
}

# ---------------------------------------------------------------------------
# PsychBootGNM_vibro() -- bootstrap statistic for the vibrotactile GNM
#
# Differs from PsychBootGNM in gnlm_functions_psejnd.R only in using the
# vibrotactile column names (faster, slower, speed).
# ---------------------------------------------------------------------------
PsychBootGNM_vibro <- function(data, indices, p = 0.75){
  data_resampled <- data[indices, ]
  y   <- with(data_resampled, cbind(faster, slower))
  x   <- assign("x", data_resampled$speed, envir = .GlobalEnv)
  pmu <- pstart(x, y)
  gnm <- tryCatch(
    {
      model <- gnlr(y = y, distribution = "binomial",
                    mu = mu, pmu = pmu, iterlim = 10000)
      PsychParametersGNM(model, p)
    },
    error = function(e) rep(NA, 4)
  )
  # Restore x to full subject data
  x <- assign("x", data$speed, envir = .GlobalEnv)
  return(gnm)
}

# ---------------------------------------------------------------------------
# process_subject_vibro() -- fit GNM + bootstrap for one subject x condition
#
# Input:  sub_data  -- data frame for one subject in one vibration condition
#                      (columns: speed, faster, slower, subject, vibration)
# Output: one-row tibble with PSE, slope, gamma, lambda, their bootstrap SEs,
#         and 95% bootstrap CIs (2.5th / 97.5th percentiles)
# ---------------------------------------------------------------------------
process_subject_vibro <- function(sub_data) {
  y <- with(sub_data, cbind(faster, slower))
  x <- assign("x", sub_data$speed, envir = .GlobalEnv)

  pmu <- pstart(x, y)

  # Point estimates from GNM
  fit_gnm <- gnlr(y = y,
                  distribution = "binomial",
                  mu = mu, pmu = pmu, iterlim = 10000)

  # 500 bootstrap replicates for SE and 95% CI
  boot_results <- boot::boot(
    data      = sub_data,
    statistic = PsychBootGNM_vibro,
    R         = 500,
    p         = 0.75
  )

  estimates <- PsychParametersGNM(fit_gnm, p = 0.75)
  ses  <- apply(boot_results$t, 2, sd, na.rm = TRUE)
  low  <- apply(boot_results$t, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  high <- apply(boot_results$t, 2, function(x) quantile(x, 0.975, na.rm = TRUE))

  param_tibble <- as_tibble(t(estimates))

  se_tibble <- as_tibble(t(ses))
  colnames(se_tibble) <- c("pse_se", "slope_se", "gamma_se", "lambda_se")

  low_tibble <- as_tibble(t(low))
  colnames(low_tibble) <- c("pse_low", "slope_low", "gamma_low", "lambda_low")

  high_tibble <- as_tibble(t(high))
  colnames(high_tibble) <- c("pse_high", "slope_high", "gamma_high", "lambda_high")

  bind_cols(
    param_tibble, se_tibble, low_tibble, high_tibble,
    tibble(subject   = unique(sub_data$subject),
           vibration = unique(sub_data$vibration))
  )
}
