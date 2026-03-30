# gnlm_functions_psejnd.R
#
# GNM (Generalised Nonlinear Model) helper functions for fitting psychometric
# functions and extracting PSE and JND.
#
# The psychometric function used here is the cumulative normal (probit) with
# two extra free parameters:
#   gamma  -- lower asymptote (guessing rate): the probability of responding
#             "longer" even when the stimulus intensity is far below threshold.
#   lambda -- lapse rate: the probability of an incorrect response at
#             intensities well above threshold (upper asymptote = 1 - lambda).
#
# The full function is:
#   P(x) = gamma + (1 - gamma - lambda) * Phi((x - mu) / sigma)
#
# where Phi is the standard normal CDF, mu is the PSE, and sigma is the
# spread (inversely related to the slope / JND).
#
# NOTE: the variable `x` (stimulus levels for the current subject) must exist
# in the calling environment before gnlr() is called, because gnlr() evaluates
# mu() in that environment. This is handled by the assign() calls in
# process_subject() and PsychBootGNM().

# ---------------------------------------------------------------------------
# mu() -- psychometric function evaluated at the current x values
#
# Parameters p are on an unconstrained scale so that gnlr() can optimise
# freely. gamma and lambda are recovered via the atan2 transformation:
#   gamma  = atan2(p[3], 1) / pi + 0.5   (maps real line -> (0, 1))
#   lambda = (1 - gamma) * (atan2(p[4], 1) / pi + 0.5)
# ---------------------------------------------------------------------------
mu <- function(p){
  mu_i    <- p[1]   # PSE (location parameter, same units as x)
  sigma_i <- p[2]   # spread (sigma of the underlying normal; JND = qnorm(0.75)*sigma)
  # Back-transform guessing and lapse parameters from unconstrained scale
  gamma_i  <- atan2(p[3], 1) / pi + 0.5
  lambda_i <- (1 - gamma_i) * (atan2(p[4], 1) / pi + 0.5)
  gamma_i + (1 - gamma_i - lambda_i) * pnorm(x, mean = mu_i, sd = sigma_i)
}

# ---------------------------------------------------------------------------
# pstart() -- generate sensible starting values for gnlr() optimisation
#
# Fits a probit GLM first, then converts its coefficients to the
# (mu, sigma, p3, p4) parameterisation used by mu() above.
# ---------------------------------------------------------------------------
pstart <- function(x, y){
  glm_start  <- glm(y ~ x, family = binomial(link = "probit"))
  coef_start <- coef(glm_start)
  # Guard against degenerate GLM fits where the slope is zero or missing
  if (is.na(coef_start[2]) || coef_start[2] == 0) {
    coef_start[2] <- 1e-6
  }
  prop_correct    <- y[, 1] / rowSums(y)
  # Guessing rate start: minimum observed proportion (floored at 0.99)
  gamma_i_start   <- min(min(prop_correct, na.rm = TRUE), 0.99)
  # Lapse rate start: distance of maximum observed proportion from 1
  lambda_i_start  <- max(1 - max(prop_correct, na.rm = TRUE), 0.01)
  mu_i_start      <- -coef_start[1] / coef_start[2]   # PSE from GLM intercept/slope
  sigma_i_start   <- 1 / coef_start[2]                # spread from GLM slope
  # Transform gamma and lambda to unconstrained p3/p4 starting values
  p3_start <- tan(pi * (gamma_i_start - 0.5))
  p4_start <- tan(pi * ((lambda_i_start / (1 - gamma_i_start)) - 0.5))
  return(c(mu_i_start, sigma_i_start, p3_start, p4_start))
}

# ---------------------------------------------------------------------------
# PsychParametersGNM() -- extract PSE, JND, gamma, lambda from a fitted gnlr object
#
# p (default 0.75) sets the quantile used for JND:
#   JND = qnorm(p) * sigma
# At p = 0.75, JND equals the 75th-percentile spread of the underlying normal.
# ---------------------------------------------------------------------------
PsychParametersGNM <- function(gnm, p = 0.75){
  pout    <- gnm$coefficients
  pse     <- pout[1]                                  # location = PSE
  jnd     <- qnorm(p) * pout[2]                       # spread -> JND
  gamma_i  <- atan2(pout[3], 1) / pi + 0.5            # back-transform guessing rate
  lambda_i <- (1 - gamma_i) * (atan2(pout[4], 1) / pi + 0.5)  # back-transform lapse rate
  return(c(pse = pse, jnd = jnd, gamma = gamma_i, lambda = lambda_i))
}

# ---------------------------------------------------------------------------
# PsychBootGNM() -- bootstrap statistic function (passed to boot::boot)
#
# Resamples rows of `data`, re-fits the GNM, and returns the four parameters.
# Returns NA for all parameters if the fit fails (e.g. convergence failure).
# After each call, x is reset to the full subject dataset so that subsequent
# calls see the correct stimulus values.
# ---------------------------------------------------------------------------
PsychBootGNM <- function(data, indices, p = 0.75){
  data_resampled <- data[indices, ]
  y   <- with(data_resampled, cbind(Longer, Total - Longer))
  # Expose resampled x to the enclosing environment (required by gnlr/mu)
  x   <- assign("x", data_resampled$X, envir = .GlobalEnv)
  pmu <- pstart(x, y)
  gnm <- tryCatch(
    {
      model <- gnlr(y = y, distribution = "binomial",
                    mu = mu, pmu = pmu, iterlim = 10000)
      PsychParametersGNM(model, p)
    },
    error = function(e) rep(NA, 4)
  )
  # Restore x to the full subject data before returning
  x <- assign("x", data$X, envir = .GlobalEnv)
  return(gnm)
}

# ---------------------------------------------------------------------------
# process_subject() -- fit GNM + bootstrap for a single subject's data
#
# Input:  sub_data  -- data frame for one subject (columns: X, Longer, Total,
#                      Subject)
# Output: one-row tibble with PSE, JND, gamma, lambda and their bootstrap SEs
# ---------------------------------------------------------------------------
process_subject <- function(sub_data) {
  y <- with(sub_data, cbind(Longer, Total - Longer))
  # Expose stimulus values to the global environment for mu() evaluation
  x <- assign("x", sub_data$X, envir = .GlobalEnv)

  pmu <- pstart(x, y)

  # Fit the GNM to point estimates
  fit_gnm <- gnlr(y = y,
                  distribution = "binomial",
                  mu = mu, pmu = pmu, iterlim = 10000)

  # Bootstrap 500 replicates to estimate parameter SEs
  boot_results <- boot::boot(
    data      = sub_data,
    statistic = PsychBootGNM,
    R         = 500,
    p         = 0.75   # quantile for JND; change to adjust definition
  )

  estimates <- PsychParametersGNM(fit_gnm, p = 0.75)
  ses       <- apply(boot_results$t, 2, sd, na.rm = TRUE)

  # Assemble into a single-row tibble with estimates and bootstrap SEs
  param_tibble <- as_tibble(t(estimates))
  se_tibble    <- as_tibble(t(ses))
  colnames(se_tibble) <- c("pse_se", "jnd_se", "gamma_se", "lambda_se")

  bind_cols(param_tibble, se_tibble, tibble(Subject = unique(sub_data$Subject)))
}
