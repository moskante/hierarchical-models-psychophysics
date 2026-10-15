# k_correction.R
#
# Utility: convert proportion correct (p) to the underlying sensitivity
# index k, after correcting for guessing (gamma) and lapsing (lambda).
#
# Background:
#   When a psychometric function includes lower (gamma) and upper (1-lambda)
#   asymptotes, the raw proportion correct is a mixture of true discriminative
#   performance and floor/ceiling effects. To isolate the sensitivity component,
#   we invert the psychometric function:
#
#     P(x) = gamma + (1 - gamma - lambda) * Phi(k)
#
#   Solving for k:
#     k = Phi^{-1}( (p - gamma) / (1 - gamma - lambda) )
#
#   k is expressed in standard normal units (z-score equivalent), so k = 0
#   corresponds to chance after correction, and larger values indicate higher
#   sensitivity.

compute_k <- function(p, gamma, lambda) {
  k <- qnorm((p - gamma) / (1 - gamma - lambda))
  return(k)
}

# Example: proportion correct = 0.50, guessing rate = 0.027, lapse rate = 0.024
compute_k(0.5, 0.027, 0.024)
