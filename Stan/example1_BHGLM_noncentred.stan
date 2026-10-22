// Alternative BH-GLM for Example 1: non-centred parameterisation of sigma.
// Same model and priors as example1_BHGLM.stan, but the subject-level scale
// parameters are written as (Eq. 28 in the paper)
//   sigma[i] = exp(log_SIGMA + tau_sigma * z_scale[i]),  z_scale[i] ~ normal(0, 1)
// This implies the same prior as the centred form
//   sigma[i] ~ lognormal(log_SIGMA, tau_sigma)
// so the posterior is the same, but the sampler avoids the funnel between
// sigma[i] and tau_sigma (fewer divergent transitions, better convergence
// when tau_sigma is small). pse[i] keeps the centred form.

data {
  int<lower=0> nobs;
  int<lower=0> nsubj;
  array[nobs] int<lower=0> n;
  array[nobs] int subject;
  array[nobs] int<lower=0> y;
  array[nobs] real x;
}

parameters {
  array[nsubj] real pse;
  vector[nsubj] z_scale;      // standardised subject deviations, log(sigma)
  real PSE;
  real log_SIGMA; // Latent parameter on log scale
  real<lower=0> tau_pse;
  real<lower=0> tau_sigma;    // Standard deviation on log scale
  real<lower=0> tau_PSE;
  real<lower=0> tau_SIGMA;    // Prior scale for log_SIGMA
}

transformed parameters {
  array[nobs] real<lower=0, upper=1> PI;
  array[nsubj] real sigma;
  array[nsubj] real jnd;
  real<lower=0> SIGMA = exp(log_SIGMA + 0.5 * square(tau_sigma));  // population mean of sigma_i on the original scale (E[sigma_i])
  real JND = 0.6745 * SIGMA;  // population JND from the mean SIGMA

  // subject-level scale (non-centred)
  for (i in 1:nsubj) {
    sigma[i] = exp(log_SIGMA + tau_sigma * z_scale[i]);  // exp() keeps sigma positive
    jnd[i]   = 0.6745 * sigma[i];
  }
  for (i in 1:nobs) {
    int s = subject[i];
    PI[i] = Phi(-pse[s]/sigma[s] + x[i]/sigma[s]);
  }
}

model {
  // Priors
  PSE ~ normal(0, tau_PSE);
  log_SIGMA ~ normal(0, tau_SIGMA); // Prior on log scale

  // Hyperpriors
  tau_pse ~ cauchy(0, 2.5);
  tau_sigma ~ cauchy(0, 2.5);
  tau_PSE ~ cauchy(0, 2.5);
  tau_SIGMA ~ cauchy(0, 2.5);

  // Hierarchical structure
  for (i in 1:nsubj) {
    pse[i] ~ normal(PSE, tau_pse);
  }
  z_scale ~ std_normal();     // non-centred: implies sigma[i] ~ lognormal(log_SIGMA, tau_sigma)

  // Likelihood
  for (i in 1:nobs) {
    y[i] ~ binomial(n[i], PI[i]);
  }
}

generated quantities {
  vector[nobs] predProb;
  array[nobs] int y_sim;

  for (i in 1:nobs) {
    predProb[i] = PI[i];
    y_sim[i] = binomial_rng(n[i], PI[i]);
  }
}
