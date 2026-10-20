// Alternative BH-GLM for Example 2: lognormal hierarchical prior on the slopes.
// Non-centred parameterisation: subject-level effects are written as
//   b0[i,h] = beta0[h] + tau_b0 * z_b0[i,h]
//   b1[i,h] = exp(log_beta1[h] + tau_b1 * z_b1[i,h])
// with z ~ normal(0, 1). This implies the same priors as the centred form
//   b0[i,h] ~ normal(beta0[h], tau_b0),  b1[i,h] ~ lognormal(log_beta1[h], tau_b1)
// but samples efficiently when the between-subject SDs are small.

data {
  int<lower=0> nobs;                    // number of observations
  int<lower=0> nsubj;                   // number of subjects
  array[nobs] int<lower=0> n;           // number of trials
  array[nobs] int subject;              // subject id
  array[nobs] int<lower=0> y;           // observed counts (successes)
  array[nobs] real x;                   // stimulus intensity
  array[nobs] int<lower=0,upper=1> vibration; // vibration condition (0 or 1)
}

parameters {
  matrix[nsubj, 2] z_b0;         // standardised subject deviations, intercept
  matrix[nsubj, 2] z_b1;         // standardised subject deviations, log(slope)

  vector[2] beta0;
  vector[2] log_beta1;           // Population-level log(beta1)

  real<lower=0> tau_b0;
  real<lower=0> tau_b1;
  real<lower=0> tau_beta0;       // Prior scale for beta0
  real<lower=0> tau_beta1;       // Prior scale for log_beta1
}

transformed parameters {
  matrix[nsubj, 2] b0;
  matrix<lower=0>[nsubj, 2] b1;
  array[nobs] real<lower=0, upper=1> PI;
  vector<lower=0>[2] beta1;

  // subject-level intercepts and slopes (non-centred)
  for (h in 1:2) {
    b0[, h] = beta0[h] + tau_b0 * z_b0[, h];
    b1[, h] = exp(log_beta1[h] + tau_b1 * z_b1[, h]);
  }

  // population mean of b1 on the original scale (E[b1])
  beta1 = exp(log_beta1 + 0.5 * square(tau_b1));
  // alternative: population median, beta1 = exp(log_beta1)

  // response probabilities for the likelihood
  for (i in 1:nobs){
    int s = subject[i];
    int v = vibration[i] + 1;
    PI[i] = Phi(b0[s,v] + b1[s,v] * x[i]);
  }
}

model {
  // Priors
  beta0 ~ normal(0, tau_beta0);
  log_beta1 ~ normal(0, tau_beta1);

  // Hyperpriors (half-normal, since the scales are constrained to be positive)
  tau_b0 ~ normal(0, 2.5);
  tau_b1 ~ normal(0, 1);       // scale 1 because b1 is on the log scale
  tau_beta0 ~ normal(0, 2.5);
  tau_beta1 ~ normal(0, 1);    // scale 1 to keep exp(log_beta1) stable

  // Hierarchical structure (non-centred)
  to_vector(z_b0) ~ std_normal();
  to_vector(z_b1) ~ std_normal();

  // Likelihood
  for (i in 1:nobs){
    y[i] ~ binomial(n[i], PI[i]);
  }
}

generated quantities {
  matrix[nsubj, 2] pse;
  matrix[nsubj, 2] jnd;
  vector[2] PSE;
  vector[2] JND;
  real diffSlope;

  vector[nobs] predProb;
  array[nobs] int y_sim;

  // derived quantities (computed here for efficiency)
  for (i in 1:nsubj) {
    for (h in 1:2) {
      pse[i, h] = -b0[i, h] / b1[i, h];
      jnd[i, h] = 0.6745 / b1[i, h];
    }
  }

  for (h in 1:2) {
    PSE[h] = -beta0[h] / beta1[h];
    JND[h] = 0.6745 / beta1[h];
  }

  diffSlope = beta1[2] - beta1[1];

  // Posterior predictive simulations
  for (i in 1:nobs) {
    predProb[i] = PI[i];
    y_sim[i] = binomial_rng(n[i], PI[i]);
  }
}
