data {
  int<lower=0> nobs;     // number of observations
  int<lower=0> nsubj;    // number of subjects
  int<lower=0> n[nobs];  // number of trials for each condition (x) and subject
  int<lower=1, upper=nsubj> subject[nobs];  // subject id 
  int<lower=0> y[nobs];  // number of successes for each condition (x) and subject
  real<lower=0> x[nobs]; // condition (x)
}

parameters {
  // 1. Latent Sensory Parameters (Sensory space)
  real mu[nsubj];       // Individual latent sensory mean (unadjusted location)
  real sigma[nsubj];    // Individual sensory scale (standard deviation)
  real MU;              // Population latent sensory mean
  real SIGMA;           // Population sensory scale
  real<lower=0> tau_mu;
  real<lower=0> tau_sigma;
  real<lower=0> tau_MU;
  real<lower=0> tau_SIGMA;  

  // 2. Asymptotic Parameters (Probability space [0, 1])
  real<lower=0, upper=1> gamma[nsubj];       // Individual guess rate
  real<lower=0, upper=1> lambda_raw[nsubj];  // Individual remaining lapse fraction
  
  real<lower=0, upper=1> GAMMA;              // Population mean guess rate
  real<lower=0, upper=1> LAMBDA_tilde;       // Population mean lapse fraction
  
  real<lower=0> tau_gamma;                   // Between-subjects guess rate variability
  real<lower=0> tau_lambda_raw;              // Between-subjects lapse fraction variability
}

transformed parameters {
  real<lower=0, upper=1> PI[nobs];  
  
  // Realized behavioral lapse rates matching the geometric bounds
  real<lower=0, upper=1> lambda[nsubj];
  real<lower=0, upper=1> LAMBDA;

  // Final corrected psychophysical metrics (directly exported to R)
  real jnd[nsubj];      // Subject-level adjusted JND
  real JND;             // Population-level adjusted JND
  real pse[nsubj];      // Subject-level adjusted PSE
  real PSE;             // Population-level adjusted PSE

  // Algebraic coupling: mathematically ensures gamma + lambda < 1
  for (i in 1:nsubj) {
    lambda[i] = (1.0 - gamma[i]) * lambda_raw[i];
  }
  
  // Population-level realized lapse rate
  LAMBDA = (1.0 - GAMMA) * LAMBDA_tilde;

  // 1. Subject-level Parameter Adjustments
  for (i in 1:nsubj) {
    // Numeric safety guardrails for the inverse cumulative normal
    real val_pse = (0.5 - gamma[i]) / (1.0 - gamma[i] - lambda[i]);
    real val_pse_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_pse));
    real k = inv_Phi(val_pse_bounded) * sigma[i];

    pse[i] = mu[i] + k;

    real val_jnd = (0.75 - gamma[i]) / (1.0 - gamma[i] - lambda[i]);
    real val_jnd_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_jnd));
    jnd[i] = inv_Phi(val_jnd_bounded) * sigma[i] - k;
  }

  // 2. Population-level Parameter Adjustments
  {
    real val_pse_pop = (0.5 - GAMMA) / (1.0 - GAMMA - LAMBDA);
    real val_pse_pop_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_pse_pop));
    real k_pop = inv_Phi(val_pse_pop_bounded) * SIGMA;

    PSE = MU + k_pop;

    real val_jnd_pop = (0.75 - GAMMA) / (1.0 - GAMMA - LAMBDA);
    real val_jnd_pop_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_jnd_pop));
    JND = inv_Phi(val_jnd_pop_bounded) * SIGMA - k_pop;
  }
  
  // Binomial response probability mapping
  for (i in 1:nobs) {
    int s = subject[i];
    PI[i] = gamma[s] + (1.0 - gamma[s] - lambda[s]) * Phi(-mu[s]/sigma[s] + x[i]/sigma[s]);
  }
}

model {
  // Global Sensory Priors
  MU ~ normal(0, tau_MU);
  SIGMA ~ normal(0, tau_SIGMA);
  
  // Global Asymptotic Priors (Truncated Normals on the probability scale)
  GAMMA ~ normal(0.05, 0.1) T[0, 1];
  LAMBDA_tilde ~ normal(0.05, 0.1) T[0, 1];
  
  // Scale Hyper-priors (Weakly informative half-Cauchy)
  tau_mu ~ cauchy(0, 2.5);
  tau_sigma ~ cauchy(0, 2.5);
  tau_MU ~ cauchy(0, 2.5);
  tau_SIGMA ~ cauchy(0, 2.5);
  
  tau_gamma ~ cauchy(0, 0.5);
  tau_lambda_raw ~ cauchy(0, 0.5);

  // Full Hierarchical Structure (with safe bounding via T[0, 1])
  for (i in 1:nsubj) {
    mu[i] ~ normal(MU, tau_mu);
    sigma[i] ~ normal(SIGMA, tau_sigma);
    
    gamma[i] ~ normal(GAMMA, tau_gamma) T[0, 1];
    lambda_raw[i] ~ normal(LAMBDA_tilde, tau_lambda_raw) T[0, 1];
  }
  
  // Likelihood
  for (i in 1:nobs) {
    y[i] ~ binomial(n[i], PI[i]);
  }
}

generated quantities {
  vector[nobs] predProb;  // store predicted probability (mean)
  int y_sim[nobs];        // posterior predictive simulations

  for (i in 1:nobs) {
    predProb[i] = PI[i]; 
    y_sim[i] = binomial_rng(n[i], PI[i]); 
  }
}
