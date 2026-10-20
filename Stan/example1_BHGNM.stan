data {
  int<lower=0> nobs;     //number of observations
  int<lower=0> nsubj;    //number of subjects
  array[nobs] int<lower=0> n;  //number of trials for each condition (x) and subject
  array[nobs] int<lower=1, upper=nsubj> subject;  //subject id 
  array[nobs] int<lower=0> y;  // number of successes for each condition (x) and subject
  array[nobs] real<lower=0> x; // condition (x)
}

parameters {
  array[nsubj] real mu;       // Latent sensory mean (unadjusted location)
  array[nsubj] real<lower=0> sigma;    // Sensory scale (standard deviation)
  real MU;              // Latent population sensory mean
  real log_SIGMA;             // Latent population scale on log scale
  real<lower=0> tau_mu;
  real<lower=0> tau_sigma;
  real<lower=0> tau_MU;
  real<lower=0> tau_SIGMA;  
  array[nsubj] real<lower=0, upper=1> gamma;
  array[nsubj] real<lower=0, upper=1> lambda;
}

transformed parameters{
  array[nobs] real<lower=0, upper=1> PI;  
  
  // Adjusted behavioral parameters
  array[nsubj] real jnd;      // Subject-level adjusted JND
  real JND;             // Population-level adjusted JND
  array[nsubj] real pse;      // Subject-level adjusted PSE
  real PSE;             // Population-level adjusted PSE
  
  // True population arithmetic mean scale E[sigma_i]
  real<lower=0> SIGMA = exp(log_SIGMA + 0.5 * square(tau_sigma));
  
  // Empirical population asymptotes
  real GAMMA;           
  real LAMBDA;
  
  GAMMA = mean(gamma);
  LAMBDA = mean(lambda);

  // 1. Subject-level Adjusted parameters
  for (i in 1:nsubj){
    // Calculate Subject-level Shift k_i
    real val_pse = (0.5 - gamma[i]) / (1.0 - gamma[i] - lambda[i]);
    real val_pse_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_pse));
    real k = inv_Phi(val_pse_bounded) * sigma[i];

    pse[i] = mu[i] + k;

    // Calculate Subject-level Adjusted JND(0.75)
    real val_jnd = (0.75 - gamma[i]) / (1.0 - gamma[i] - lambda[i]);
    real val_jnd_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_jnd));
    jnd[i] = inv_Phi(val_jnd_bounded) * sigma[i] - k;
  }

  // 2. Population-level Adjusted parameters
  {
    // Calculate Population-level Shift k_pop
    real val_pse_pop = (0.5 - GAMMA) / (1.0 - GAMMA - LAMBDA);
    real val_pse_pop_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_pse_pop));
    real k_pop = inv_Phi(val_pse_pop_bounded) * SIGMA;

    PSE = MU + k_pop;

    // Calculate Population-level Adjusted JND(0.75)
    real val_jnd_pop = (0.75 - GAMMA) / (1.0 - GAMMA - LAMBDA);
    real val_jnd_pop_bounded = fmax(1e-5, fmin(1.0 - 1e-5, val_jnd_pop));
    JND = inv_Phi(val_jnd_pop_bounded) * SIGMA - k_pop;
  }

  // 3. Probabilities mapping
  for (i in 1:nobs){
    int s = subject[i];
    PI[i] = gamma[s] + (1 - gamma[s] - lambda[s]) * Phi(-mu[s]/sigma[s] + x[i]/sigma[s]);
  }
}

model {
  // Priors
  MU ~ normal(0, tau_MU);
  log_SIGMA ~ normal(0, tau_SIGMA);
  
  // Hyperpriors
  tau_mu ~ cauchy(0,2.5);
  tau_sigma ~ cauchy(0,2.5);
  tau_MU ~ cauchy(0,2.5);
  tau_SIGMA ~ cauchy(0,2.5);
  
  // Hierarchical structure
  for (i in 1:nsubj){
    mu[i] ~ normal(MU, tau_mu);
    sigma[i] ~ lognormal(log_SIGMA, tau_sigma); // Log-normal distribution for scale
    gamma[i] ~ uniform(0, 1);
    lambda[i] ~ uniform(0, 1 - gamma[i]);
  }
  
  // Likelihood
  for (i in 1:nobs){
    y[i] ~ binomial(n[i], PI[i]);
  }
}

generated quantities {
  vector[nobs] predProb;  // store predicted probability (mean)
  array[nobs] int y_sim;      // posterior predictive simulations

  for (i in 1:nobs) {
    predProb[i] = PI[i]; // mean prediction for DHARMa
    y_sim[i] = binomial_rng(n[i], PI[i]); // simulated observation
  }
}
