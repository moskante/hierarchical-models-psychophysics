// changed to compute the adjusted values of pse and jnd

data {
  int<lower=0> nobs;     //number of observation
  int<lower=0> nsubj;    //number of sujects
  int<lower=0> n[nobs];  //number of trials for each condition (x) and subject
  int<lower=1, upper=nsubj> subject[nobs];  //subject id 
  int<lower=0> y[nobs];  // number of successes for each condition (x) and subject
  real<lower=0> x[nobs]; // condition (x)
}

parameters {
  real mu[nsubj];       // Latent sensory mean (unadjusted location)
  real sigma[nsubj];    // Sensory scale (standard deviation)
  real MU;              // Latent population sensory mean
  real SIGMA;           // Population sensory scale
  real<lower=0> tau_mu;
  real<lower=0> tau_sigma;
  real<lower=0> tau_MU;
  real<lower=0> tau_SIGMA;  
  real<lower=0, upper=1> gamma[nsubj];
  real<lower=0, upper=1> lambda[nsubj];
}

transformed parameters{
  real<lower=0, upper=1> PI[nobs];  
  
  // Adjusted behavioral parameters (the ones we actually care about)
  real jnd[nsubj];      // Subject-level adjusted JND
  real JND;             // Population-level adjusted JND
  real pse[nsubj];      // Subject-level adjusted PSE
  real PSE;             // Population-level adjusted PSE
  
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
  SIGMA ~ normal(0, tau_SIGMA);
  
  // Hyperpriors
  tau_mu ~ cauchy(0,2.5);
  tau_sigma ~ cauchy(0,2.5);
  tau_MU ~ cauchy(0,2.5);
  tau_SIGMA ~ cauchy(0,2.5);
  
  // Hierarchical structure
  for (i in 1:nsubj){
    mu[i] ~ normal(MU, tau_mu);
    sigma[i] ~ normal(SIGMA, tau_sigma);
    gamma[i] ~ uniform(0, 1);
    lambda[i] ~ uniform(0, 1-gamma[i]);
  }
  
  // Likelihood
  for (i in 1:nobs){
    y[i] ~ binomial(n[i], PI[i]);
  }
}

generated quantities {
  vector[nobs] predProb;  // store predicted probability (mean)
  int y_sim[nobs];      // posterior predictive simulations

  for (i in 1:nobs) {
    predProb[i] = PI[i]; // mean prediction for DHARMa
    y_sim[i] = binomial_rng(n[i], PI[i]); // simulated observation
  }
}