data {
  int<lower=0> nobs;          
  int<lower=0> nsubj;         
  int<lower=0> n[nobs];       
  int subject[nobs];          
  array[nobs] int<lower=0> y; 
  real x[nobs];               
}

parameters {
  real pse[nsubj];
  real<lower=0> sigma[nsubj];
  real PSE;
  real log_SIGMA; // Latent parameter on log scale             
  real<lower=0> tau_pse;
  real<lower=0> tau_sigma;    // Standard deviation on log scale
  real<lower=0> tau_PSE;
  real<lower=0> tau_SIGMA;    // Prior scale for log_SIGMA
}

transformed parameters {
  real<lower=0, upper=1> PI[nobs];
  real jnd[nsubj];
  
  // 1. Calculate true population mean on original scale (E[sigma_i])
  real<lower=0> SIGMA = exp(log_SIGMA + 0.5 * square(tau_sigma));
  
  // 2. Derive JND using corrected SIGMA
  real JND = 0.6745 * SIGMA;
  
  for (i in 1:nsubj) {
    jnd[i] = 0.6745 * sigma[i];
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
    sigma[i] ~ lognormal(log_SIGMA, tau_sigma); // Smooth sampling
  }
  
  // Likelihood
  for (i in 1:nobs) {
    y[i] ~ binomial(n[i], PI[i]);
  }
}

generated quantities {
  vector[nobs] predProb;  
  int y_sim[nobs];        
  
  for (i in 1:nobs) {
    predProb[i] = PI[i]; 
    y_sim[i] = binomial_rng(n[i], PI[i]); 
  }
}
