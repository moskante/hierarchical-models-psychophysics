data {
  int<lower=0> nobs;          //number of observation
  int<lower=0> nsubj;         //number of sujects
  int<lower=0> n[nobs];       //number of trials for each condition (x) and subject
  int subject[nobs];          //subject id          
  array[nobs] int<lower=0> y; // number of successes for each condition (x) and subject
  real x[nobs];               // condition (x)
}

parameters {
  real pse[nsubj];
  real sigma[nsubj];
  real PSE;
  real SIGMA;
  real<lower=0> tau_pse;
  real<lower=0> tau_sigma;
  real<lower=0> tau_PSE;
  real<lower=0> tau_SIGMA;
}

transformed parameters{
  real<lower=0, upper=1> PI[nobs];
  real jnd[nsubj];
  real JND = 0.6745*SIGMA;
  
  for (i in 1:nsubj){
    jnd[i] = 0.6745*sigma[i];
  }
  
  for (i in 1:nobs){
    int s = subject[i];
    PI[i]= Phi(-pse[s]/sigma[s]+x[i]/sigma[s]);
  }
}

model {
  // Priors
  PSE ~ normal(0, tau_PSE);
  SIGMA ~ normal(0, tau_SIGMA);
  
  // Hyperpriors
  tau_pse ~ cauchy(0,2.5);
  tau_sigma ~ cauchy(0,2.5);
  tau_PSE ~ cauchy(0,2.5);
  tau_SIGMA ~ cauchy(0,2.5);
  
  // Hyerarchical structure  
  for (i in 1:nsubj){
    pse[i] ~ normal(PSE, tau_pse);
    sigma[i] ~ normal(SIGMA, tau_sigma);
  }
  
  // Likelihood
  for (i in 1:nobs){
    y[i] ~ binomial(n[i],PI[i]);
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


