data {
  int<lower=0> nobs;
  int<lower=0> nsubj;
  int<lower=0> n[nobs];
  int subject[nobs];          
  array[nobs] int<lower=0> y;
  real x[nobs];
}
parameters {
  real PSE;
  real SIGMA;
  real pse[nsubj];
  real sigma[nsubj];
  real<lower=0> tau_pse;
  real<lower=0> tau_sigma;
  real<lower=0> tau_PSE;
  real<lower=0> tau_SIGMA;
}

transformed parameters{
  real mu[nobs]; 
  real<lower=0, upper=1> pi1[nobs];
  real jnd[nsubj];
  real JND = 0.6745*SIGMA;
  
  for (i in 1:nsubj){
    jnd[i] = 0.6745*sigma[i];
  }
  for (i in 1:nobs){
    int s = subject[i];
    mu[i] = -pse[s]/sigma[s]+x[i]/sigma[s];
    pi1[i]= Phi(mu[i]);
  }
}

model {
  
  PSE ~ normal(0, tau_PSE);
  SIGMA ~ normal(0, tau_SIGMA);
  tau_pse ~ cauchy(0,2.5);
  tau_sigma ~ cauchy(0,2.5);
  tau_PSE ~ cauchy(0,2.5);
  tau_SIGMA ~ cauchy(0,2.5);
    
  for (i in 1:nsubj){
    pse[i] ~ normal(PSE, tau_pse);
    sigma[i] ~ normal(SIGMA, tau_sigma);
  }
  
  for (i in 1:nobs){
    y[i] ~ binomial(n[i],pi1[i]);
  }
}

generated quantities {
  vector[nobs] predProb;  // store predicted probability (mean)
  int y_sim[nobs];      // posterior predictive simulations
  
  for (i in 1:nobs) {
    predProb[i] = pi1[i]; // mean prediction for DHARMa
    y_sim[i] = binomial_rng(n[i], pi1[i]); // simulated observation
  }
}


