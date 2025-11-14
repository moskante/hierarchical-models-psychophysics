data {
  int<lower=0> nobs;     //number of observation
  int<lower=0> nsubj;    //number of sujects
  int<lower=0> n[nobs];  //number of trials for each condition (x) and subject
  int<lower=1, upper=nsubj> subject[nobs];  // id 
  int<lower=0> y[nobs];  // number of successes for each condition (x) and subject
  real<lower=0> x[nobs]; // condition (x)
}

parameters {
  real pse[nsubj];
  real sigma[nsubj];
  real PSE;
  real SIGMA; 
  real<lower=0, upper=1> gamma[nsubj];
  real<lower=0, upper=1> lambda[nsubj];
  real<lower=0> tau_PSE;
  real<lower=0> tau_SIGMA;  
  real<lower=0> tau_pse;
  real<lower=0> tau_sigma;
}

transformed parameters{
  real<lower=0, upper=1> pi1[nobs];  
  real<lower=0, upper=1> pi2[nobs]; 
  real mu[nobs];  
  real jnd[nsubj];
  real JND = 0.6745*SIGMA;
  
  for (i in 1:nsubj){
    jnd[i] = 0.6745*sigma[i];
  }
  
  for (i in 1:nobs){
    int s = subject[i];
    mu[i] = -pse[s]/sigma[s]+x[i]/sigma[s];
    pi2[i]=(1 - gamma[s] - lambda[s])* Phi(mu[i]);
    pi1[i]=gamma[s]+pi2[i];
  }
}

model {
  
  for (i in 1:nsubj){
    pse[i] ~ normal(PSE, tau_pse);
    sigma[i] ~ normal(SIGMA, tau_sigma);
    
    gamma[i] ~ uniform(0, 1);
    lambda[i]~ uniform(0, 1-gamma[i]);
  }
  
  
  for (i in 1:nobs){
    y[i] ~ binomial(n[i],pi1[i]);
  }
  
  tau_pse ~ cauchy(0,2.5);
  tau_sigma ~ cauchy(0,2.5);
  tau_PSE ~ cauchy(0,2.5);
  tau_SIGMA ~ cauchy(0,2.5);

  PSE ~ normal(0, tau_PSE);
  SIGMA ~ normal(0, tau_SIGMA);
  
}

generated quantities {
  vector[nobs] predProb;  // store predicted probability (mean)
  int y_sim[nobs];      // posterior predictive simulations

  for (i in 1:nobs) {
    predProb[i] = pi1[i]; // mean prediction for DHARMa
    y_sim[i] = binomial_rng(n[i], pi1[i]); // simulated observation
  }
}








