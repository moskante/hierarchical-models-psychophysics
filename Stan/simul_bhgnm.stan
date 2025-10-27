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
  real beta[nsubj];
  real PSE;
  real bb; 
  real<lower=0, upper=1> gamma[nsubj];
  real<lower=0, upper=1> lambda[nsubj];
  real<lower=0> tauaa;
  real<lower=0> taubb;  
  real<lower=0> taua;
  real<lower=0> taub;
}

transformed parameters{
  real<lower=0, upper=1> pi1[nobs];  
  real<lower=0, upper=1> pi2[nobs]; 
  real mu[nobs];  
  real jnd[nsubj];
  real JND = 0.6745*bb;
  
  for (i in 1:nsubj){
    jnd[i] = 0.6745*beta[i];
  }
  
  for (i in 1:nobs){
    int s = subject[i];
    mu[i] = -pse[s]/beta[s]+x[i]/beta[s];
    //    pi2[i]=(1-gamma[subject[i]]-lambda[subject[i]])*Phi(mu[i]);
    pi2[i]=(1 - fmin(gamma[s] + lambda[s], 0.999))* Phi(mu[i]);
    pi1[i]=gamma[s]+pi2[i];
  }
}

model {
  
  for (i in 1:nsubj){
    pse[i] ~ normal(PSE, taua);
    beta[i] ~ normal(bb, taub);
    
    gamma[i] ~ uniform(0, 1);
    lambda[i]~ uniform(0, 1-gamma[i]);
  }
  
  
  for (i in 1:nobs){
    y[i] ~ binomial(n[i],pi1[i]);
  }
  
  taua~cauchy(0,2.5);
  taub~cauchy(0,2.5);
  tauaa~cauchy(0,2.5);
  taubb~cauchy(0,2.5);

  PSE ~ normal(0, tauaa);
  bb ~ normal(0, taubb);
  
}

generated quantities {
  vector[nobs] predProb;  // store predicted probability (mean)
  int y_sim[nobs];      // posterior predictive simulations

  for (i in 1:nobs) {
    predProb[i] = pi1[i]; // mean prediction for DHARMa
    y_sim[i] = binomial_rng(n[i], pi1[i]); // simulated observation
  }
}








