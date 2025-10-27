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
  real bb;
  real pse[nsubj];
  real beta[nsubj];
  real<lower=0> taua;
  real<lower=0> taub;
  real<lower=0> tauaa;
  real<lower=0> taubb;
}

transformed parameters{
  real mu[nobs]; 
  real<lower=0, upper=1> pi1[nobs];
  real jnd[nsubj];
  real JND = 0.6745*bb;
  
  for (i in 1:nsubj){
    jnd[i] = 0.6745*beta[i];
  }
  for (i in 1:nobs){
    int s = subject[i];
    mu[i] = -pse[s]/beta[s]+x[i]/beta[s];
    pi1[i]= Phi(mu[i]);
  }
}

model {
  
  PSE ~ normal(0, tauaa);
  bb ~ normal(0, taubb);
  taua~cauchy(0,2.5);
  taub~cauchy(0,2.5);
  tauaa~cauchy(0,2.5);
  taubb~cauchy(0,2.5);
    
  for (i in 1:nsubj){
    pse[i] ~ normal(PSE, taua);
    beta[i] ~ normal(bb, taub);
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


