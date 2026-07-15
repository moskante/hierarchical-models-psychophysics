data {
  int<lower=0> nobs;                    // number of observations
  int<lower=0> nsubj;                   // number of subjects
  int<lower=0> n[nobs];                 // number of trials
  int subject[nobs];                    // subject id
  array[nobs] int<lower=0> y;           // observed counts (successes)
  real x[nobs];                         // stimulus intensity
  int<lower=0,upper=1> vibration[nobs]; // vibration condition (0 or 1)
  
}

parameters {
  matrix[nsubj, 2] b0;
  matrix[nsubj, 2] b1;
  vector[2] beta0;             
  vector[2] beta1;  
  real<lower=0> tau_b0;       
  real<lower=0> tau_b1;     
  real<lower=0> tau_beta0;      
  real<lower=0> tau_beta1;      
}

transformed parameters {
  real<lower=0, upper=1> PI[nobs];
  matrix[nsubj, 2] pse;
  matrix[nsubj, 2] jnd;
  vector[2] PSE;
  vector[2] JND;
  
  for (i in 1:nsubj) {
    for (h in 1:2) {
      pse[i, h] = -b0[i, h]/b1[i, h]; 
      jnd[i, h]  = 0.6745/b1[i, h];        
    }
  }
  
  for (h in 1:2) {
    PSE[h] = -beta0[h]/beta1[h]; 
    JND[h]  = 0.6745/beta1[h];        
  }
  
  real diffSlope = beta1[2] - beta1[1];
  
  for (i in 1:nobs){
    int s = subject[i];
    int v = vibration[i] + 1;
    PI[i]= Phi(b0[s,v]+b1[s,v]*x[i]);
  }
}

model {
  // Priors
  beta0 ~ normal(0, tau_beta0);
  beta1 ~ normal(0, tau_beta1);
  
  // Hyperpriors
  tau_b0 ~ cauchy(0,2.5);
  tau_b1 ~ cauchy(0,2.5);
  tau_beta0 ~ cauchy(0,2.5);
  tau_beta1 ~ cauchy(0,2.5);
  
  // Hierarchical structure
  for (i in 1:nsubj){
    for (h in 1:2){
      b0[i,h] ~ normal(beta0[h], tau_b0);
      b1[i,h] ~ normal(beta1[h], tau_b1);}
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