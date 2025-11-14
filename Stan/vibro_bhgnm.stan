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
  // Hyperparameters
  vector[2] beta_0;         // mean of b_0 by vibration
  vector[2] beta_1;         // mean of b_1 by vibration
  real<lower=0> tau_beta0;  // SD for beta_0
  real<lower=0> tau_beta1;  // SD for beta_1
  real<lower=0> tau_b0;     // subject-level SD for b_0
  real<lower=0> tau_b1;     // subject-level SD for b_1
 
  // Subject-level parameters
  matrix[nsubj, 2] b_0;
  matrix[nsubj, 2] b_1;
  matrix<lower=0, upper=1>[nsubj, 2] gamma;
  matrix<lower=0, upper=1>[nsubj, 2] lambda;
}

transformed parameters {
  real<lower=0, upper=1> pi1[nobs];  
  real<lower=0, upper=1> pi2[nobs];  
  real mu[nobs]; 
  matrix[nsubj, 2] pse;
  matrix[nsubj, 2] jnd;
  vector[2] PSE;
  vector[2] JND;
  for (h in 1:2) {
    PSE[h] = beta_0[h]/beta_1[h]; 
    JND[h]  = 0.6745/beta_1[h];        
  }
  real diffPSE = PSE[2] - PSE[1];
  real diffJND = JND[2] - JND[1];
  real diffSlope = beta_1[2] - beta_1[1];
  
  for (i in 1:nsubj) {
    for (h in 1:2) {
      pse[i, h] = -b_0[i, h]/b_1[i, h]; 
      jnd[i, h]  = 0.6745/b_1[i, h];        
    }
  }
  
  for (i in 1:nobs){
    int s = subject[i];
    int v = vibration[i] + 1;
    mu[i] = b_0[s,v]+b_1[s,v]*x[i];
    pi2[i]=(1 - gamma[s,v] - lambda[s,v]) * Phi(mu[i]);
    pi1[i]=gamma[s,v]+pi2[i];
  }
}

model {
  for (i in 1:nsubj){
    for (h in 1:2){
      b_0[i,h] ~ normal(beta_0[h], tau_b0);
      b_1[i,h] ~ normal(beta_1[h], tau_b1);
      gamma[i,h] ~ uniform(0, 1);
      lambda[i,h]~ uniform(0, 1-gamma[i,h]);
    }
  }
  
  tau_beta0 ~ cauchy(0,2.5);
  tau_beta1 ~ cauchy(0,2.5);
  tau_b0 ~ cauchy(0,2.5);
  tau_b1 ~ cauchy(0,2.5);
  
  beta_0 ~ normal(0, tau_beta0);
  beta_1 ~ normal(0, tau_beta1);

  
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
