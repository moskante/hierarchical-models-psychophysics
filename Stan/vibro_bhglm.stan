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
  vector[2] aa;             // mean of alpha by vibration
  vector[2] bb;             // mean of beta by vibration
  real<lower=0> tauaa;      // SD for aa
  real<lower=0> taubb;      // SD for bb
  real<lower=0> taua;       // subject-level SDs for alpha
  real<lower=0> taub;       // subject-level SDs for beta
  
  // Subject-level parameters
  matrix[nsubj, 2] alpha;
  matrix[nsubj, 2] beta;
}

transformed parameters {
  real mu[nobs]; 
  real<lower=0, upper=1> pi1[nobs];
  matrix[nsubj, 2] pse;
  matrix[nsubj, 2] jnd;
  vector[2] PSE;
  vector[2] JND;
  for (h in 1:2) {
    PSE[h] = -aa[h]/bb[h]; 
    JND[h]  = 0.6745/bb[h];        
  }
  real diffPSE = PSE[2] - PSE[1];
  real diffJND = JND[2] - JND[1];
  real diffSlope = bb[2] - bb[1];
  
  for (i in 1:nsubj) {
    for (h in 1:2) {
      pse[i, h] = -alpha[i, h]/beta[i, h]; 
      jnd[i, h]  = 0.6745/beta[i, h];        
    }
  }
  
  for (i in 1:nobs){
    int s = subject[i];
    int v = vibration[i] + 1;
    mu[i] = alpha[s,v]+beta[s,v]*x[i];
    pi1[i]= Phi(mu[i]);
  }
}

model {
  for (i in 1:nsubj){
    for (h in 1:2){
      alpha[i,h] ~ normal(aa[h], taua);
      beta[i,h] ~ normal(bb[h], taub);}
  }
  taua~cauchy(0,2.5);
  taub~cauchy(0,2.5);
  
  tauaa~cauchy(0,2.5);
  taubb~cauchy(0,2.5);
  
  aa ~ normal(0, tauaa);
  bb ~ normal(0, taubb);
  
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