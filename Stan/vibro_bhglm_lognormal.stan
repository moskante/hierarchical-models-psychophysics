data {
  int<lower=0> nobs;                    // number of observations
  int<lower=0> nsubj;                   // number of subjects
  array[nobs] int<lower=0> n;           // number of trials
  array[nobs] int subject;              // subject id
  array[nobs] int<lower=0> y;           // observed counts (successes)
  array[nobs] real x;                   // stimulus intensity
  array[nobs] int<lower=0,upper=1> vibration; // vibration condition (0 or 1)
}

parameters {
  matrix[nsubj, 2] b0;
  matrix<lower=0>[nsubj, 2] b1;  // Corretto: b1 dichiarato al posto di log_b1
  
  vector[2] beta0;     
  vector[2] log_beta1;           // Population-level log(beta1)
 
  real<lower=0> tau_b0;       
  real<lower=0> tau_b1;     
  real<lower=0> tau_beta0;       // Prior scale for beta0
  real<lower=0> tau_beta1;       // Prior scale for log_beta1
}

transformed parameters {
  array[nobs] real<lower=0, upper=1> PI;
  vector<lower=0>[2] beta1;
  
  // 1. Calculate true population mean on original scale (E[sigma_i])
  beta1 = exp(log_beta1 + 0.5 * square(tau_b1));//
  //beta1 = exp(log_beta1 );//Median
  
  // Calcolo probabilità per la verosimiglianza
  for (i in 1:nobs){
    int s = subject[i];
    int v = vibration[i] + 1;
    PI[i] = Phi(b0[s,v] + b1[s,v] * x[i]);
  }
}

model {
  // Priors
  beta0 ~ normal(0, tau_beta0);
  log_beta1 ~ normal(0, tau_beta1);
  
  // Hyperpriors
   // tau_b0 ~ cauchy(0, 2.5);
   // tau_b1 ~ cauchy(0, 2.5);
   // tau_beta0 ~ cauchy(0, 2.5);
   // tau_beta1 ~ cauchy(0, 2.5);
  tau_b0 ~ normal(0, 2.5);
  tau_b1 ~ normal(0, 1);       // Ridotto a 1 perché è sulla scala logaritmica
  tau_beta0 ~ normal(0, 2.5);
  tau_beta1 ~ normal(0, 1);    // Ridotto a 1 per evitare esplosioni nell'esponenziale
  
  // Hierarchical structure
  for (i in 1:nsubj){
    for (h in 1:2){
      b0[i,h] ~ normal(beta0[h], tau_b0);
      b1[i,h] ~ lognormal(log_beta1[h], tau_b1);
    }
  }
  
  // Likelihood
  for (i in 1:nobs){
    y[i] ~ binomial(n[i], PI[i]);
  }
}

generated quantities {
  matrix[nsubj, 2] pse;
  matrix[nsubj, 2] jnd;
  vector[2] PSE;
  vector[2] JND;
  real diffSlope;
  
  vector[nobs] predProb;  
  array[nobs] int y_sim;      

  // Valori derivati calcolati a posteriori (spostati per efficienza)
  for (i in 1:nsubj) {
    for (h in 1:2) {
      pse[i, h] = -b0[i, h] / b1[i, h]; 
      jnd[i, h] = 0.6745 / b1[i, h];        
    }
  }
  
  for (h in 1:2) {
    PSE[h] = -beta0[h] / beta1[h]; 
    JND[h] = 0.6745 / beta1[h];        
  }
  
  diffSlope = beta1[2] - beta1[1];
  
  // Posterior predictive simulations
  for (i in 1:nobs) {
    predProb[i] = PI[i]; 
    y_sim[i] = binomial_rng(n[i], PI[i]); 
  }
}