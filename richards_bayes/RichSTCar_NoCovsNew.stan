
functions {
  real log_RichLin(real lr, real lh, real ci, real ls, real ti) {
    real h = exp(lh);
    real s = exp(ls);
    real lout;
    
    real lhct = log1p_exp(h * (ci-ti)-ls);
    if (is_inf(lhct))
    {
      lhct = h * (ci-ti) - ls;
    }
    if (exp(h * (ci-ti) - ls) < 10^(-14))
    {
      lhct = exp(h * (ci-ti) - ls);
    }
    
    lout = lr + lh + h*(ci-ti) - (s+1) * lhct;
    
    return lout;
  }
  
  real sparse_carar_lpdf(vector phi, real alpha, 
  int[,] W_sparse, vector W_weight, vector D_sparse, vector lambda, int n, int W_n) {
    row_vector[n] phit_D; // phi' * D
    row_vector[n] phit_W; // phi' * W
    vector[n] ldet_terms;
    
    phit_D = (phi .* D_sparse)';
    phit_W = rep_row_vector(0, n);
    for (i in 1:W_n) {
      phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + W_weight[i]*phi[W_sparse[i, 2]];
      phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + W_weight[i]*phi[W_sparse[i, 1]];
    }
    
    for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
    return 0.5 * (sum(ldet_terms)
    - (phit_D * phi - alpha * (phit_W * phi)));
  }
  
}


data {
  // Data block
  int<lower=0> N; // Number of obs
  int Y[N]; // Vector of Obs
  
  int<lower=0> Nreg;  // Number of regions
  int<lower=0> Ntimes;  // Number of times
  real<lower=0> t[Ntimes]; // Times
  
  matrix<lower = 0>[Nreg, Nreg] W; // Adjacency
  int W_n;                         // Number of adjacent region pairs
  
  vector[N] lOff;   //Offset
  
  int rId[N]; // Associate obs to region
  int tId[N]; // Associate obs to time
  
}

transformed data {
  // Vector of zeros
  vector[Nreg] zeros;
  // Data for sparse car
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[W_n] W_weight;     // Connection weights
  vector[Nreg] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[Nreg] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  
  zeros = rep_vector(0, Nreg);
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
  for (i in 1:(Nreg - 1)) {
    for (j in (i + 1):Nreg) {
      if (W[i, j] > 0) {
        W_sparse[counter, 1] = i;
        W_sparse[counter, 2] = j;
        W_weight[counter] = W[i, j];
        counter = counter + 1;
      }
    }
  }
  }
  
  for (i in 1:Nreg) D_sparse[i] = sum(W[i]);
  
  {
    vector[Nreg] invsqrtD;  
    for (i in 1:Nreg) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}

parameters {  
  // Parameters block
  
  // Baseline parameter
  real logbase;
  // Richard's parameter
  real logr;
  real logh;
  real c;
  real logs;
  
  // CAR parameters
  vector[Nreg] phis[Ntimes];
  real<lower=0> sigma;
  real<lower = 0, upper = 1> alpha;
  
  // AR parameters
  real<lower=0, upper=1> rho;
  
}

transformed parameters {
  // Transformed parameters block
  vector[N] lm;
  vector[Nreg] phi[Ntimes];
  
  phi[1] = sigma*phis[1];
  for (i in 2:Ntimes){
    phi[i] = rho*phis[i] + sigma*phis[i];
  }
  
  // Mean function
  for (i in 1:N) {
    lm[i] = phi[tId[i]][rId[i]] + 
    log_sum_exp(logbase, log_RichLin(logr, logh, c, logs, t[tId[i]]));
  }
  lm = lOff + lm;
} 


model {
  // Model block
  // Baseline priors
  logbase ~ normal(0, 1);
  // Richards Priors
  logr ~ normal(0, 1);
  logh ~ normal(0, 1);
  c ~ normal(t[Ntimes], t[Ntimes]);
  logs ~ normal(0, 1);
  
  // CAR prior
  phis[1] ~ sparse_carar(alpha, W_sparse, W_weight, D_sparse, lambda, Nreg, W_n);
  for (i in 2:Ntimes){
    phis[i] ~ sparse_carar(alpha, W_sparse, W_weight, D_sparse, lambda, Nreg, W_n);
  }
  alpha ~ beta(.5, .5);
  sigma ~ normal(0, 1);
  
  // likelihood
  Y ~ poisson_log(lm); 
}

generated quantities{
  // Output block
  real Y_pred[N];
  vector[N] log_lik;  
  
  Y_pred = poisson_log_rng(lm);  // Posterior predictive distribution
  for (i in 1:N){
    log_lik[i] = poisson_log_lpmf(Y[i] | lm[i]);
  }

}

