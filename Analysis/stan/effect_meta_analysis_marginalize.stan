//
// Performs correlation meta-analysis


data {
  int<lower=0> N;    // number of observations
  int<lower=0> J;    // number of covariates
  int<lower=0> K;    // number of studies
  int<lower=0> N_full;      // number of observations for which all covariates are available
  int<lower=0> N_missing;   // number of observations with missing covariates
  int<lower=0> M;    // number of unique covariates that can be missing
  matrix[N, J] X;
  vector[N] y;
  
  array[N] int<lower=1, upper=K> map_to_study;    // map from observation to study
  // vector[J]<lower=1> n_covar_levels;    // number of levels per covariate
  array[N_full] int<lower=1, upper=N> ind_full;
  array[N_missing] int<lower=1, upper=N> ind_missing;
  matrix[M, J] X_missing;    
  array[N_missing] int<lower=1, upper=M> missing_covar_starts;
  array[N_missing] int<lower=1, upper=M> missing_covar_ends;
  vector<lower=0,upper=1>[M] p_missing_covar;
}

transformed data {
  array[N_missing] int len_missing_comb;
  
  for (i in 1:N_missing) {
    len_missing_comb[i] = missing_covar_ends[i] - missing_covar_starts[i] + 1;
  }
}

parameters {
  vector[J] beta;
  real<lower=0> sigma;
  vector[K] eta_raw;
  real<lower=0> sigma_re;
}

transformed parameters {
  vector[N_full] mu_full;
  vector[N_missing] ll_missing;
  
  if (N_full > 0) {
    // mean values for observations with all covariates
    mu_full = X[ind_full, ] * beta + eta_raw[map_to_study[ind_full]] * sigma_re;
  }
  
  if (N_missing > 0) {
    for (i in 1:N_missing) {
      int l = len_missing_comb[i];
      int k = ind_missing[i];
      vector[l] ll_tmp;
      int counter = 1;
      
      for (j in missing_covar_starts[i]:missing_covar_ends[i]) {
        real mu_tmp  = X_missing[j, ] * beta + eta_raw[map_to_study[k]] * sigma_re;
        ll_tmp[counter] = normal_lpdf(y[k]| mu_tmp, sigma) + log(p_missing_covar[j]);
        counter += 1;
      }
      ll_missing[i] = log_sum_exp(ll_tmp);
    }
  }  
  
}
model {
  
  y[ind_full] ~ normal(mu_full, sigma);
  target += sum(ll_missing);
  
  beta ~ normal(0, 1);
  eta_raw ~ std_normal();
  sigma ~ normal(0, 2);
  sigma_re ~ normal(0, 2);
}
