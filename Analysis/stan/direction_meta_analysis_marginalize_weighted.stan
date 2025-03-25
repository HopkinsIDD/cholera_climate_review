// Performs direction meta-analysis


data {
  int<lower=0> N;    // number of observations
  int<lower=0> J;    // number of covariates
  int<lower=0> K;    // number of studies
  int<lower=0> N_full;      // number of observations for which all covariates are available
  int<lower=0> N_missing;   // number of observations with missing covariates
  matrix[N, J] X;
  int<lower=0> L;    // number of classes of precision
  int<lower=0> M;    // number of missing covariates
  
  array[N] int<lower=0, upper=1> y;    // observations
  
  array[N] int<lower=1, upper=K> map_to_study;    // map from observation to study
  array[N] int<lower=1, upper=L> map_to_weight;    // map from observation to study
  array[N_full] int<lower=1, upper=N> ind_full;
  array[N_missing] int<lower=1, upper=N> ind_missing;
  matrix[M, J] X_missing;    
  array[N_missing] int<lower=1, upper=M> missing_covar_starts;
  array[N_missing] int<lower=1, upper=M> missing_covar_ends;
  vector<lower=0,upper=1>[M] p_missing_covar;
  
  // post-stratification
  int U_poststrat;
  int U_covar;
  matrix[U_poststrat, U_covar] poststrat_mat_comb;
  matrix[U_covar, J] u_covar_mat;
}

transformed data {
  array[N_missing] int len_missing_comb;
  
  for (i in 1:N_missing) {
    len_missing_comb[i] = missing_covar_ends[i] - missing_covar_starts[i] + 1;
  }
}

parameters {
  vector[J] beta;
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
      real alpha_tmp  = eta_raw[map_to_study[k]] * sigma_re;
      
      for (j in missing_covar_starts[i]:missing_covar_ends[i]) {
        ll_tmp[counter] =  bernoulli_logit_glm_lpmf(y[k]|  to_matrix(X_missing[j, ]), alpha_tmp, beta) + log(p_missing_covar[j]);
        counter += 1;
      }
      ll_missing[i] = log_sum_exp(ll_tmp);
    }
    
  }  
  
}
model {
  vector[N_full] alpha_full;
  alpha_full = eta_raw[map_to_study[ind_full]] * sigma_re;
  target += bernoulli_logit_glm_lpmf(y[ind_full]| X[ind_full, ], alpha_full, beta);
  target += sum(ll_missing);
  
  beta ~ normal(0, .5);
  eta_raw ~ std_normal();
  sigma_re ~ normal(0, 2);
}

generated quantities {
  vector[U_poststrat] poststrat;
  
  poststrat = inv_logit(to_vector(poststrat_mat_comb * u_covar_mat * beta));
}
