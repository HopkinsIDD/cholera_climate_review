//
// Performs correlation meta-analysis


data {
  int<lower=0> N;    // number of observations
  int<lower=0> J;    // number of covariates
  int<lower=0> K;    // number of studies
  int<lower=0> N_full;      // number of observations for which all covariates are available
  int<lower=0> N_missing;   // number of observations with missing covariates
  matrix[N, J] X;
  int<lower=0> L;    // number of classes of precision
  
  vector[N] y;    // observations
  
  array[N] int<lower=1, upper=K> map_to_study;    // map from observation to study
  array[N] int<lower=1, upper=L> map_to_weight;    // map from observation to study
  array[N_full] int<lower=1, upper=N> ind_full;
  array[N_missing] int<lower=1, upper=N> ind_missing;
  array[N_missing, J] int<lower=0, upper=J> missing_covar_levels;    // map to index of covariates to marginalize over
  
  vector<lower=0,upper=1>[J] p_missing_covar;
}


parameters {
  vector[J] beta;
  positive_ordered[L] sigma;
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
      real mu;   // mean
      int k = ind_missing[i];
      vector[J] ll_tmp = rep_vector(0, J);
      int counter = 1;
      // Value without missing covariates (X is 0 in places where covariate is missing)
      mu = X[k, ] * beta + eta_raw[map_to_study[k]] * sigma_re;
      
      ll_tmp[1] = normal_lpdf(y[k]| mu, sigma[map_to_weight[k]]);
      
      for (j in 1:J) {
        real mu_tmp = mu;
        int l = missing_covar_levels[i, j];
        real tot_p = 0;
        
        if (l > 0) {
          mu_tmp += beta[l];
          counter += 1;
          tot_p += p_missing_covar[l];
          ll_tmp[counter] += normal_lpdf(y[k]| mu_tmp, sigma[map_to_weight[k]]) + log(p_missing_covar[l]);
        }
        
        ll_tmp[1] += log1m(tot_p);
      }
      
      ll_missing[i] = log_sum_exp(ll_tmp[1:counter]);
    }
  }  
  
}
model {
  
  y[ind_full] ~ normal(mu_full, sigma[map_to_weight[ind_full]]);
  target += sum(ll_missing);
  
  beta ~ normal(0, 1);
  eta_raw ~ std_normal();
  sigma ~ normal(0, 2);
  sigma_re ~ normal(0, 2);
}

generated quantities {
  vector[L] I_squared;
  
  I_squared = sigma_re^2 / (sigma_re^2 + sigma.^2);
}
