// Performs direction meta-analysis

// To integrate out variability
functions {
  real to_integrate(
    real x,          // Function argument
    real xc,         // Complement of function argument
    //  on the domain (defined later)
    array[] real theta,    // parameters
    array[] real x_r,      // data (real)
    array[] int x_i) {     // data (integer)
    
    real mu = theta[1];
    real sigma = theta[2];
    real res;
    
    res = inv_logit(mu + sigma * inv_Phi(x));
    return res;
    }
}
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
  vector[N] epsilon_raw;
  real<lower=0> sigma_epsilon;
}

transformed parameters {
  vector[N_missing] ll_missing;
  
  if (N_missing > 0) {
    
    for (i in 1:N_missing) {
      int l = len_missing_comb[i];
      int k = ind_missing[i];
      vector[l] ll_tmp;
      int counter = 1;
      real alpha_tmp  = eta_raw[map_to_study[k]] * sigma_re + epsilon_raw[k] * sigma_epsilon;
      
      for (j in missing_covar_starts[i]:missing_covar_ends[i]) {
        ll_tmp[counter] = bernoulli_logit_glm_lpmf(y[k]|  to_matrix(X_missing[j, ]), alpha_tmp, beta) + log(p_missing_covar[j]);
        counter += 1;
      }
      ll_missing[i] = log_sum_exp(ll_tmp);
    }
    
  }  
  
}
model {
  vector[N_full] alpha_full;
  alpha_full = eta_raw[map_to_study[ind_full]] * sigma_re + epsilon_raw[ind_full] * sigma_epsilon;
  target += bernoulli_logit_glm_lpmf(y[ind_full]| X[ind_full, ], alpha_full, beta);
  target += sum(ll_missing);
  
  beta ~ normal(0, .5);
  eta_raw ~ std_normal();
  epsilon_raw ~ std_normal();
  sigma_re ~ normal(0, 1);
  sigma_epsilon ~ normal(0, 1);
}

generated quantities {
  vector[U_poststrat] poststrat;    // postsratified mean
  // vector[U_poststrat] poststrat_pred_mean;    // expectation of postsratified prediction interval mean 
  vector[U_poststrat] poststrat_pred;    // postsratified prediction interval mean
  vector[J] eff_size = exp(beta);
  // Note: with one Bernoulli draw per estimate, sigma_epsilon is weakly
  // identified, so i_squared largely reflects the priors on the two sds.
  // L and map_to_weight are passed in but not used by this model.
  real i_squared = sigma_re^2 / (sigma_re^2 + sigma_epsilon^2);
  
  {
    vector[U_covar] eta;          // logit-linear perdictor 
    vector[U_covar] eta_re;
    
    eta = u_covar_mat * beta;
    poststrat = poststrat_mat_comb * inv_logit(eta);
    
    // for (i in 1:U_covar) {
      //   array[0] real x_r;
      //   array[0] int x_i;
      //   array[2] real par;
      //   par[1] = eta[i];
      //   par[2] = sigma_re;
      //   eta_re[i] = integrate_1d(to_integrate, 0, 1, par, x_r, x_i, 1e-8);
      // }
      
      poststrat_pred = poststrat_mat_comb * inv_logit(eta + normal_rng(0, sigma_re));
      // poststrat_pred_mean = poststrat_mat_comb * eta_re;
  }
}

