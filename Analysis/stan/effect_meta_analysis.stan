//
// Performs correlation meta-analysis


data {
  int<lower=0> N;    // number of observations
  int<lower=0> J;    // number of covariates
  int<lower=0> K;    // number of studies
  matrix[N, J] X;
  vector[N] y;
  array[N] int<lower=1, upper=K> map_to_study;    // map from observation to study
}

parameters {
  vector[J] beta;
  real<lower=0> sigma;
  vector[K] eta_raw;
  real<lower=0> sigma_re;
}

model {
  vector[N] mu = X * beta + eta_raw[map_to_study] * sigma_re;
  
  y ~ normal(mu, sigma);
  beta ~ normal(0, 1);
  eta_raw ~ std_normal();
  sigma ~ normal(0, 2);
  sigma_re ~ normal(0, 2);
}
