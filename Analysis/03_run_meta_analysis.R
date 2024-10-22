# This script runs the stan code for the meta-analysis

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(here)
library(optparse)
library(cmdstanr)

source("Analysis/utils.R")

# User-supplied options
option_list <- list(
  make_option(c("-e", "--effect"), 
              default = "corr", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-c", "--exp_class"), 
              default = "precipitation", action ="store", type = "character", help = "Make intermediate plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Load stan data ----------------------------------------------------------

stan_data <- readRDS(make_stan_data_filename(effect = opt$effect, exposure = opt$exp_class))

# Compile and run stan model ----------------------------------------------

stan_model <- cmdstan_model("Analysis/stan/correlation_meta_analysis_marginalize_weighted.stan")

stan_fit <-  stan_model$sample(
  data = stan_data$stan_data,
  iter_warmup = 100,
  iter_sampling = 100,
  parallel_chains = 4,
  init = .1,
  max_treedepth = 10,
  refresh = 10,
  save_warmup = FALSE
)

stan_fit$save_object(make_output_filename(effect = opt$effect, exposure = opt$exp_class))
