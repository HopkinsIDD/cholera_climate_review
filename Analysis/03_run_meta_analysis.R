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
              default = "direction_positive", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-c", "--exp_class"), 
              default = "precipitation", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-o", "--out_class"), 
              default = "incidence", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-l", "--level"), 
              default = "Africa", action ="store", type = "character", help = "Stratification level"),
  make_option(c("-s", "--stratification"), 
              default = "continent", action ="store", type = "character", help = "Stratification column"),
  make_option(c("-r", "--redo"), 
              default = FALSE, action = "store", type = "logical", help = "Redo fit")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)


# Load stan data ----------------------------------------------------------

if (!file.exists(make_stan_data_filename(opt = opt))) {
  # 02 found no data for this stratum
  cat("-- No stan data for", make_stan_data_filename(opt = opt), "\n")
  quit(save = "no", status = 0)
}

stan_data <- readRDS(make_stan_data_filename(opt = opt))

# Inverse Fisher z-transform of the post-stratified estimates for correlations
stan_data$stan_data$transform <- as.integer(opt$effect == "corr")

# Compile and run stan model ----------------------------------------------

if (!str_detect(opt$effect, "direction")) {
  stan_model <- cmdstan_model("Analysis/stan/effect_meta_analysis_marginalize_weighted.stan")
} else {
  stan_model <- cmdstan_model("Analysis/stan/direction_meta_analysis_marginalize_weighted.stan")
}

if (file.exists(make_stan_fit_filename(opt = opt)) & !opt$redo) {
  cat("-- Stan fit already done for", make_stan_fit_filename(opt = opt), "\n")
  
  stan_fit <- readRDS(make_stan_fit_filename(opt = opt))
} else {
  
  stan_fit <-  stan_model$sample(
    data = stan_data$stan_data,
    iter_warmup = 250,
    iter_sampling = 500,
    parallel_chains = 4,
    init = .1,
    max_treedepth = 10
  )
  
  stan_fit$save_object(make_stan_fit_filename(opt = opt))
} 


# Genquant ----------------------------------------------------------------

genquant <- stan_model$generate_quantities(fitted_params = stan_fit,
                                           data = stan_data$stan_data,
                                           parallel_chains = 4)

genquant$save_object(make_stan_genquant_filename(opt = opt))

# Scraps ------------------------------------------------------------------

# 
# stan_fit$summary(c("sigma_re", "beta"))
# stan_fit$draws(c("beta")) %>%
#   bayesplot::mcmc_intervals()
# stan_fit$draws(c("poststrat_pred", "poststrat_pred_mean")) %>% 
#   bayesplot::mcmc_intervals()

# p_corr <- stan_fit$summary("poststrat") %>%
#   mutate(var = stan_data$postrat_names,
#          what = str_split(var, "-") %>% map_chr(~.[1]),
#          strata = str_split(var, "-") %>% map_chr(~.[2]),
#          strata = factor(strata, levels = stan_data$postrat_names %>% str_split("-") %>% map_chr(~.[2]))) %>%
#   ggplot(aes(x = strata, y = mean)) +
#   geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
#   geom_point() +
#   geom_hline(aes(yintercept = 0), lty = 2, lwd = .5) +
#   facet_grid(. ~ what, scales = "free", space = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   labs(y = "correlation")
# 
# ggsave(p_corr, filename = "figures/poststrat_corrlation.png", width = 12, height = 5) 

# stan_fit$summary("eta_re")
# stan_fit$summary("sigma_re")
# 
# pd <- position_dodge(.2)
# p_pos <- stan_fit$summary("poststrat") %>%
#   mutate(set = "mean",
#          var = stan_data$postrat_names) %>%
#   bind_rows(stan_fit$summary("poststrat_pred") %>%
#               mutate(set = "pred",
#                      var = stan_data$postrat_names)) %>%
#   mutate(what = str_split(var, "-") %>% map_chr(~.[1]),
#          strata = str_split(var, "-") %>% map_chr(~.[2]),
#          strata = factor(strata, levels = stan_data$postrat_names %>% str_split("-") %>% map_chr(~.[2]))) %>%
#   ggplot(aes(x = strata, y = mean)) +
#   geom_errorbar(aes(ymin = q5, ymax = q95, color = set, group = set), width = 0, position = pd) +
#   geom_point(position = pd, aes(color = set)) +
#   geom_hline(aes(yintercept = .5), lty = 2, lwd = .5) +
#   facet_grid(. ~ what, scales = "free", space = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   labs(y = "probability of positive direction")
# 
# p_pos

# ggsave(p_pos, filename = "figures/poststrat_positive_direction.png", width = 12, height = 5) 


