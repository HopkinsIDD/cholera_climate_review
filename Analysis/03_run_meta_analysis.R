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
              default = "temperature", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-o", "--out_class"), 
              default = "incidence", action ="store", type = "character", help = "Make intermediate plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Load stan data ----------------------------------------------------------

stan_data <- readRDS(make_stan_data_filename(opt = opt))

# Compile and run stan model ----------------------------------------------

if (!str_detect(opt$effect, "direction")) {
  stan_model <- cmdstan_model("Analysis/stan/effect_meta_analysis_marginalize_weighted.stan")
} else {
  stan_model <- cmdstan_model("Analysis/stan/direction_meta_analysis_marginalize_weighted.stan")
}

stan_fit <-  stan_model$sample(
  data = append(
    stan_data$stan_data, 
    list(transform = 1)
  ),
  iter_warmup = 100,
  iter_sampling = 100,
  parallel_chains = 4,
  init = .1,
  max_treedepth = 10
)

stan_fit$save_object(make_stan_fit_filename(opt = opt))

# p_corr <- stan_fit$summary("poststrat") %>% 
  # mutate(var = stan_data$postrat_names,
  #        what = str_split(var, "-") %>% map_chr(~.[1]),
  #        strata = str_split(var, "-") %>% map_chr(~.[2]),
  #        strata = factor(strata, levels = stan_data$postrat_names %>% str_split("-") %>% map_chr(~.[2]))) %>%
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
# 
# p_pos <- stan_fit$summary("poststrat") %>% 
#   mutate(var = stan_data$postrat_names,
#          what = str_split(var, "-") %>% map_chr(~.[1]),
#          strata = str_split(var, "-") %>% map_chr(~.[2]),
#          strata = factor(strata, levels = stan_data$postrat_names %>% str_split("-") %>% map_chr(~.[2]))) %>%
#   ggplot(aes(x = strata, y = mean)) +
#   geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) + 
#   geom_point() +
#   geom_hline(aes(yintercept = .5), lty = 2, lwd = .5) +
#   facet_grid(. ~ what, scales = "free", space = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   labs(y = "probability of positive direction")
# 
# ggsave(p_pos, filename = "figures/poststrat_positive_direction.png", width = 12, height = 5) 
