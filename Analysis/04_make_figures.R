# This script makes figures


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)


# Load results ------------------------------------------------------------

res_files <- dir("generated_data/", pattern = "stan_output", full.names = T)

effect_sizes <- map_df(res_files, function(x) {
  stan_fit <- readRDS(x)
  effect <- str_extract(x, get_all_effects() %>% str_c(collapse = "|"))
  exp_class <- str_extract(x, get_all_exp_class() %>% str_c(collapse = "|"))
  stan_data <- readRDS(make_stan_data_filename(effect = effect, exposure = exp_class))
  
  stan_fit$summary("beta") %>%
    mutate(var = colnames(stan_data$X),
           what = str_extract(var, get_all_covar_patterns()),
           strata = str_remove(var, what),
           strata = factor(strata, levels = unlist(map(stan_data$all_covars, ~ levels(stan_data$data[[.]]))))) %>%
    filter(!is.na(what)) %>% 
    mutate(what = map_chr(what, function(x) {
      ref <- levels(data[[x]])[1]
      str_c(x, "\nref:", ref)
    })) %>% 
    select(what, strata, mean, q5, q95) %>% 
    mutate(effect = effect,
           exp_class = exp_class)
})


pd <- position_dodge(width = .2)
p_effects <- effect_sizes %>%  
  filter(effect != "diff") %>%
  ggplot(aes(x = strata, y = mean)) +
  geom_hline(aes(yintercept = 0), lwd = .3, lty = 2) +
  geom_point(aes(color = effect), position = pd) +
  geom_errorbar(aes(ymin = q5, ymax = q95, color = effect), 
                position = pd, width = 0) +
  facet_grid(exp_class ~ what, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggthemes::scale_color_few("Dark") +
  labs(y = "effect size")


ggsave(p_effects, filename = "figures/meta_analysis_effect_sizes.png",
       width = 12, height = 6, dpi = 300)
