# This script makes figures


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)

source("Analysis/utils.R")

# Load results ------------------------------------------------------------

res_files <- dir("generated_data/", pattern = "stan_output", full.names = T)

es_effect <- get_effect_sizes(str_subset(res_files, "direction", negate = T))
es_direction <- get_effect_sizes(str_subset(res_files, "direction"))

plot_effect_sizes <- function(df) {
  pd <- position_dodge(width = .3)
  
  df %>%  
    filter(effect != "diff") %>%
    ggplot(aes(x = strata, y = mean)) +
    geom_hline(aes(yintercept = 0), lwd = .3, lty = 2) +
    geom_point(aes(color = effect), position = pd) +
    geom_errorbar(aes(ymin = q5, ymax = q95, color = effect), 
                  position = pd, width = 0) +
    facet_grid(exp_class ~ what, scales = "free", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    # scale_color_manual(values = c("blue", "red")) +
    ggthemes::scale_color_few("Dark") +
    labs(y = "effect size")
}


p_effects <- plot_effect_sizes(es_effect)

ggsave(p_effects, filename = "figures/meta_analysis_effect_sizes.png",
       width = 12, height = 6, dpi = 300)

p_dir <- plot_effect_sizes(es_direction) +
  scale_color_manual(values = c("blue", "red"))

ggsave(p_dir, filename = "figures/meta_analysis_effect_sizes_direction.png",
       width = 12, height = 6, dpi = 300)
