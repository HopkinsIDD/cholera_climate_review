# This script makes figures


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)

source("Analysis/utils.R")

# Load results ------------------------------------------------------------

res_files <- dir("generated_data/", pattern = "stan_fit", full.names = T)

es_effect <- get_effect_sizes(str_subset(res_files, "corr"))
# es_direction <- get_effect_sizes(str_subset(res_files, "direction"))

plot_effect_sizes <- function(df) {
  pd <- position_dodge(width = .3)
  
  df %>%  
    ggplot(aes(x = strata, y = mean)) +
    geom_hline(aes(yintercept = 0), lwd = .3, lty = 2) +
    geom_point(aes(color = out_class), position = pd) +
    geom_errorbar(aes(ymin = q5, ymax = q95, color = out_class), 
                  position = pd, width = 0) +
    facet_grid(exp_class ~ what, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    # scale_color_manual(values = c("blue", "red")) +
    ggthemes::scale_color_few("Dark") +
    labs(y = "effect size")
}


p_effects <- plot_effect_sizes(es_effect)

ggsave(p_effects, filename = "figures/meta_analysis_effect_sizes.png",
       width = 12, height = 6, dpi = 300)

# p_dir <- plot_effect_sizes(es_direction) +
#   scale_color_manual(values = c("blue", "red"))
# 
# ggsave(p_dir, filename = "figures/meta_analysis_effect_sizes_direction.png",
#        width = 12, height = 6, dpi = 300)



plot_poststrat <- function(df, 
                           yintercept = 0,
                           y_title = "") {
  pd <- position_dodge(width = .4)
  
  df %>%  
    ggplot(aes(x = strata, y = mean)) +
    geom_hline(aes(yintercept = yintercept), lwd = .3, lty = 2) +
    geom_point(aes(color = exp_class, pch = out_class), position = pd) +
    geom_errorbar(aes(ymin = q5, ymax = q95, color = exp_class, lty = out_class), 
                  position = pd, width = 0, alpha = .7, lwd = .3) +
    facet_grid(exp_class ~ what, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    # scale_color_manual(values = c("blue", "red")) +
    ggthemes::scale_color_few("Dark") +
    scale_linetype_manual(values = c(1, 6)) +
    labs(y = y_title)
}

# Correlation
corr_poststrat <- get_poststrat(str_subset(res_files, "corr"))
p_poststrat <- plot_poststrat(corr_poststrat, yintercept = 0, y_title = "post-stratified correlation")

ggsave(p_poststrat, filename = "figures/meta_analysis_poststrat_correlation.png",
       width = 12, height = 6, dpi = 300)

# Direction
dir_poststrat <- get_poststrat(str_subset(res_files, "direction"))
p_poststrat <- plot_poststrat(dir_poststrat, yintercept = 0.5, y_title = "post-stratified probability\nof positive association")

ggsave(p_poststrat, filename = "figures/meta_analysis_poststrat_direction.png",
       width = 12, height = 6, dpi = 300)
