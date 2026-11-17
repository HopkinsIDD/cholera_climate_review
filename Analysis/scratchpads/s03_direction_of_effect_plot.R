# This script attempts to reproduce the effect of direction plot in
# Boon and Thomson 2021 https://doi.org/10.1002/jrsm.1458
# The manuscript versions (figure_3.png, supfig_eff_dir_plot_*.png) are made in
# Analysis/04_make_figures.R; this scratchpad writes to figures/s03_* only.


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(gt)
library(here)
library(patchwork)

source("Analysis/utils.R")

# Load data ---------------------------------------------------------------

est_data <- readRDS(here("generated_data/processed_estimate_data.rds"))
scores <- get_score_categories()
dir_data <- make_effdir_data(est_data)

# Try by exposure ---------------------------------------------------------

# Make overall plot
compose_effdir_plot(data = dir_data,
                    scores = scores,
                    filename = "figures/s03_eff_dir_plot.png")

# Make plots by continent
walk(c("Asia", "Africa", "Other"), function(x) {
  compose_effdir_plot(data = dir_data %>% filter(continent == x),
                      scores = scores,
                      filename = str_glue("figures/s03_eff_dir_plot_{x}.png"))
})
