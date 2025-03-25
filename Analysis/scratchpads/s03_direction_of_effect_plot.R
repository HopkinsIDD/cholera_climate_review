# This script attempts to reproduce the effect of direction plot in 
# Boon and Thomson 2021 https://doi.org/10.1002/jrsm.1458


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(gt)    # first try using gt
library(here)

# Load data ---------------------------------------------------------------

# Estimate data
est_data <- readRDS(here("generated_data/processed_estimate_data.rds"))

# Study quality scores. For now taking the min score for a given study, 
# but may need to update
scores <- readRDS(here("generated_data/paper_scores.rds")) %>% 
  group_by(studykey) %>% 
  slice_min(tot_score, with_ties = F) %>%
  ungroup() %>% 
  distinct(studykey, tot_score) %>% 
  mutate(score_cat = cut(tot_score, c(0, 8, 10, Inf), include.lower = T, right = F),
         score_cat = factor(score_cat, labels = c("low", "medium", "high")))

thresh <- 0.7

# Compute statistics for direction of effect plot
dir_data <- est_data %>% 
  group_by(studykey, out_class_metric, exp_class, continent, admin_level) %>%
  # Apply the algorithm in Box 1 of Boon and Thomson
  summarise(N = n(),    # Number of estimates
            n_pos = sum(a_direction == "positive"),
            n_neg = sum(a_direction == "negative"),
            n_dir = n_pos + n_neg,
            overall_dir = case_when(n_pos >= thresh * n_dir ~ "positive",
                                    n_neg >= thresh * n_dir ~ "negative",
                                    TRUE ~ "uncertain") %>% 
              factor(levels = c("positive", "negative", "uncertain"),
                     labels = c("higher", "lower", "uncertain"))) %>% 
  ungroup() 


all_gts <- dir_data %>% 
  group_by(exp_class, out_class_metric) %>% 
  group_split() %>% 
  map(function(x) {
    
    dir_data2 <- x %>% 
      mutate(studykey = str_replace(studykey, "\xe6", "e"),
             dir_icon = case_when(overall_dir == "higher" ~ "caret-up",
                                  overall_dir == "lower" ~ "caret-down",
                                  TRUE ~ "caret-left,caret-right"),
             pubyear = str_extract(studykey, "[0-9]{4}") %>% as.numeric()) %>% 
      select(studykey, pubyear, exp_class, out_class_metric, continent, admin_level, dir_icon) %>% 
      inner_join(scores) %>% 
      arrange(dir_icon, score_cat, pubyear) 
    
    gt_list <- dir_data2 %>% 
      group_by(continent) %>% 
      group_split() %>% 
      map(function(y) {
        y %>% 
          select(studykey, admin_level, score_cat, dir_icon) %>% 
          gt() %>% 
          fmt_icon("dir_icon") %>% 
          tab_header(
            title = md(str_glue("exposure: **{y$exp_class[1]}**\ncontinent: **{y$continent[1]}**\noutcome: **{y$out_class_metric[1]}**")),
          ) %>% 
          tab_style(
            style = list(cell_fill(color = c("#EB9696"))),
            locations = cells_body(rows = score_cat == "low")
          ) %>% 
          tab_style(
            style = list(cell_fill(color = c("#E8AD7E"))),
            locations = cells_body(rows = score_cat == "medium")
          ) %>% 
          tab_style(
            style = list(cell_fill(color = c("#A5C295"))),
            locations = cells_body(rows = score_cat == "high")
          ) 
      })
    
    gt_group(.list = gt_list)  %>% 
      gtsave(filename = str_glue("figures/eff_dir_plot_exp-{x$exp_class[1]}_outcome-{x$out_class_metric[1]}.html"))
    
    gt_list
  })


all_gts %>% 
  unlist(recursive = F) %>% 
  gt_group(.list = .)  %>% 
  gtsave(filename = str_glue("figures/eff_dir_plot_all.html"))
