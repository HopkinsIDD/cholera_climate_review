# This script attempts to reproduce the effect of direction plot in 
# Boon and Thomson 2021 https://doi.org/10.1002/jrsm.1458


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(gt)
library(here)

source("Analysis/utils.R")

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
         score_cat = factor(score_cat, labels = c("low", "medium", "high")) %>% 
           forcats::fct_relevel(rev(c("low", "medium", "high"))))


# Compute statistics for direction of effect plot
dir_data <- est_data %>% 
  mutate(out_class_case_presentation = case_when(
    out_class_case_presentation == "not stated" ~ "unknown",
    !str_detect(out_class_case_presentation, "confirmed") ~ "suspected",
    !str_detect(out_class_case_presentation, "suspected") ~ "confirmed",
    TRUE ~ "mixed"
  )) %>% 
  group_by(studykey, out_class_metric, exp_class, continent, admin_level, out_class_case_presentation) %>%
  summarise_directions() %>% 
  ungroup()  %>% 
  mutate(studykey = str_replace(studykey, "\xe6", "e"),
         dir_icon = case_when(overall_dir == "higher" ~ "caret-up",
                              overall_dir == "lower" ~ "caret-down",
                              TRUE ~ "caret-left,caret-right"),
         overall_dir = factor(overall_dir, labels = c("Positive", "Negative", "Conflicting")), 
         pubyear = str_extract(studykey, "[0-9]{4}") %>% as.numeric()) 

# all_gts <- dir_data %>% 
#   group_by(exp_class, out_class_metric) %>% 
#   group_split() %>% 
#   map(function(x) {
#     
#     dir_data2 <- x %>% 
#       select(studykey, pubyear, exp_class, out_class_metric, continent, admin_level, dir_icon, overall_dir) %>% 
#       mutate(dir_icon_lvls = factor(dir_icon, levels = c("caret-up", "caret-left,caret-right", "caret-down"))) %>% 
#       inner_join(scores) %>% 
#       arrange(dir_icon_lvls, score_cat, pubyear) %>% 
#       group_by(dir_icon) %>% 
#       mutate(bottom_line = row_number() == n()) %>% 
#       ungroup()
#     
#     gt_list <- dir_data2 %>% 
#       group_by(continent) %>% 
#       group_split() %>% 
#       map(function(y) {
#         y %>% 
#           select(studykey, admin_level, score_cat, dir_icon, bottom_line, overall_dir) %>% 
#           gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
#           fmt_icon("dir_icon") %>% 
#           tab_header(
#             title = md(str_glue("exposure: **{y$exp_class[1]}**\ncontinent: **{y$continent[1]}**\noutcome: **{y$out_class_metric[1]}**")),
#           ) %>% 
#           add_doe_styling()
#       })
#     
#     gt_group(.list = gt_list)  %>% 
#       gtsave(filename = str_glue("figures/eff_dir_plot_exp-{x$exp_class[1]}_outcome-{x$out_class_metric[1]}.html"))
#     
#     # Also save with all combined
#     gt_comb <- dir_data2 %>% 
#       select(studykey, continent, admin_level, score_cat, dir_icon, bottom_line, overall_dir) %>% 
#       gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
#       fmt_icon("dir_icon") %>% 
#       tab_header(
#         title = md(str_glue("exposure: **{x$exp_class[1]}**\noutcome: **{x$out_class_metric[1]}**")),
#       ) %>% 
#       add_doe_styling()
#     
#     gt_comb  %>% 
#       gtsave(filename = str_glue("figures/eff_dir_plot_exp-{x$exp_class[1]}_outcome-{x$out_class_metric[1]}_combined.html"))
#     
#     gt_list
#   })
# 
# 
# all_gts %>% 
#   unlist(recursive = F) %>% 
#   gt_group(.list = .)  %>% 
#   gtsave(filename = str_glue("figures/eff_dir_plot_all.html"))
# 
# 
# # Attempt at combining things
# dir_dat3 <- dir_data %>% 
#   select(studykey, pubyear, exp_class, out_class_metric, continent, admin_level, dir_icon, overall_dir) %>% 
#   mutate(dir_icon_lvls = factor(dir_icon, levels = c("caret-up", "caret-left,caret-right", "caret-down"))) %>% 
#   inner_join(scores) %>% 
#   arrange(dir_icon_lvls, score_cat, pubyear) %>% 
#   group_by(dir_icon) %>% 
#   mutate(bottom_line = row_number() == n()) %>% 
#   ungroup() %>% 
#   select(studykey, pubyear, exp_class, out_class_metric, continent, admin_level, dir_icon, overall_dir, score_cat, bottom_line) %>% 
#   filter(out_class_metric == "occurrence") %>% 
#   pivot_wider(names_from = "exp_class",
#               values_from = "dir_icon") %>% 
#   select(studykey, pubyear, continent, admin_level, overall_dir, score_cat, bottom_line, 
#          precipitation, temperature, flood, drought) %>% 
#   arrange(score_cat, desc(precipitation), desc(temperature), desc(flood), desc(drought),pubyear) %>% 
#   select(-pubyear) 
# 
# gt1 <- dir_dat3 %>%
#   filter(overall_dir == "Positive") %>% 
#   gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
#   fmt_icon(c("precipitation", "temperature", "flood", "drought")) %>% 
#   add_doe_styling()  
# 
# 
# gt2 <- dir_dat3 %>%
#   filter(overall_dir != "Positive") %>% 
#   gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
#   fmt_icon(c("precipitation", "temperature", "flood", "drought")) %>% 
#   add_doe_styling()
# 
# 
# gtExtras::gt_two_column_layout(tables = list(gt1, gt2), 
#                                output = "save",
#                                filename = "figures/eff_dir_plot_all_exp_classes_occurrence.png",
#                                vwidth = 1450, vheight = 820)
# 
# # gtsave(filename = str_glue("figures/eff_dir_plot_all_exp_classes.html"))

# Attempt at combining two exposures
make_dir_data <- function(df, 
                          outcome, 
                          covars = c("continent", "admin_level", "out_class_case_presentation")) {
  df %>% 
    # select(studykey, pubyear, exp_class, out_class_metric, continent, admin_level, dir_icon, overall_dir) %>% 
    mutate(dir_icon_lvls = factor(dir_icon, levels = c("caret-up", "caret-left,caret-right", "caret-down"))) %>% 
    inner_join(scores) %>% 
    arrange(dir_icon_lvls, score_cat, pubyear) %>% 
    group_by(dir_icon) %>% 
    mutate(bottom_line = row_number() == n()) %>% 
    # ungroup() %>% 
    # group_by(studykey, out_class_metric) %>%
    # mutate(overall_dir = overall_dir[!is.na(overall_dir)][1]) %>% 
    # ungroup() %>% 
    select(studykey, pubyear, exp_class, out_class_metric, all_of(covars), dir_icon, overall_dir, score_cat, bottom_line) %>% 
    filter(out_class_metric == outcome) %>% 
    pivot_wider(names_from = "exp_class",
                values_from = "dir_icon") %>% 
    select(studykey, pubyear, all_of(covars), overall_dir, score_cat, bottom_line, 
           precipitation, temperature, flood, drought) %>% 
    mutate(overall_dir = overall_dir %>% forcats::fct_relevel("Positive", "Conflicting", "Negative")) %>% 
    arrange(overall_dir, score_cat, desc(precipitation), desc(temperature), desc(flood), desc(drought),pubyear) %>% 
    select(-pubyear)
}

dir_dat_occ <- make_dir_data(df = dir_data, outcome = "occurrence")
dir_dat_incid <- make_dir_data(df = dir_data, outcome = "incidence")


gt1 <- dir_dat_incid %>%
  gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
  fmt_icon(c("precipitation", "temperature", "flood", "drought")) %>% 
  add_doe_styling() %>% 
  tab_header(title = "Incidence")

gt2 <- dir_dat_occ %>%
  gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
  fmt_icon(c("precipitation", "temperature", "flood", "drought")) %>% 
  add_doe_styling() %>% 
  tab_header(title = "Occurrence")

gtExtras::gt_two_column_layout(tables = list(gt1, gt2), 
                               output = "save",
                               filename = "figures/eff_dir_plot_all_exp_classes_outcomes.png",
                               vwidth = 860,
                               zoom = 1)#,
# vheight = 500)


# Try by exposure ---------------------------------------------------------

# Attempt at combining two exposures
make_dir_data2 <- function(df, 
                           exposure, 
                           covars = c("continent", "admin_level", "out_class_case_presentation")) {
  df %>% 
    # select(studykey, pubyear, exp_class, out_class_metric, continent, admin_level, dir_icon, overall_dir) %>% 
    mutate(dir_icon_lvls = factor(dir_icon, levels = c("caret-up", "caret-left,caret-right", "caret-down"))) %>% 
    inner_join(scores) %>% 
    arrange(dir_icon_lvls, score_cat, pubyear) %>% 
    group_by(dir_icon) %>% 
    mutate(bottom_line = row_number() == n()) %>% 
    # ungroup() %>% 
    # group_by(studykey, out_class_metric) %>%
    # mutate(overall_dir = overall_dir[!is.na(overall_dir)][1]) %>% 
    # ungroup() %>% 
    select(studykey, pubyear, exp_class, out_class_metric, all_of(covars), dir_icon, overall_dir, score_cat, bottom_line) %>% 
    filter(exp_class == exposure) %>% 
    pivot_wider(names_from = "out_class_metric",
                values_from = "dir_icon") %>% 
    select(studykey, pubyear, all_of(covars), overall_dir, score_cat, bottom_line, 
           incidence, occurrence) %>% 
    mutate(overall_dir = overall_dir %>% forcats::fct_relevel("Positive", "Conflicting", "Negative")) %>% 
    arrange(overall_dir, score_cat, desc(incidence), desc(occurrence), pubyear) %>% 
    select(-pubyear)
}

dir_dat_precip <- make_dir_data2(df = dir_data, exposure = "precipitation")
gt_p <- dir_dat_precip %>%
  gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
  fmt_icon(c("incidence", "occurrence")) %>% 
  add_doe_styling2() %>% 
  tab_header(title = "Precipitation")

dir_dat_precip <- make_dir_data2(df = dir_data, exposure = "precipitation")
gt_p <- dir_dat_precip %>%
  gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
  fmt_icon(c("incidence", "occurrence")) %>% 
  add_doe_styling2() %>% 
  tab_header(title = "Precipitation")
gtsave(gt_p, filename = "figures/dirplot_precipitation.png")


dir_dat_temp <- make_dir_data2(df = dir_data, exposure = "temperature")
gt_t <- dir_dat_temp %>%
  gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
  fmt_icon(c("incidence", "occurrence")) %>% 
  add_doe_styling2() %>% 
  tab_header(title = "Temperature")
gtsave(gt_t, filename = "figures/dirplot_temperature.png")


dir_dat_flood <- make_dir_data2(df = dir_data, exposure = "flood")
gt_f <- dir_dat_flood %>%
  gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
  fmt_icon(c("incidence", "occurrence")) %>% 
  add_doe_styling2() %>% 
  tab_header(title = "Flood")
gtsave(gt_f, filename = "figures/dirplot_flood.png")

dir_dat_drought <- make_dir_data2(df = dir_data, exposure = "drought")
gt_d <- dir_dat_drought %>%
  gt(groupname_col = "overall_dir", rowname_col = "studykey") %>% 
  fmt_icon(c("incidence", "occurrence")) %>% 
  add_doe_styling2() %>% 
  tab_header(title = "Drougth")
gtsave(gt_d, filename = "figures/dirplot_drought.png")

library(patchwork)
library(magick)

# p_exposures <- wrap_table(gt_p) | (wrap_table(gt_t) / wrap_table(gt_f) / wrap_table(gt_d, panel = "full"))

pp <- image_read("figures/dirplot_precipitation.png") %>% image_ggplot() + theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"))
pt <- image_read("figures/dirplot_temperature.png") %>% image_ggplot() + theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"))
pf <- image_read("figures/dirplot_flood.png") %>% image_ggplot() + theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"))
pd <- image_read("figures/dirplot_drought.png") %>% image_ggplot()+ theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"))

p_exposures <- (pp + pt + pf + pd) + 
  plot_layout(
    design = 
      "AB
    AB
    AC
    AD")

ggsave(p_exposures, filename = "figures/eff_dir_plot_exposures_by_outcome.png",
       width = 11,
       height = 18)

