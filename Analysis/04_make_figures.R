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
    geom_point(aes(color = exp_class, pch = out_class), position = pd) +
    geom_errorbar(aes(ymin = q5, ymax = q95, color = exp_class, lty = out_class), 
                  position = pd, width = 0) +
    facet_grid(exp_class ~ what, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    # scale_color_manual(values = c("blue", "red")) +
    ggthemes::scale_color_few("Dark") +
    scale_linetype_manual(values = c(1, 6)) +
    labs(y = "effect size")
}


p_effects <- plot_effect_sizes(es_effect)

ggsave(p_effects, filename = "figures/meta_analysis_effect_sizes.png",
       width = 12, height = 6, dpi = 300)


es_effect_dir <- get_effect_sizes(str_subset(res_files, "dir"))
p_effects_dir <- plot_effect_sizes(es_effect_dir)

ggsave(p_effects_dir, filename = "figures/meta_analysis_effect_sizes_direction.png",
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

# Fig 2. Spatial-temporal coverage of studies ----
## Fig 2a Spatial coverage map ----
plot_data <- read.csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv") %>% 
  rename(admin_level = new_admin) %>% 
  mutate(admin_level = case_when(
    admin_level == "subnational" ~ "sub-national",
    admin_level == "supranational" ~ "regional",
    TRUE ~ "national"
  ))

### Identify cholera endemic countries ----
report<-read.csv("Data/endemicity/who_annual_report.csv") %>% 
  tidyr::separate(obs,c('time_left','time_right','suspected_cases','deaths','location_id'),sep = ';',remove=T) %>% 
  merge(.,read.csv("Data/endemicity/country_location_id_name.csv"),by.x="location_id",by.y="id",all.x=TRUE) %>% 
  filter(is.na(region)==F) %>% 
  arrange(region) %>% 
  mutate(region = as.factor(region),suspected_cases=as.numeric(suspected_cases))

plot_data <- plot_data %>%
  rowwise() %>% 
  mutate(
    year_tl = ifelse(nchar(eff_tl) == 5, as.numeric(lubridate::year(as.Date(as.numeric(eff_tl),origin="1899-12-30"))), as.numeric(substr(eff_tl,1,4))),
    year_tr = ifelse(nchar(eff_tr) == 5, as.numeric(lubridate::year(as.Date(as.numeric(eff_tr),origin="1899-12-30"))), as.numeric(substr(eff_tr,1,4)))
  ) %>% 
  mutate(year_tl = lubridate::ymd(paste0(year_tl, "-01-01")),
         year_tr = lubridate::ymd(paste0(year_tr, "-12-31"))) %>%
  left_join(report %>% select(region) %>% distinct(), by = c("eff_country" = "region")) %>% 
  mutate(eff_country = as.factor(eff_country))

report <- report %>% 
  filter(region %in% plot_data$eff_country) %>% 
  mutate(region = factor(region,levels = levels(plot_data$eff_country))) %>% 
  group_by(region,time_left,time_right,location_id) %>%
  summarize(suspected_cases = mean(suspected_cases,na.rm = T),
            deaths = mean(deaths,na.rm=T))

plot_data_with_endemicity <- identify_endemicity(
  processed_data = plot_data %>%
    mutate(year_tl = 2025) %>% 
    separate_rows(eff_country, sep = "\\|"), # to examine the endemicity since 2000
  study_year_col = "year_tl",
  study_year_lb = 2000,
  study_country_col = "eff_country",
  who_data = report %>% mutate(year = lubridate::year(time_left)),
  who_year_col = "year",
  who_country_col = "region",
  who_case_col = "suspected_cases",
  prior_year_with_sCh = 11
)

fig2a <- plot_spatial_distribution_map(plot_data = plot_data_with_endemicity,
                                       country_col="eff_country",
                                       group_col = "admin_level",
                                       group_levels = c("regional","national","sub-national"),
                                       endemicity_col = "endemicity_status")

## Fig 2b Temporal ----
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

plot_data <- read.csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv") %>% 
  rowwise() %>% 
  mutate(
    year_tl = ifelse(nchar(eff_tl) == 5, as.character(as.Date(as.numeric(eff_tl),origin="1899-12-30")), as.character(paste0(as.numeric(substr(eff_tl,1,4)),"-01-01"))),
    year_tr = ifelse(nchar(eff_tr) == 5, as.character(as.Date(as.numeric(eff_tr),origin="1899-12-30")),  as.character(paste0(as.numeric(substr(eff_tr,1,4)),"-01-01")))
  ) %>%
  mutate(
    year_tl = ymd(year_tl),
    year_tr = ymd(year_tr),
    duration_days = as.numeric(year_tr - year_tl),
    start_year = year(year_tl),
    start_date = lubridate::make_date(start_year, 1, 1)  
  ) %>% 
  rename(admin_level = new_admin) %>% 
  mutate(admin_level = case_when(
    admin_level == "subnational" ~ "sub-national",
    admin_level == "supranational" ~ "regional",
    TRUE ~ "national"
  )) %>% 
  add_continent(country_col="eff_country")


continent_colors <- c(
  "Africa" = "#FFA378",
  "Asia" = "#A8B545",
  "Other" = "#8C796D"
)

timeline_plot <- ggplot(plot_data %>% arrange(year_tl)) +
  geom_segment(aes(x = year_tl, xend = year_tr, 
                   y = reorder(studykey, year_tl), 
                   yend = reorder(studykey, year_tl), 
                   color = continent),
               linewidth = 0.8) +
  scale_color_manual(values = continent_colors) +
  scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
  labs(x = "Year", y = "Study", color = "Continent") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0)
  )

hist_top <- ggplot(plot_data %>% distinct(studykey, start_date, continent)) +
  geom_histogram(aes(x = start_date, fill = continent,color = continent), 
                 binwidth = 365, position = "stack") +
  scale_fill_manual(values = continent_colors) +
  scale_color_manual(values = continent_colors) +
  scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
  labs(x = NULL, y = "Study count", fill = "Continent") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0)
  )

hist_bottom <- ggplot(plot_data %>% distinct(studykey,duration_days,continent)) +
  geom_histogram(aes(x = duration_days, fill = continent, color = continent), binwidth = 365, position = "stack") +
  scale_fill_manual(values = continent_colors) +
  scale_color_manual(values = continent_colors) +
  labs(x = "Study Duration (days)", y = "Study Count", fill = "Continent") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0))

## Assemble Fig2 ----
shared_legend <- get_legend(
  timeline_plot + theme(legend.position = "right")
)

hist_top_clean <- hist_top + guides(fill = "none") + theme(legend.position = "none")
hist_bottom_clean <- hist_bottom + guides(fill = "none") + theme(legend.position = "none")
timeline_clean <- timeline_plot + guides(color = "none") + theme(legend.position = "none")

fig2a_labeled <- fig2a + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
hist_top_labeled <- hist_top_clean + theme(plot.margin = unit(c(0.5, 0.5, 0.02, 0.5), "lines"))
hist_bottom_labeled <- hist_bottom_clean + theme(plot.margin = unit(c(0.02, 0.5, 0.5, 0.5), "lines"))
timeline_labeled <- timeline_clean + theme(plot.margin = unit(c(0.02, 0.5, 0.02, 0.5), "lines"))

right_column <- plot_grid(
  hist_top_labeled, 
  timeline_labeled, 
  hist_bottom_labeled,
  ncol = 1,
  rel_heights = c(0.5, 3, 0.5),
  labels = c("b", "c", "d"),
  label_size = 14,
  label_fontface = "bold",
  align = "v",
  axis = "lr"    
)

fig2_main <- plot_grid(
  fig2a_labeled,
  right_column,
  ncol = 2,
  rel_heights = c(0.8, 1),
  rel_widths = c(1, 1),
  labels = c("a", ""),
  label_size = 14,
  label_fontface = "bold",
  align = "hv",      
  axis = "tb" 
)

fig2_combined <- plot_grid(
  fig2_main,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.08)
)

ggsave("Figures/fig2.png", fig2_combined,
       width = 16, height = 12, dpi = 300, bg = "white")
