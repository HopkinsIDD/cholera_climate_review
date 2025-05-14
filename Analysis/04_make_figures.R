# This script makes figures


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)
library(dplyr)
library(RColorBrewer)

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
  )) %>% 
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
  add_continent(country_col="eff_country") %>% 
  filter_final_analysis_dataset(., criteria = expr(start_year > 1961))

### Identify cholera endemic countries ----
## sql code to extract WHO annual report data from Cholera Taxonomy Database
#select * from observations where 1=2;
#\copy (select time_left,time_right,suspected_cases,deaths,location_id, "primary" from observations where (observation_collection_id <=480) and (observation_collection_id>=300) ) to 'who_annual_report.csv' (format csv, delimiter ';');

who_data<-read.csv("Data/endemicity/who_annual_report.csv") %>%
  tidyr::separate(obs,c('time_left','time_right','suspected_cases','deaths','location_id','primary'),sep = ';',remove=T) %>%
  filter(primary == 't') %>% 
  merge(.,read.csv("Data/endemicity/country_location_id_name.csv"),by.x="location_id",by.y="id",all.x=TRUE) %>%
  filter(is.na(region)==F) %>%
  arrange(region) %>%
  mutate(country = as.factor(region),suspected_cases=as.numeric(suspected_cases)) %>% 
  filter(time_right >= as.Date("2000-01-01") & time_right <= as.Date("2023-12-31")) %>% 
  mutate(year = lubridate::year(time_right)) %>% 
  select(country,suspected_cases,year)


fig2a <- plot_spatial_distribution_map(plot_data = plot_data,
                                       country_col="eff_country",
                                       who_data = who_data,
                                       group_col = "admin_level",
                                       group_levels = c("regional","national","sub-national"))

## Fig 2b Temporal ----
continent_colors <- c(
  "Africa" = "#FFA378",
  "Asia" = "#A8B545",
  "Other" = "#8C796D"
)

# Make the plot
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

# Set up decade ranges and colors
decade_df <- data.frame(
  xmin = seq(as.Date("1971-01-01"), as.Date("2011-01-01"), by = "10 years"),
  xmax = seq(as.Date("1981-01-01"), as.Date("2021-01-01"), by = "10 years"),
  decade = paste0(seq(1970, 2010, by = 10), "s"),
  fill_color = c("mistyrose", "lightblue", "honeydew", "lavender", "wheat")
)#, "lightcyan"

custom_breaks <- seq.Date(from = as.Date("1970-01-01"), to = max(plot_data$year_tr, na.rm = TRUE), by = "5 years")
custom_labels <- c(format(custom_breaks, "%Y"))

hist_top <- ggplot(
  plot_data %>%
    distinct(studykey, year_tl, year_tr, continent) %>%
    mutate(start_year = lubridate::year(year_tl),
           end_year = lubridate::year(year_tr)) %>%
    mutate(start_year = ifelse(start_year < 1961, 1961, start_year)) %>%
    mutate(end_year = ifelse(start_year <= 1961, 1961, end_year)) %>%
    uncount(end_year - start_year + 1, .id = "year") %>%
    mutate(year = start_year + as.numeric(year) - 1), 
  aes(x = as.Date(paste0(year, "-01-01")), fill = continent, color = continent)
) +
  # Add colored rectangles first (background)
  geom_rect(
    data = decade_df,
    aes(xmin = xmin, xmax = xmax, ymin = 40, ymax = Inf, fill = decade),
    inherit.aes = FALSE,
    alpha = 1  # Transparent so histogram still visible
  ) +
  # Main histogram
  geom_histogram(binwidth = 365, position = "stack") +
  # Your continent colors
  scale_fill_manual(
    values = c(continent_colors, setNames(decade_df$fill_color, decade_df$decade))
  ) +
  scale_color_manual(values = continent_colors) +
  scale_x_date(
    breaks = custom_breaks,
    labels = custom_labels,
    expand = expansion(mult = c(0.01, 0.01)),
    sec.axis = dup_axis(
      breaks = seq(as.Date("1961-01-01"), as.Date("2021-01-01"), by = "10 years"),
      labels = paste0(seq(1960, 2020, by = 10), "s")
    )
  ) +
  labs(x = NULL, y = "Study count", fill = "Continent") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x.top = element_text(angle = 0, vjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0)
  )

## Assemble Fig2 ----
shared_legend <- get_legend(
  timeline_plot + theme(legend.position = "right")
)

# Remove shared legend
hist_top_clean <- hist_top + guides(fill = "none") + theme(legend.position = "none")
timeline_clean <- timeline_plot + guides(color = "none") + theme(legend.position = "none")

fig2a_labeled <- fig2a + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
hist_top_labeled <- hist_top_clean + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
timeline_labeled <- timeline_clean + theme(plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "lines"))

# Create left column of the figure
left_column <- plot_grid(
  fig2a_labeled,
  hist_top_labeled, 
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("a", "c"),
  label_size = 14,
  label_fontface = "bold",
  align = "v",
  axis = "lr"    
)

fig2_main <- plot_grid(
  left_column, 
  timeline_labeled,
  ncol = 2,
  rel_heights = c(1,1),
  rel_widths = c(1, 1),
  labels = c("", "b"),
  label_size = 14,
  label_fontface = "bold",
  align = "hv",      
  axis = "tb" 
)

# Add shared legend
fig2_combined <- plot_grid(
  fig2_main,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.08)
)

ggsave("Figures/fig2.png", fig2_combined,
       width = 18, height = 10, dpi = 300, bg = "white")

# Supplement figures ----
## Figure S1.Map of all location shapefiles ----
library(tidyverse)
library(rgeoboundaries)
library(sf)
sf::sf_use_s2(FALSE)
world <- gb_adm0()  # ADM0 = national boundaries
shapefiles <- readRDS("Data/shapefiles/location_shapefiles.rds")
shapefiles <- st_transform(shapefiles, crs = st_crs(world))

loc_subnational <- plot_data %>% 
  distinct(lctn_pr,admin_level) %>% 
  filter(admin_level == 'sub-national')

loc_national <- plot_data %>% 
  distinct(lctn_pr,admin_level) %>% 
  filter(admin_level == 'national')

loc_regional <- plot_data %>% 
  distinct(lctn_pr,admin_level) %>% 
  filter(admin_level == 'regional')

blues <- brewer.pal(n = 4, name = "Blues")[-1]  # Remove lightest
names(blues) <- c("Regional", "National", "Subnational")  # assign manually to match levels

shapefiles_tagged <- bind_rows(
  shapefiles %>% filter(lctn_pr %in% unique(loc_regional$lctn_pr)) %>% mutate(spatial_scale = "Regional"),
  shapefiles %>% filter(lctn_pr %in% unique(loc_national$lctn_pr)) %>% mutate(spatial_scale = "National"),
  shapefiles %>% filter(lctn_pr %in% unique(loc_subnational$lctn_pr)) %>% mutate(spatial_scale = "Subnational")
)

shapefiles_tagged$spatial_scale <- factor(
  shapefiles_tagged$spatial_scale,
  levels = c("Regional", "National", "Subnational")
)

figS1 <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", size = 0.2) +
  geom_sf(data = shapefiles_tagged, aes(fill = spatial_scale, color = spatial_scale), size = 0.4) +
  scale_fill_manual(
    name = "Spatial Scale",
    values = blues
  ) +
  scale_color_manual(
    values = blues,
    guide = "none"
  ) +
  theme_bw() +
  coord_sf(expand = FALSE)

ggsave('Figures/figS1.png',plot=figS1,width = 10, height = 8, dpi = 300, bg = "white")

## Figure S2.Spatial distribution of extracted covariates ----

## Figure S3. Histogram of time period duration and median by continent ----
plot_data_binned <- plot_data %>%
  mutate(duration_years = floor(duration_days / 365),
         bin_group = floor(duration_years / 5) * 5,
         bin_midpoint = bin_group + 2.5)  # for bar center on numeric x

plot_data_summary <- plot_data_binned %>%
  group_by(continent, bin_midpoint) %>%
  summarise(n = n(), .groups = "drop")

figS3 <- ggplot() +
  geom_col(data = plot_data_summary,
           aes(x = bin_midpoint, y = n, fill = continent, color = continent),
           position = position_dodge(width = 4.5), width = 4) +  # spacing control
  geom_vline(data = duration_medians, 
             aes(xintercept = median_duration / 365, color = continent, linetype = line_label),
             linewidth = 0.8) +
  scale_fill_manual(values = continent_colors) +
  scale_color_manual(values = continent_colors) +
  scale_linetype_manual(values = c("Median duration" = "dashed")) +
  scale_x_continuous(
    name = "Study Duration (years, grouped every 5 years)",
    breaks = seq(0, max(plot_data$duration_days, na.rm = TRUE)/365, by = 5)
  ) +
  labs(
    y = "Study Count",
    fill = "Continent",
    color = "Continent",
    linetype = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0)
  )

ggsave('Figures/figS3.png',plot=figS3,width = 12, height = 6, dpi = 300, bg = "white")

## Figure S4. Histogram of time period duration and median by country ----
plot_data_country <- plot_data %>%
  separate_rows(eff_country, sep = "\\|") %>%
  filter(!is.na(duration_days), duration_days > 0)

duration_medians <- plot_data_country %>%
  group_by(eff_country) %>%
  summarise(median_duration = median(duration_days, na.rm = TRUE)) %>%
  mutate(line_label = "Median duration")

figS4 <- ggplot(plot_data_country) +
  geom_histogram(aes(x = duration_days), 
                 binwidth = 365, position = "stack") +
  geom_vline(data = duration_medians, 
             aes(xintercept = median_duration, linetype = line_label),
             linewidth = 0.4, color = "black") +
  scale_linetype_manual(values = c("Median duration" = "dashed")) +
  labs(
    x = "Study Duration (days)",
    y = "Study Count",
    fill = "Continent",
    color = "Continent",
    linetype = "",
    title = "Figure S4. Histogram of study duration and median by country"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0),
    strip.text = element_text(size = 8)
  ) +
  facet_wrap(~eff_country)

ggsave('Figures/figS4.png',plot=figS4,width = 16, height = 8, dpi = 300, bg = "white")

# Exploratory figures ----
## Fig.E1: distributions of % of years with reported cases between 2000-2023 ----
summary <- identify_endemicity(
  processed_data = plot_data,
  who_data = who_data,
  group_col = "exp_class"
)

ggplot(summary, aes(x = percent_year_with_cases)) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "white", boundary = 0, closed = "left") +
  facet_wrap(~ exp_class) +
  labs(
    title = "Distribution of % of Years with Reported Cholera Cases (2000-2023)",
    x = "% of Years with Reported Cases",
    y = "Number of Countries"
  ) +
  theme_bw()
