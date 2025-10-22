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
library(here)
library(gt)

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
  mutate(score_cat = cut(tot_score, c(0, 5, 9, Inf), labels = c("low", "medium", "high"), include.lower = T, right = F),
         score_cat = forcats::fct_relevel(score_cat, rev( c("low", "medium", "high"))))


# Compute statistics for direction of effect plot
dir_data <- est_data %>% 
  group_by(studykey, out_class_case_presentation, out_class_metric, exp_class, continent, admin_level) %>%
  summarise_directions() %>% 
  ungroup()  %>% 
  mutate(studykey = str_replace(studykey, "\xe6", "e"),
         dir_icon = case_when(overall_dir == "higher" ~ "caret-up",
                              overall_dir == "lower" ~ "caret-down",
                              TRUE ~ "caret-left,caret-right"),
         overall_dir = factor(overall_dir, labels = c("Positive", "Negative", "Conflicting")), 
         pubyear = str_extract(studykey, "[0-9]{4}") %>% as.numeric()) 


# Figure 2: Spatial-temporal coverage of studies --------------------------

## Fig 2a Spatial coverage map ----
plot_data <- readRDS("generated_data/processed_estimate_data.rds") %>% 
  rowwise() %>% 
  mutate(
    year_tl = ymd(eff_tl),
    year_tr = ymd(eff_tr),
    duration_days = as.numeric(year_tr - year_tl),
    start_year = year(year_tl),
    start_date = lubridate::make_date(start_year, 1, 1)
  ) %>% 
  filter_final_analysis_dataset(., criteria = expr(start_year > 1961))  %>% 
  mutate(admin_level_simple = if_else(admin_level == "regional", "regional", "national or sub-national"))


fig2a <- plot_spatial_distribution_map(
  plot_data = plot_data,
  country_col = "eff_country",
  who_data = get_WHO_case_data(),
  group_col = "admin_level_simple",
)

## Fig 2b Temporal ----

# Make the plot
timeline_plot <- ggplot(plot_data %>% arrange(year_tl)) +
  geom_segment(aes(x = year_tl, xend = year_tr,
                   y = reorder(studykey, year_tl),
                   yend = reorder(studykey, year_tl),
                   color = continent),
               linewidth = 2) +
  scale_color_manual(values = colors_continent()) +
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
    values = c(colors_continent(), setNames(decade_df$fill_color, decade_df$decade))
  ) +
  scale_color_manual(values = colors_continent()) +
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

# ggsave(filename = "figures/fig2.png", fig2_combined,
#        width = 18, height = 10, dpi = 300, bg = "white")


## Second pass figure 2 ----
timeline_plot2 <- timeline_plot + 
  geom_segment(aes(x = year_tl, xend = year_tr,
                   y = reorder(studykey, year_tl),
                   yend = reorder(studykey, year_tl)),
               color = "white",
               linewidth = .3) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  scale_x_date(breaks = seq.Date("1970-01-02", "2020-01-02", by = "10 years"), 
               date_labels = "%Y", limits = as.Date(c("1965-01-01", "2020-12-31"))) +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("Temporal distribution of studies") +
  guides(color = "none")


hist_top2 <- ggplot(
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
  # Main histogram
  geom_histogram(binwidth = 365, position = "stack") +
  # Your continent colors
  scale_fill_manual(
    values = c(colors_continent())
  ) +
  scale_color_manual(values = colors_continent()) +
  labs(x = NULL, y = "Study count", fill = "Continent", color = "Continent") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.x.top = element_text(angle = 0, vjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0)
  ) +
  scale_x_date(breaks = seq.Date("1970-01-02", "2020-01-02", by = "5 years"), 
               date_labels = "%Y", limits = as.Date(c("1965-01-01", "2020-12-31"))) +
  theme(panel.grid.minor = element_blank()) + 
  coord_flip()  +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 0))

p_fig2_lower <- cowplot::plot_grid(
  timeline_plot2,
  hist_top2,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(1, .1)
)


fig2_v2 <- (fig2a + timeline_plot2 + hist_top2) + 
  plot_layout(
    design = "
    AAAAAAA
    AAAAAAA
    AAAAAAA
    AAAAAAA
    BBBBBBC
    BBBBBBC
    ")

ggsave(fig2_v2, filename = "figures/figure_2_v2.png",
       width = 16.5, height = 11, dpi = 300)

# Figure 3: direction of effect -------------------------------------------

# Make overall plot
compose_effdir_plot(data = dir_data,
                    filename = "figures/figure_3.png")


# Figure 4: meta-regression -----------------------------------------------

# Load 
res_files <- dir("generated_data/", pattern = "stan_genquant", full.names = T) %>%
  # Get new runs with full set of covariates
  str_subset("strat") %>%
  str_subset("direction")

# Get effect sizes
es_effect_dir <- get_effect_sizes(res_files) %>%
  mutate(covar_type = if_else(str_detect(what, "resolution|lag|admin"), 
                              "analytical approach", "contextual factors") %>% 
           factor(levels = c("contextual factors", "analytical approach")))


# Save
saveRDS(es_effect_dir, file = "generated_data/metareg_effect_sizes.rds")

pd <- position_dodge(width = .5)

fig4 <- es_effect_dir %>% 
  filter(stratification == "none") %>% 
  ggplot(aes(x = strata, y = mean)) +
  geom_hline(aes(yintercept = 1), lwd = .3, lty = 2) +
  geom_point(aes(pch = out_class, color = out_class), position = pd) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5, color = out_class, lty = out_class), 
                position = pd, width = 0) +
  ggh4x::facet_nested(exp_class ~ covar_type + what, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggthemes::scale_color_few("Dark") +
  scale_linetype_manual(values = c(1, 6)) +
  labs(
    x = "covariate levels",
    y = "odds ratio of positive direction of association", 
    shape = "Outcome type",
    lty = "Outcome type",
    color = "Outcome type"
  ) +
  theme(
    strip.text.x = element_text(size = 6)
  )

ggsave(fig4, filename = "figures/figure_4.png",
       width = 12.5, height = 7, dpi = 300)

# Supplement figures ------------------------------------------------------

## Figure S1.Map of all location shapefiles ----
library(tidyverse)
library(rgeoboundaries)
library(sf)
sf::sf_use_s2(FALSE)
world <- gb_adm0()  # ADM0 = national boundaries

# fix the regional shapefiles in the location_shapefiles.rds
# data <- readRDS('Data/shapefiles/location_shapefiles.rds')
# 
# lp = "new-45"
# ctry_list = "AGO|BDI|BEN|BFA|BWA|CAF|CIV|CMR|COD|COG|DJI|ERI|ETH|GAB|GHA|GIN|GMB|GNB|GNQ|KEN|LBR|LSO|MLI|MOZ|MRT|MWI|NAM|NER|NGA|RWA|SEN|SLE|SOM|SWZ|TCD|TGO|TZA|UGA|ZAF|ZMB|ZWE"
# 
# country_vec <- str_split(ctry_list, pattern = "\\|", simplify = TRUE) %>% as.vector()
# tmp <- st_union(world %>% filter(shapeGroup %in% country_vec))
# plot(tmp)
# shp <- data %>% filter(lctn_pr == lp)
# plot(shp$geometry)
# 
# data[data$lctn_pr == lp,]$geometry = tmp
# 
# plot(data[data$lctn_pr == lp,]$geometry)
# 
# #data <- data %>% filter(lctn_pr != lp)
# 
# saveRDS(data,"Data/shapefiles/location_shapefiles.rds")

###
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

library(colorspace)
colors <- colorspace::sequential_hcl(3, palette = "BluYl")
names(colors) <- c( "Subnational","National", "Regional")

shapefiles_tagged <- bind_rows(
  #shapefiles %>% filter(lctn_pr %in% unique(loc_regional$lctn_pr)) %>% mutate(spatial_scale = "Regional"),
  shapefiles %>% filter(lctn_pr %in% unique(loc_national$lctn_pr)) %>% mutate(spatial_scale = "National"),
  shapefiles %>% filter(lctn_pr %in% unique(loc_subnational$lctn_pr)) %>% mutate(spatial_scale = "Subnational")
)

shapefiles_tagged$spatial_scale <- factor(
  shapefiles_tagged$spatial_scale,
  levels = c("Regional", "National", "Subnational")
)
# if a region contains both national and subnational level data, show the subnational data.
bbox <- sf::st_bbox(shapefiles_tagged)

figS1 <- ggplot(data = shapefiles_tagged) +
  geom_sf(data = world, fill = "grey95", color = "grey70", size = 0.1) +
  geom_sf(data = shapefiles_tagged, aes(fill = spatial_scale),color = "grey70", size = 0.1) +
  scale_fill_manual(
    name = "Finest Spatial\nScale Reported",
    values = colors
  ) +
  scale_color_manual(
    values = colors,
    guide = "none"
  ) +
  theme_bw() +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  )+
  facet_grid(~spatial_scale)

ggsave('figures/figS1.png',plot=figS1,width = 10, height = 4, dpi = 300, bg = "white")

## Figure S2.Spatial distribution of extracted covariates ----
shapefiles <- readRDS("Data/shapefiles/location_shapefiles.rds")
shapefiles <- st_transform(shapefiles, crs = st_crs(world))

external_sf <- plot_data %>% 
  filter(admin_level %in% c("national", "sub-national")) %>% 
  distinct(average_hdi, new_climatic_zone, inland_coast, mean_percent_year_with_cases,lctn_pr, admin_level) %>%
  left_join(
    shapefiles,
    by = "lctn_pr"
  ) %>%
  mutate(
    hdi_cat = case_when(
      average_hdi >= 0     & average_hdi < 0.55 ~ "0-0.55",
      average_hdi >= 0.55  & average_hdi < 0.7  ~ "0.55-0.7",
      average_hdi >= 0.7   & average_hdi < 0.8  ~ "0.7-0.8",
      average_hdi >= 0.8   & average_hdi <= 1   ~ "0.8-1.0",
      TRUE ~ NA_character_
    ),
    hdi_cat = factor(hdi_cat, levels = c("0-0.55", "0.55-0.7", "0.7-0.8", "0.8-1.0"))
  ) %>% 
  mutate(endemicity_group = case_when(
    mean_percent_year_with_cases <= 50 ~ "0-50%",
    mean_percent_year_with_cases > 50 & mean_percent_year_with_cases <= 95 ~ "50-95%",
    mean_percent_year_with_cases > 95 ~ ">95%",
  ),
  endemicity_group = factor(endemicity_group, levels = c("0-50%","50-95%",">95%"))
  ) %>% 
  sf::st_as_sf()

bbox <- sf::st_bbox(external_sf)

map_hdi <- ggplot(data = external_sf) +
  geom_sf(data = world, fill = "grey95", color = "grey90", size = 0.1) +
  geom_sf(data = external_sf, aes(fill = hdi_cat), color = "grey60", size = 0.1) +
  scale_fill_brewer(
    palette = "Blues",
    name = "HDI Category",
    na.value = "white"
  ) +
  theme_bw() +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  ) +
  facet_grid(~ admin_level)

map_climate <- ggplot(data = external_sf) +
  geom_sf(data = world, fill = "grey95", color = "grey90", size = 0.1) +
  geom_sf(data = external_sf, aes(fill = new_climatic_zone), color = "grey60", size = 0.1) +
  scale_fill_brewer(
    palette = "YlGnBu",
    name = "Climatic Zone",
    na.value = "white"
  ) +
  theme_bw() +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  )+
  facet_grid(~ admin_level)

map_inland <- ggplot(data = external_sf) +
  geom_sf(data = world, fill = "grey95", color = "grey90", size = 0.1) +
  geom_sf(data = external_sf, aes(fill = inland_coast), color = "grey60", size = 0.1) +
  scale_fill_manual(
    values = c("coastal" = "#6baed6", "inland" = "#08519c"),
    name = "Proximity to ocean",
    na.value = "white"
  ) +
  theme_bw() +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  )+
  facet_grid(~ admin_level)


map_endemicity <- ggplot(data = external_sf) +
  geom_sf(data = world, fill = "grey95", color = "grey90", size = 0.1) +
  geom_sf(data = external_sf, aes(fill = endemicity_group), color = "grey60", size = 0.1) +
  scale_fill_brewer(
    palette = "YlGnBu",
    name = "% of years with\nWHO reported cases\nbewteen 2000-2023",
    na.value = "white"
  ) +
  theme_bw() +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  )+
  facet_grid(~ admin_level)


figS2 <- map_hdi / map_climate / map_inland / map_endemicity +
  plot_layout(ncol = 1)

ggsave('figures/figS2.png',plot=figS2,width = 12, height = 12, dpi = 300, bg = "white")

## Figure S3. Histogram of time period duration and median by continent ----
plot_data_binned <- plot_data %>%
  mutate(duration_years = floor(duration_days / 365),
         bin_group = floor(duration_years / 5) * 5,
         bin_midpoint = bin_group + 2.5)  # for bar center on numeric x

plot_data_summary <- plot_data_binned %>%
  group_by(continent, bin_midpoint) %>%
  summarise(n = n(), .groups = "drop")

plot_data_country <- plot_data %>% 
  group_by(studykey,eff_country,continent) %>%
  summarise(duration_days = mean(duration_days)) %>% # average study duration at the study level
  separate_rows(eff_country, sep = "\\|") %>%
  filter(!is.na(duration_days), duration_days > 0)

duration_medians_continent <- plot_data_country %>%
  group_by(continent) %>% 
  summarise(median_duration = median(duration_days, na.rm = TRUE)) %>%
  mutate(line_label = "Median duration")

duration_medians <- plot_data_country %>%
  group_by(eff_country) %>% 
  summarise(median_duration = median(duration_days, na.rm = TRUE)) %>%
  mutate(line_label = "Median duration")


figS3 <- ggplot() +
  geom_col(data = plot_data_summary,
           aes(x = bin_midpoint, y = n, fill = continent, color = continent),
           position = position_dodge(width = 4.5), width = 4) +  # spacing control
  geom_vline(data = duration_medians_continent, 
             aes(xintercept = median_duration / 365, color = continent, linetype = line_label),
             linewidth = 0.8) +
  scale_fill_manual(values = colors_continent()) +
  scale_color_manual(values = colors_continent()) +
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
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10), 
    plot.title = element_text(face = "bold", hjust = 0),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

ggsave('figures/figS3.png',plot=figS3,width = 12, height = 6, dpi = 300, bg = "white")

## Figure S4. Histogram of time period duration and median by country ----

figS4 <- ggplot(plot_data_country) +
  geom_histogram(aes(x = duration_days/365,fill = continent), 
                 binwidth = 1, position = "stack") +
  geom_vline(data = duration_medians, 
             aes(xintercept = median_duration/365, linetype = line_label),
             linewidth = 0.4, color = "black") +
  scale_linetype_manual(values = c("Median duration" = "dashed")) +
  scale_fill_manual(values = colors_continent()) +
  scale_x_continuous(
    name = "Study Duration (years)",
    breaks = seq(0, max(plot_data$duration_days, na.rm = TRUE)/365, by = 5)
  ) +
  labs(
    y = "Study Count",
    fill = "Continent",
    color = "Continent",
    linetype = "",
    title = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0),
    strip.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  facet_wrap(~eff_country)

ggsave('figures/figS4.png',plot=figS4,width = 16, height = 16, dpi = 300, bg = "white")

## Figure S5. Tile plot by study for each exposure variable ----
quality_data <- plot_data %>%
  make_scores() %>% 
  group_by(studykey,exp_class) %>%
  summarise(
    exp_temp_perfect = as.integer(all(exp_temp_tf == 1)),
    exp_spatial_perfect = as.integer(all(exp_spatial_tf == 1)),
    exp_source_perfect = as.integer(all(exp_source_tf == 1)),
    
    out_class_case_presentation_perfect = as.integer(all(out_class_case_presentation_tf == 1)),
    out_temp_perfect = as.integer(all(out_temp_tf == 1)),
    out_spatial_perfect = as.integer(all(out_spatial_tf == 1)),
    out_case_df_perfect = as.integer(all(out_case_df_tf == 1)),
    
    a_class_perfect = as.integer(all(a_class_tf == 1)),
    a_uncert_def_perfect = as.integer(all(a_uncert_def_tf == 1)),
    
    .groups = "drop"
  ) %>% 
  select(studykey, 
         exp_class,
         exp_temp_perfect,
         exp_spatial_perfect,
         exp_source_perfect,
         out_class_case_presentation_perfect,
         out_temp_perfect,
         out_spatial_perfect,
         out_case_df_perfect,
         a_class_perfect,
         a_uncert_def_perfect)

plot_long <- quality_data %>%
  pivot_longer(
    cols = -c(studykey, exp_class),
    names_to = "variable",
    values_to = "is_perfect"
  ) %>%
  mutate(
    # Grouping type for facetting
    type = case_when(
      str_starts(variable, "exp_") ~ "Exposure",
      str_starts(variable, "out_") ~ "Outcome",
      str_starts(variable, "a_") ~ "Analysis"
    ),
    
    # Descriptive label without prefixes
    variable_label = recode(variable,
                            "exp_source_perfect" = "Q1",
                            "exp_spatial_perfect" = "Q2",
                            "exp_temp_perfect" = "Q3",
                            "out_case_df_perfect" = "Q4",
                            "out_class_case_presentation_perfect" = "Q5",
                            "out_spatial_perfect" = "Q6",
                            "out_temp_perfect" = "Q7",
                            "a_class_perfect" = "Q8",
                            "a_uncert_def_perfect" = "Q9"
    )
  ) %>%
  mutate(
    exp_class = factor(exp_class,levels = c("precipitation","temperature","flood","drought")),
    type = factor(type, levels = c("Exposure","Outcome","Analysis")),
    studykey = fct_rev(factor(studykey, levels = sort(unique(studykey))))
  ) %>% 
  mutate(
    variable_label = factor(variable_label),
    studykey = factor(studykey, levels = sort(unique(studykey)))
  )

figS5_prep <- ggplot(plot_long %>% filter(exp_class == "precipitation"), aes(x = variable_label, y = studykey, fill = factor(is_perfect))) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("0" = "#1f78b4", "1" = "#ffdf00"),  # Blue (0), Yellow (1)
    name = "Reporting status"
  ) +
  labs(x = "Reporting Completeness Questions", y = "Study") +
  facet_grid(exp_class ~ type, scales = "free", space = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

figS5_other3 <- ggplot(plot_long %>% filter(exp_class != "precipitation"), aes(x = variable_label, y = studykey, fill = factor(is_perfect))) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("0" = "#1f78b4", "1" = "#ffdf00"),  # Blue (0), Yellow (1)
    name = "Reporting status"
  ) +
  labs(x = "Reporting Completeness Questions", y = "Study") +
  facet_grid(exp_class ~ type, scales = "free", space = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )


figS5 <- (figS5_prep / figS5_other3) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave('figures/figS5_1022.png',plot=figS5,width = 10, height = 12, dpi = 300, bg = "white")


## Figure S6-S8: effect size stratification by continent, climate and endemicity ----

walk(
  c("new_climatic_zone", "continent", "pct_years_cholera_cat"),
  function(x) {
    pd <- position_dodge(width = .4)
    
    p_strat <- es_effect_dir %>% 
      filter(stratification == x) %>% 
      ggplot(aes(x = strata, y = mean)) +
      geom_hline(aes(yintercept = 1), lwd = .3, lty = 3) +
      geom_errorbar(aes(ymin = q5, ymax = q95, color = level, lty = out_class), 
                    position = pd, width = 0, alpha = .8, lwd = .3) +
      geom_point(aes(fill = level), pch = 21, color = "black", position = pd, stroke = .2) +
      ggh4x::facet_nested(out_class + exp_class ~ covar_type + what, scales = "free_x", space = "free_x") +
      theme_bw() +
      scale_linetype_manual(values = c(1, 6)) +
      scale_color_manual(values = get_stratification_colors(stratification = x)) +
      scale_fill_manual(values = get_stratification_colors(stratification = x)) +
      labs(y = "odds ratio", color = get_covar_dict()[x]) +
      guides(shape = "none") +
      scale_y_continuous(limits = c(0, 4))+
      theme(
        strip.text.x = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )
    
    ggsave(p_strat, filename = str_glue("figures/supfig_effect_sizes_stratified_{x}.png"),
           width = 12.5, height = 8)
    
  })

## Effect direction stratification ----

walk(c("Asia", "Africa", "Other"), function(x) {
  compose_effdir_plot(data = dir_data %>% filter(continent == x),
                      filename = str_glue("figures/supfig_eff_dir_plot_{x}.png"))
})


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

## ICC for all fits ----

# Load 
fit_files <- dir("generated_data/", pattern = "stan_genquant", full.names = T) %>%
  # Get new runs with full set of covariates
  str_subset("strat") %>%
  str_subset("direction") %>% 
  str_subset("none")



get_i_squared(fit_files) %>% 
  ggplot(aes(x = exp_class, y = mean)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0) +
  geom_point() + 
  facet_grid(. ~ out_class) +
  theme_bw()

# Scraps ------------------------------------------------------------------


#  Stratification by continent 
# 
# supfig_strat_continent <- es_effect_dir %>% 
#   filter(stratification == "continent") %>% 
#   ggplot(aes(x = strata, y = mean)) +
#   geom_hline(aes(yintercept = 0), lwd = .3, lty = 2) +
#   geom_point(aes(color = level, pch = out_class), position = pd) +
#   geom_errorbar(aes(ymin = q5, ymax = q95, color = level, lty = out_class), 
#                 position = pd, width = 0) +
#   ggh4x::facet_nested(exp_class + out_class ~ covar_type + what, scales = "free_x", space = "free_x") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   scale_color_manual(values = colors_continent()) +
#   scale_linetype_manual(values = c(1, 6)) +
#   labs(y = "effect size")
# 
# ggsave(supfig_strat_continent, filename = "figures/supfig_effect_sizes_continent.png",
#        width = 12, height = 8, dpi = 300)
# 
#  Stratification by climate zone
# 
# supfig_strat_climate <- es_effect_dir %>% 
#   filter(stratification == "new_climatic_zone") %>% 
#   ggplot(aes(x = strata, y = mean)) +
#   geom_hline(aes(yintercept = 0), lwd = .3, lty = 2) +
#   geom_point(aes(color = level, pch = out_class), position = pd) +
#   geom_errorbar(aes(ymin = q5, ymax = q95, color = level, lty = out_class), 
#                 position = pd, width = 0) +
#   ggh4x::facet_nested(exp_class + out_class ~ covar_type + what, scales = "free_x", space = "free_x") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   ggthemes::scale_color_few("Dark") +
#   scale_linetype_manual(values = c(1, 6)) +
#   labs(y = "effect size")
# 
# ggsave(supfig_strat_climate, filename = "figures/supfig_effect_sizes_climate.png",
#        width = 12, height = 8, dpi = 300)
