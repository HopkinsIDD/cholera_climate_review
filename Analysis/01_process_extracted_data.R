# This script is processes the extracted systematic review data

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(here)
library(optparse)
library(patchwork)
library(lubridate)
library(purrr)

source("Analysis/utils.R")

# User-supplied options
option_list <- list(
  make_option(c("-p", "--do_plots"), 
              default = FALSE, action ="store", type = "logical", help = "Make intermediate plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))

print(opt)

# Load data ---------------------------------------------------------------

dat_path <- here::here("Data/extracted_data/cholera_climate_systematic_review_extract_into_sheets.xlsx")

# Load effect extraction data
est_data <- #readxl::read_excel(dat_path, sheet = "estimates_0329") %>%
  read_csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv") %>%
  # Fix dates
  mutate(across(c("eff_tl", "eff_tr"), ~ char_to_date(.) %>% as.Date())) %>% 
  # Filter out estimates before 7th pandemic 
  filter_final_analysis_dataset(
    ., criteria = list(expr(year(eff_tl) > 1961))
  ) %>% 
  # Unpack countries for multi-country studies
  mutate(country_list = map(eff_country, ~ str_split(., pattern = "\\|")[[1]])) %>% 
  # Fix case definition 
  mutate(out_class_case_presentation = case_when(
    out_class_case_presentation == "bacteriologically confirmed cases only" ~ "Confirmed",
    out_class_case_presentation == "suspected cases only" ~ "Suspected",
    out_class_case_presentation == "deaths" ~ "Deaths",
    # out_class_case_presentation == "sero-incidence infections and bacteriologically confirmed cases" ~ "Sero-incidence and Confirmed cases",
    # out_class_case_presentation == "sero-incidence infections, bacteriologically confirmed cases and suspected cases" ~ "Sero-incidence and Confirmed and Suspected cases",
    out_class_case_presentation == "not stated" ~ "Unknown",
    TRUE ~ "Mixed"
  )) %>% 
  # Add new categorical variable for duration 
  mutate(duration = as.numeric(eff_tr - eff_tl)) %>%
  mutate(study_duration = case_when(
    duration <= 366 ~ "<1 year",
    duration <= 366*2 ~ "1-2 years",
    duration <= 366*5 ~ "2-5 years",
    duration <= 366*10 ~ "5-10 years",
    TRUE ~ ">10 years"
  )) %>% 
  # Fix estimates
  mutate(across(c("a_effect", "a_uncert_lb", "a_uncert_ub"), ~ as.numeric(.))) %>% 
  # Fix classes of estimates
  mutate(a_class_backup = a_class) %>% 
  mutate(a_class = case_when(str_detect(a_raw, "correlation") ~ "correlation",
                             a_class == "correlation" & !str_detect(a_class, "correlation") ~ "other",
                             TRUE ~ a_class),
         a_class = case_when(studykey == "Hulland2019" & exp_class == "precipitation" ~ "risk difference",
                             T ~ a_class))  %>%
  # Determine if effect is already on log scale (important for ratios)
  group_by(studykey, exp_class, a_raw) %>% 
  ungroup() %>% 
  # Se ratio to log10 ratio for plots
  mutate(
    a_effect_raw = a_effect,
    a_effect_natural_scale = case_when(
      # Update Feb 13 2025: We now have a column specifying the transformation to apply if any
      a_transform == "exp" ~ exp(a_effect),
      a_transform == "base10" ~ 10^a_effect,
      TRUE ~ a_effect
    ),
    a_uncert_lb_natural_scale = case_when(
      a_transform == "exp" ~ exp(a_uncert_lb),
      a_transform == "base10" ~ 10^a_uncert_lb,
      TRUE ~ a_uncert_lb
    ),
    a_uncert_ub_natural_scale = case_when(
      a_transform == "exp" ~ exp(a_uncert_ub),
      a_transform == "base10" ~ 10^a_uncert_ub,
      TRUE ~ a_uncert_ub
    ),
    # Now apply log10 to ratios
    a_effect = case_when(
      a_class == "ratio" ~ log10(a_effect_natural_scale),
      TRUE ~ a_effect_natural_scale
    ),
    # Now apply log10 to ratios
    a_uncert_lb = case_when(
      a_class == "ratio" ~ log10(a_uncert_lb_natural_scale),
      TRUE ~ a_uncert_lb_natural_scale
    ),
    # Now apply log10 to ratios
    a_uncert_ub = case_when(
      a_class == "ratio" ~ log10(a_uncert_ub_natural_scale),
      TRUE ~ a_uncert_ub_natural_scale
    ),
    # Update class to reflect change in ratios
    a_class = case_when(a_class == "ratio" ~ "log10_ratio",
                        TRUE ~ a_class)
  ) %>% 
  mutate(
    a_class = a_class_backup
  ) %>% 
  select(!a_class_backup) %>% 
  # Spatial scales
  mutate(
    admin_level = case_when(new_admin == "supranational" ~ "regional",
                            new_admin == "subnational" ~ "sub-national",
                            TRUE ~ new_admin),
    admin_level = factor(admin_level, 
                         levels = c("regional", "national", "sub-national")),
  ) %>% 
  # Add continent
  add_continent(country_col = "eff_country") %>% 
  # Modification of factor levels and labels
  mutate(
    # Levels of direction
    a_direction = factor(a_direction, levels = c("negative", "none", "positive", "not stated")),
    # Levels of exposure variable
    exp_class = factor(exp_class, 
                       levels = c("precipitation", "temperature", "flood", "drought")),
    # Define lags between exposures and outcomes
    lag_days = ifelse(lag_days == "missing", NA_character_, lag_days) %>% as.numeric(),
    lag_class = cut(lag_days, c(0, 1, 8, 15, 32, 62, 92, Inf),
                    include.lowest = TRUE, right = FALSE),
    lag_class = factor(lag_class, 
                       levels = c("[0,1)", "[1,8)", "[8,15)", "[15,32)", "[32,62)", "[62,92)", "[92,Inf)"),
                       labels = c("no lag", "1 week", "2 weeks", "1 month", "2 months", "3 months", "> 3 months")),
    # Define rainy season
    rainy_season = factor(rainy_season, 
                          levels = c("out_of_rainy_season", "mixed", "in_rainy_season"), 
                          labels = c("out of season", "both", "in season")),
    # Define urba/rural
    rural_urban = factor(rural_urban, 
                               levels = c("rural", "mixed", "urban")),
    urban_rural_class = rural_urban,
    continent = factor(continent),
    # Levels of uncertainty
    a_uncert_significant = factor(a_uncert_significant, 
                                  levels = c("TRUE", "FALSE", "NA"),
                                  labels = c("TRUE", "FALSE", "unspecified")),
    sig_report = case_when(
      a_uncert_def == "not stated" ~ FALSE,
      TRUE ~ TRUE
    ),
    # Levels of temporal scale
    out_temp_scale = factor(get_out_tempscale_dict()[as.character(out_temp_scale)], levels = get_out_tempscale_dict()),
    # Define categories of pct years with cholera
    pct_years_cholera_cat = cut(mean_percent_year_with_cases, c(0, 50, 95, Inf), include.lowest = T, right = FALSE) %>% 
      factor(levels = c("[0,50)", "[50,95)", "[95,Inf]"),
             labels = c("<50%", "50-95%", ">95%")),
    # Levels of climate zones
    new_climatic_zone = factor(new_climatic_zone, 
                               levels = c("Tropical", "Temperate", "Arid", "Polar"),
                               labels = c("tropical", "temperate", "arid", "polar")),
    # Levels of climate zones
    inland_coast = factor(inland_coast, 
                          levels = c("Coastal", "Inland"),
                          labels = c("coastal", "inland")),
    # Levels of HDI categories
    hdi_cat = factor(hdi_cat, 
                     levels = c("low hdi", "medium hdi", "high hdi", "very high hdi"),
                     labels = c("low", "medium", "high", "very high"))
  )

# Save
saveRDS(est_data, file = "generated_data/processed_estimate_data.rds")

# Metadata
meta_data <- readxl::read_excel(dat_path, sheet = "metadata_0329") %>% 
  filter(remove_keep == "keep")

# Timelines
est_timings <- est_data %>% 
  group_by(studykey, continent, eff_loc, admin_level) %>% 
  summarise(tl = min(eff_tl, na.rm = T),
            tr = max(eff_tr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(dt = as.numeric(difftime(tr, tl, units = "days")))  %>% 
  group_by(studykey) %>% 
  mutate(min_tl = min(tl)) %>% 
  ungroup() %>% 
  mutate(studykey = factor(studykey) %>% 
           forcats::fct_reorder(min_tl)) %>% 
  filter(tl > "1900-01-01")


# Save
saveRDS(est_timings, file = "generated_data/paper_timings.rds")


if (opt$do_plots) {
  
  # Number of papers by exposure class
  est_data %>% 
    group_by(studykey, exp_class) %>% 
    slice(1) %>% 
    ungroup() %>% 
    count(exp_class)
  
  # Plot of number of effect estimates by paper
  est_data %>% 
    count(studykey) %>% 
    ggplot() +
    geom_histogram(aes(x = n)) +
    theme_bw() +
    ggtitle("Number of effects per paper")
  
  # Plot of number of effect estimates by exposure type
  est_data %>% 
    count(studykey, exp_class) %>% 
    ggplot() +
    geom_histogram(aes(x = n, fill = exp_class), position = "dodge") +
    facet_grid(exp_class ~ .) +
    theme_bw() +
    ggtitle("Number of effects per paper and exposure class")
  
  
  # Map of countries [this takes time the first time it is run]
  ne_sf <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
    mutate(country = iso_a3)
  
  est_data_bbox <- sf::st_bbox(ne_sf %>% inner_join(est_data %>% distinct(country = eff_country)))
  
  # All isocodes for which to do extractions
  all_isocodes <- unique(ne_sf$iso_a3) %>% sort()
  
  # Count paper by admin level
  p_paper_admin_lev <- est_data %>%
    group_by(studykey, admin_level) %>% 
    slice(1)  %>% 
    count_countries(all_isocodes = all_isocodes,
                    group_col = c("admin_level")) %>% 
    inner_join(ne_sf, .) %>% 
    mutate(count = ifelse(count == 0, NA, count)) %>% 
    ggplot(aes(fill = count)) +
    geom_sf(color = "white", lwd = .05) +
    scale_fill_viridis_c(direction = -1, na.value = "lightgray") +
    facet_grid(admin_level ~ ., switch = "y") +
    theme_bw() +
    theme(plot.background = element_rect(fill = 'white')) +
    coord_sf(xlim = est_data_bbox[c("xmin", "xmax")], ylim = est_data_bbox[c("ymin", "ymax")])
  
  ggsave(p_paper_admin_lev, filename = "figures/map_papers_by_admin_level.png",
         width = 10, height = 7, dpi = 300)
  
  # Count paper by exp class and admin level
  p_paper_exp_class <- est_data %>%
    group_by(studykey, exp_class, admin_level) %>% 
    slice(1)  %>% 
    count_countries(all_isocodes = all_isocodes,
                    group_col = c("exp_class", "admin_level")) %>% 
    inner_join(ne_sf, ., by = c("iso_a3" = "country")) %>% 
    mutate(count = ifelse(count == 0, NA, count)) %>% 
    ggplot(aes(fill = count)) +
    geom_sf(color = "white", lwd = .05) +
    scale_fill_viridis_c(direction = -1, na.value = "lightgray") +
    facet_grid(admin_level ~ exp_class) +
    theme_bw()
  
  ggsave(p_paper_exp_class, filename = "figures/map_papers_by_exp_class.png",
         width = 10, height = 7, dpi = 300)
  
  # Count effects by country and admin level
  p_est_exp_class <- est_data %>%
    count_countries(all_isocodes = all_isocodes,
                    group_col = c("exp_class", "admin_level")) %>% 
    inner_join(ne_sf, ., by = c("iso_a3" = "country")) %>% 
    mutate(count = ifelse(count == 0, NA, count)) %>% 
    ggplot(aes(fill = count)) +
    geom_sf(color = "white", lwd = .05) +
    scale_fill_viridis_c() +
    facet_grid(admin_level ~ exp_class) +
    ggthemes::theme_map() +
    theme(plot.background = element_rect(fill = 'white'))
  
  ggsave(p_est_exp_class, filename = "figures/map_estimates_by_exp_class.png",
         width = 15, height = 5, dpi = 300)
}

# Compute scores ----------------------------------------------------------

# scores <- est_data[,colnames(est_data) %in% c("studykey","eff_country","exp_class","exp_temp_scale","exp_spatial_scale","exp_source","out_class_metric","out_class_case_presentation","out_medicattended","out_temp_scale","out_spatial_scale","case_def_raw","exp_out_spatial_scale_match","a_effect","a_adjusted","a_uncert_def")] %>% 
#   mutate(
#     exp_temp_tf = ifelse(exp_temp_scale == 99,0,1),
#     exp_spatial_tf = ifelse(exp_spatial_scale == "not stated",0,1),
#     exp_source_tf = ifelse(exp_source == "not stated",0,1),
#     out_class_case_presentation_tf = ifelse(out_class_case_presentation == "not stated",0,1),
#     out_temp_tf = ifelse(out_temp_scale == 99,0,1),
#     out_spatial_tf = ifelse(out_spatial_scale == "not stated",0,1),
#     out_case_df_tf = ifelse(case_def_raw == "not stated",0,1),
#     a_match_tf = ifelse(exp_out_spatial_scale_match == "unknown",0,1),
#     a_class_tf = ifelse(a_effect == "not stated",0,1),
#     a_adjusted_tf = ifelse(a_adjusted == "not stated",0,1),
#     a_uncert_def_tf = ifelse(a_uncert_def == "not stated",0,1)
#     
#   ) %>% 
#   mutate(
#     tot_score = exp_temp_tf + exp_spatial_tf + exp_source_tf + out_class_case_presentation_tf + out_spatial_tf + out_temp_tf + out_case_df_tf + a_match_tf + a_class_tf + a_adjusted_tf + a_uncert_def_tf
#   ) 

scores <- est_data %>%
  make_scores() %>%
  select(studykey, eff_country, 
         exp_class, exp_temp_scale, exp_spatial_scale, exp_source, 
         out_class_metric, out_class_case_presentation, out_medicattended, out_temp_scale, out_spatial_scale, case_def_raw,
         a_effect, a_uncert_def, 
         names(.)[str_detect(names(.), "_tf|tot_score")])

# Save
saveRDS(scores, file = "generated_data/paper_scores.rds")


# Directions of estimates -------------------------------------------------

if (opt$do_plots) {
  
  # Thresholds for ploting
  a_class_thresh <- get_estimate_thresholds()
  
  # What was sudied?
  p_exp <- est_data %>%
    # Get for each paper the unique combination of exposure and admin level
    group_by(studykey, continent, exp_class) %>%
    slice(1) %>% 
    ungroup() %>% 
    ggplot(aes(x = exp_class)) +
    geom_bar(aes(fill = continent)) +
    scale_fill_manual(values = colors_continent()) +
    guides(fill = "none") +
    theme_bw()
  
  p_admin <- est_data %>%
    # Get for each paper the unique combination of exposure and admin level
    group_by(studykey, admin_level) %>%
    slice(1) %>% 
    ungroup() %>% 
    ggplot(aes(x = admin_level)) +
    geom_bar(aes(fill = continent)) +
    scale_fill_manual(values = colors_continent()) +
    guides(fill = "none") +
    theme_bw()
  
  
  p_exp_tile <- est_data %>%
    # Get for each paper the unique combination of exposure and admin level
    group_by(studykey, exp_class, admin_level) %>%
    slice(1) %>% 
    ungroup() %>% 
    count(exp_class, admin_level) %>% 
    ggplot(aes(x = exp_class, y = admin_level)) +
    geom_tile(aes(fill = n)) +
    scale_fill_viridis_c() +
    theme_bw()
  
  
  p1 <- p_exp +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0, 0, -1, 0), units = "lines"),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(breaks = c(0, 25, 50))
  
  p2 <- p_exp_tile +
    theme(legend.position = "bottom",
          plot.margin = unit(c(0, 0, 0, 0), units = "line")) +
    scale_fill_gradient(low = "gray", high = "black") +
    labs(x = "Exposure", y = "Administrative level of association")
  
  p3 <- p_admin +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), units = "line"),
          panel.grid.minor = element_blank())  +
    scale_y_continuous(breaks = c(0, 25, 50))
  
  p_exp_comb <- p1 + plot_spacer() + p2 + p3 +
    plot_layout(ncol = 2, nrow = 2, widths = c(3, .6), heights = c(.6, 4))
  
  ggsave(p_exp_comb, filename = "figures/exp_count_by_admin_level.png", width = 8, height = 6.5)
  
  
  p_dir2 <- map(
    unique(est_data$admin_level), 
    function(x) {
      
      dat <- est_data %>%
        filter(!(a_direction %in%  c("not stated", "none")), admin_level == x) %>% 
        make_direction_counts(vars = c("admin_level", "exp_class", "continent"))
      
      if (nrow(dat) == 0) {
        return(NULL)
      }
      x_breaks_by <- case_when(
        x == "national"  ~ 20,
        x == "sub-national" ~ 25,
        T ~ 5
      )
      
      x_lim <- case_when(
        x == "national"  ~ 40,
        x == "sub-national" ~ 75,
        T ~ 5
      )
      x_breaks <- c(seq(-x_lim, x_lim, by = x_breaks_by), 0) %>% unique() %>% sort()
      
      p <- dat %>% 
        mutate(admin_level = as.character(admin_level),
               exp_class = factor(exp_class, levels = c("precipitation", "temperature", "flood", "drought"))) %>% 
        ggplot(aes(y = continent, x = n, fill = a_direction, alpha = stat_sign)) +
        geom_bar(stat = "identity", width = .5) +
        geom_vline(aes(xintercept = 0), lty = 2, lwd = .5) +
        facet_grid(admin_level ~ exp_class, switch = "y", drop = F) +
        theme_bw() +
        scale_fill_manual(values = rev(c("red", "blue"))) +
        scale_alpha_manual(values = c(.2, .4, 1), drop = F) +
        labs(x = "Number of estimates", y = "Continent", fill = "Estimate direction",
             alpha = "Statistical significance") +
        theme(panel.border = element_blank(),
              axis.line.x = element_line(),
              legend.position = "right",
              legend.direction = "vertical",
              panel.grid.minor = element_blank()) +
        scale_x_continuous(breaks = x_breaks,
                           labels = abs(x_breaks),
                           limits = c(-x_lim, x_lim))
      
      
      if (x != "regional") {
        p <- p +
          theme(strip.text.x = element_blank(),
                strip.background.x = element_blank())
      }
      
      if (x != "sub-national") {
        p <- p + 
          guides(fill = "none", alpha = "none") +
          theme(axis.title.y = element_blank())
      }
      
      if (x != "sub-national") {
        p <- p +
          theme(axis.title.x = element_blank())
      }
      
      if (x == "sub-national") {
        p <- p + theme(legend.position = c(.87, .5),
                       legend.box.background =  element_rect(fill = "white", color = "black"))
      }
      
      p
    })
  
  
  p_estimates_direction2 <- p_dir2[[3]] +
    p_dir2[[2]] +
    p_dir2[[1]] +
    plot_layout(ncol = 1, heights = c(1, 1,1.52))
  
  ggsave(p_estimates_direction2, filename = "figures/estimate_directions.png", 
         width = 11, height = 7, dpi = 300)
  
  p_timings <- est_timings  %>% 
    ggplot(aes(xmin = tl, xmax = tr, y = studykey)) +
    geom_linerange(aes(color = continent)) +
    geom_point(aes(color = continent, x = tl), size = .7) +
    geom_point(aes(color = continent, x = tr), size = .7) +
    # facet_grid(continent ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    scale_x_date(date_breaks = "5 year", 
                 date_labels = "%Y", limits = c("1970-01-01", "2022-01-01") %>% as.Date()) +
    scale_color_manual(values = colors_continent()) +
    labs(y = "study", x = "year")
  
  p_num_papers <- map_df(seq.Date(min(est_timings$tl), max(est_timings$tr), by = "1 years"), 
                         function(x) {
                           filter(est_timings, 
                                  tl <= x, tr >= x) %>% 
                             mutate(year = x)
                         }) %>% 
    mutate(continent = factor(continent, levels = c("Asia", "Africa", "Other"))) %>% 
    ggplot(aes(x = year)) +
    geom_bar(aes(fill = continent)) +
    theme_bw() +
    scale_fill_manual(values = colors_continent()) +
    guides(fill = "none") +
    scale_x_date(date_breaks = "5 year", 
                 date_labels = "%Y", limits = c("1970-01-01", "2022-01-01") %>% as.Date())
  
  p_timeline <- cowplot::plot_grid(
    p_num_papers +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = unit(c(1, 1, 0, 1), units = "lines")),
    p_timings +
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.minor.x = element_blank(),
            plot.margin = unit(c(0, 1, 0, 1), units = "lines")),
    labels = c("b", "c"),
    align = "v",
    axis = "lr",
    ncol = 1,
    rel_heights = c(.5, 3)
  ) 
  
  
  p_figure_1 <- cowplot::plot_grid(
    p_paper_admin_lev +
      theme(plot.background = element_blank(),
            legend.position = "right") +
      guides(fill = guide_colorbar("# papers")),
    p_timeline,
    rel_widths = c(1, 1),
    labels = c("a", NA),
    align = "h",
    axis = "tb",
    nrow = 1
  ) +
    theme(plot.background = element_rect(fill = "white"))
  
  ggsave(p_figure_1, filename = "figures/figure_1.png", 
         width = 14, height = 8, dpi = 400)
  
  
  # Boxplot by full score and continent
  p_corr <- est_data %>%
    # filter(a_class == "correlation") %>%
    filter(a_direction != "not stated") %>% 
    ggplot(aes(y = a_effect, x = continent)) +
    geom_boxplot(width = .4,
                 position = position_dodge(preserve = "single", width = .5), 
                 outliers = FALSE) +
    geom_jitter(aes(color = a_direction, pch = admin_level), width = .1, height = 0, alpha = .5) +
    geom_hline(data = get_estimate_thresholds(),
               aes(yintercept = y), lty = 2, lwd = .4) +
    facet_grid(a_class ~ exp_class, scales = "free_y") +
    theme_bw() +
    scale_color_manual(values = c("red", "gray", "blue")) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(y = "Correlation coefficient [-]", 
         color = "Effect direction",
         shape = "Admin level")
  
  ggsave(p_corr, filename = "figures/correlation_by_continent.png", 
         width = 14, height = 8, dpi = 300)
  
}


# List estimates on the "wrong" side of the threshold
inconsitent_direction <- est_data %>% 
  inner_join(
    get_estimate_thresholds() %>% 
      rename(threshold = y)
  ) %>% 
  mutate(a_direction_thresh = case_when(a_effect > threshold ~ "positive",
                                        a_effect < threshold ~ "negative",
                                        TRUE ~ "none")) %>% 
  filter(a_direction != "none",
         a_direction != a_direction_thresh) 

write_csv(inconsitent_direction, "generated_data/inconsistent_directions.csv")

if (opt$do_plots) {
  est_data %>% 
    ggplot(aes(y = a_effect, x = exp_class)) +
    geom_boxplot(aes(color = continent), width = .4,
                 position = position_dodge(preserve = "single", width = .5)) +
    geom_hline(data = tribble(
      ~a_class, ~y,
      "correlation", 0,
      "ratio", 0,
      "risk difference", 0
    ),
    aes(yintercept = y), lty = 2, lwd = .4) +
    facet_grid(a_class ~ exp_class, scales = "free") +
    theme_bw()
  
  
  est_data %>% 
    mutate(uncert_width = a_uncert_ub - a_uncert_lb) %>% 
    ggplot(aes(x = a_effect, y = log10(uncert_width))) +
    geom_point() +
    facet_grid(a_class ~ exp_class, scales = "free") 
}


# Scores ------------------------------------------------------------------

sum_scores_long <- scores %>% 
  group_by(studykey) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(studykey, contains("_tf")) %>% 
  pivot_longer(cols = contains("tf"),
               names_to = "what",
               values_to = "score") %>% 
  group_by(studykey) %>% 
  mutate(exp_score = sum(score[str_detect(what, "exp")], na.rm = T),
         out_score = sum(score[str_detect(what, "out")], na.rm = T),
         a_score = sum(score[str_detect(what, "a_")], na.rm = T),
         tot_score = sum(score, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(across(contains("_score"), function(x) {x == max(x)})) %>% 
  select(studykey, contains("_")) %>% 
  distinct() %>% 
  inner_join(est_data %>% distinct(studykey, continent)) %>% 
  select(studykey, continent, contains("_")) %>% 
  pivot_longer(cols = contains("_score"))  %>% 
  mutate(what = case_when(name == "tot_score" ~ "overall",
                          TRUE ~ "categories")) 


score_dict <- c("tot_score" = "total score",
                "exp_score" = "exposure",
                "out_score" = "outcome",
                "a_score" = "analysis")  


saveRDS(sum_scores_long, file = "generated_data/sum_scores_long.rds")

if (opt$do_plots) {
  
  colors_scores <- c(c("#F5C662"), "#803CE6")
  
  p_scores <- cowplot::plot_grid(
    sum_scores_long %>% 
      filter(what == "overall") %>% 
      ggplot(aes(x = continent, fill = value)) +
      geom_bar() +
      facet_grid(. ~ name, labeller = labeller(name = score_dict)) +
      theme_bw() +
      scale_fill_manual(values = colors_scores) +
      guides(fill = guide_legend("Full score")) +
      coord_flip(),
    sum_scores_long %>% 
      filter(what != "overall") %>% 
      ggplot(aes(x = continent, fill = value)) +
      geom_bar() +
      facet_grid(name ~ ., labeller = labeller(name = score_dict)) +
      theme_bw() +
      scale_fill_manual(values = colors_scores) +
      guides(fill = "none") +
      coord_flip(),
    nrow = 2,
    align = "v",
    axis = "lr",
    labels = "auto",
    rel_heights = c(.5, 1))
  
  
  ggsave(p_scores, filename = "figures/scores.png", 
         width = 7, height = 9, dpi = 300)
  
  # Full score papers
  
  full_score_papers <- sum_scores_long %>% filter(what == "overall", value)
  
  p_dir_full <- map(
    unique(est_data$admin_level), 
    function(x) {
      dat <- est_data %>%
        filter(!(a_direction %in%  c("not stated", "none")), admin_level == x,
               studykey %in% full_score_papers$studykey) %>% 
        make_direction_counts(vars = c("admin_level", "exp_class", "continent")) %>% 
        simplify_direction_levels()
      
      if (nrow(dat) == 0) {
        return(NULL)
      }
      
      x_breaks_by <- case_when(
        x == "national"  ~ 5,
        x == "sub-national" ~ 5,
        T ~ 5
      )
      
      x_lim <- case_when(
        x == "national"  ~ 15,
        x == "sub-national" ~ 25,
        T ~ 5
      )
      x_breaks <- c(seq(-x_lim, x_lim, by = x_breaks_by), 0) %>% unique() %>% sort()
      
      p <- dat %>% 
        mutate(admin_level = as.character(admin_level)) %>% 
        ggplot(aes(y = continent, x = n, fill = a_direction, alpha = stat_sign)) +
        geom_bar(stat = "identity", width = .5) +
        geom_vline(aes(xintercept = 0), lty = 2, lwd = .5) +
        facet_grid(admin_level ~ exp_class, switch = "y", drop = F) +
        theme_bw() +
        scale_fill_manual(values = rev(c("red", "blue")), drop = F) +
        scale_alpha_manual(values = c(.2, .4, 1), drop = F) +
        labs(x = "Number of estimates", y = "Continent", fill = "Estimate direction",
             alpha = "Statistical significance") +
        theme(panel.border = element_blank(),
              axis.line.x = element_line(),
              legend.position = "right",
              legend.direction = "vertical",
              panel.grid.minor = element_blank()) +
        scale_x_continuous(breaks = x_breaks,
                           labels = abs(x_breaks),
                           limits = c(-x_lim, x_lim))
      
      
      if (x != "national") {
        p <- p +
          theme(strip.text.x = element_blank(),
                strip.background.x = element_blank())
      }
      
      if (x != "sub-national") {
        p <- p + 
          guides(fill = "none", alpha = "none") +
          theme(axis.title.y = element_blank())
      }
      
      if (x != "sub-national") {
        p <- p +
          theme(axis.title.x = element_blank())
      }
      
      if (x == "sub-national") {
        p <- p + theme(legend.position = c(.87, .5),
                       legend.box.background =  element_rect(fill = "white", color = "black"))
      }
      
      p
    })
  
  
  p_estimates_direction_full_score <- p_dir_full[[2]] +
    p_dir_full[[1]] +
    plot_layout(ncol = 1, heights = c(1, 1.52))
  
  ggsave(p_estimates_direction_full_score, filename = "figures/estimate_directions_full_score.png", 
         width = 11, height = 5, dpi = 300)
  
}

# Lag -----------------------------------------------------------

dir_counts_lag <- est_data %>%
  filter(!(a_direction %in%  c("not stated", "none"))) %>% 
  make_direction_counts(vars = c("exp_class", "lag_class")) %>% 
  add_lag_class_def() %>% 
  simplify_direction_levels() 


if (opt$do_plots) {
  
  p_dir_lag <- plot_direction_by_var(df = dir_counts_lag,
                                     x_var = "lag_class",
                                     do_facet = TRUE,
                                     facet_scales = "free",
                                     facet_space = "free_x")
  
  ggsave(p_dir_lag, filename = "figures/effect_directions_lags.png", 
         width = 10, height = 6, dpi = 300)
}

# Rainy season ------------------------------------------------------------

dir_counts_rainy_season <- est_data %>%
  filter(!(a_direction %in%  c("not stated", "none"))) %>% 
  make_direction_counts(vars = c("exp_class", "rainy_season")) %>% 
  add_rainy_season_def() %>% 
  simplify_direction_levels() 


if (opt$do_plots) {
  
  p_dir_rainy <- plot_direction_by_var(df = dir_counts_rainy_season,
                                       x_var = "rainy_season",
                                       do_facet = TRUE,
                                       facet_scales = "free",
                                       facet_space = "free_x")
  
  
  
  ggsave(p_dir_rainy, filename = "figures/effect_directions_rainy_season.png", 
         width = 10, height = 6, dpi = 300)
}

# Urban/rural -------------------------------------------------------------

dir_counts_urban_rural <- est_data %>%
  filter(!(a_direction %in%  c("not stated", "none")), admin_level == "sub-national") %>% 
  make_direction_counts(vars = c("exp_class", "urban_rural_class")) %>% 
  add_urban_rural_def() %>% 
  simplify_direction_levels() 

if (opt$do_plots) {
  
  p_dir_urb <- plot_direction_by_var(df = dir_counts_urban_rural,
                                     x_var = "urban_rural_class",
                                     do_facet = TRUE,
                                     facet_scales = "free",
                                     facet_space = "free_x")
  
  
  ggsave(p_dir_urb, filename = "figures/effect_directions_urbanicity.png", 
         width = 10, height = 6, dpi = 300)
  
}


if (opt$do_plots) {
  # Combine directions
  dir_plots <- cowplot::plot_grid(
    p_dir_lag + 
      guides(fill = "none", alpha = "none"),
    p_dir_rainy + 
      guides(fill = "none", alpha = "none") +
      theme(axis.title.y = element_blank()),
    p_dir_urb  +
      theme(axis.title.y = element_blank()),
    nrow = 1,
    labels = "auto",
    rel_widths = c(1.5, 1, 1.5),
    align = "h",
    axis = "tb"
  )
  
  ggsave(dir_plots, filename = "figures/combined_directions_by_covar.png", 
         width = 14, height = 6.5, dpi = 300)
}


# Effect sizes ------------------------------------------------------------


## Correlation coefficients ----
covar_dict <- c(
  "lag_class" = "Time lag between exposure and outcome",
  "rainy_season" = "Study season",
  "urban_rural_class" = "Urbanicity"
)

corr_data <- est_data %>% 
  filter(str_detect(a_raw, "correlation")) %>% 
  correct_effect(thresh = 0) %>%
  filter(!(a_direction %in%  c("not stated", "none"))) %>% 
  define_stat_significance() %>% 
  simplify_direction_levels() %>% 
  add_lag_class_def() %>% 
  add_urban_rural_def() %>% 
  add_rainy_season_def()


if (opt$do_plots) {
  
  corr_plots <- map_plot(plot_data = corr_data,
                         y_lab = "Correlation coefficient",
                         plot_function = plot_effect_by_var,
                         plot_function_opt = list(do_facet = T,
                                                  facet_drop = T))
  
  
  p_corr <- cowplot::plot_grid(
    plotlist = corr_plots,
    nrow = 1,
    align = "h",
    axis = "tb",
    # labels = "auto",
    rel_widths = c(1.7, 1, 1.6)
  )
  
  ggsave(p_corr, filename = "figures/correlations_by_covar.png", 
         width = 14, height = 7, dpi = 300)
}


## Risk/odds ratios ----

log10_ratio_data <- est_data %>% 
  filter(str_detect(a_class, "log10_ratio")) %>% 
  correct_effect(thresh = 0) %>%
  filter(!(a_direction %in%  c("not stated", "none"))) %>% 
  define_stat_significance() %>% 
  simplify_direction_levels() %>% 
  add_lag_class_def() %>% 
  add_urban_rural_def() %>% 
  add_rainy_season_def()

if (opt$do_plots) {
  
  log10_ratio_plots <- map_plot(plot_data = log10_ratio_data,
                                y_lab = "Risk ratio (log10)",
                                plot_function = plot_effect_by_var,
                                plot_function_opt = list(do_facet = T,
                                                         facet_drop = F,
                                                         facet_scales = "free",
                                                         facet_space = "free_x",
                                                         ylim = NULL))
  
  
  p_log10_ratio <- cowplot::plot_grid(
    plotlist = log10_ratio_plots,
    nrow = 1,
    align = "h",
    axis = "tb",
    # labels = "auto",
    rel_widths = c(1.7, .7, 1.6)
  )
  
  ggsave(p_log10_ratio, filename = "figures/log10_ratios_by_covar.png", 
         width = 13, height = 7, dpi = 300)
  
}


## Risk differences ----
risk_diff_data <- est_data %>% 
  filter(str_detect(a_class, "risk difference")) %>% 
  correct_effect(thresh = 0) %>%
  filter(!(a_direction %in%  c("not stated", "none"))) %>% 
  define_stat_significance() %>% 
  simplify_direction_levels() %>% 
  add_lag_class_def() %>% 
  add_urban_rural_def() %>% 
  add_rainy_season_def()

if (opt$do_plots) {
  
  risk_diff_plots <- map_plot(vars = c("lag_class", "rainy_season"),
                              plot_data = risk_diff_data,
                              y_lab = "Risk difference",
                              plot_function = plot_effect_by_var,
                              plot_function_opt = list(do_facet = T,
                                                       facet_drop = F,
                                                       facet_scales = "free",
                                                       facet_space = "free_x",
                                                       ylim = NULL))
  
  
  p_risk_diff <- cowplot::plot_grid(
    plotlist = risk_diff_plots,
    nrow = 1,
    align = "h",
    axis = "tb",
    # labels = "auto",
    rel_widths = c(1, 1)
  )
  
  ggsave(p_risk_diff, filename = "figures/risk_diff_by_covar.png", 
         width = 9, height = 7, dpi = 300)
  
}

if (opt$do_plots) {
  # Combined plot
  p_comb <- cowplot::plot_grid(
    p_corr,
    p_log10_ratio,
    p_risk_diff,
    ncol = 1,
    labels = "auto"
  )
  
  ggsave(p_comb, filename = "figures/combined_by_covar.png", 
         width = 14, height = 20, dpi = 300)
  
}


# Save
save(corr_data, log10_ratio_data, risk_diff_data, 
     file = "generated_data/processed_effect_size_data_by_type.rdata")
