# This script is a quick look at the data for the Jun 24 presentation


# Preamble ----------------------------------------------------------------

library(tidyverse)


# Functions ---------------------------------------------------------------

add_continent <- function(df, country_col) {
  # Define countries
  asia_countries <- c("BGD","CHN","IND","IRN","MYS","PAK","VNM")
  other_countries <- c("YEM","HTI","PER")
  
  df %>% 
    mutate(continent = case_when(admin_level != "regional" & !!rlang::sym(country_col) %in% asia_countries ~ "Asia",
                                 admin_level != "regional" & !!rlang::sym(country_col) %in% other_countries ~ "Other",
                                 admin_level != "regional" ~ "Africa",
                                 admin_level == "regional" & str_detect(!!rlang::sym(country_col), str_c(asia_countries, collapse = "|")) ~ "Asia",
                                 admin_level == "regional" & str_detect(!!rlang::sym(country_col), str_c(other_countries, collapse = "|")) ~ "Other",
                                 TRUE ~ "Africa")
    )
}

count_countries_single <- function(df,
                                   all_isocodes) {
  map_df(all_isocodes, function(x) {
    count <- map_lgl(df$country_list, ~ x %in% .)
    tibble(country = x,
           count = sum(count))
  })
}

count_countries <- function(df, 
                            all_isocodes,
                            group_col = NULL) {
  
  if (!("country_list" %in% colnames(df))) {
    stop("First extract country list")
  }
  
  if (!is.null(group_col)) {
    df %>% 
      group_by_at(group_col) %>% 
      group_modify(function(x, y) {
        count_countries_single(df = x,
                               all_isocodes = all_isocodes) 
      })
  } else {
    count_countries_single(df = df,
                           all_isocodes = all_isocodes)
  }
}

char_to_date <- function(x) {
  out <- vector(length =  length(x))
  
  for (i in 1:length(x)) {
    if (!str_detect(x[i], "[0-9]{4}-[0-9]{2}-[0-9]{2}")) {
      out[i] <- as.Date(as.numeric(as.character(x[i])), origin = "1899-12-30") %>% as.character()
    } else {
      out[i] <- x[i] %>% str_replace_all("00", "01")
    }
  }
  as.Date(out)
}

# Load data ---------------------------------------------------------------

dat_path <- here::here("Data/cholera_climate_systematic_review_extract_into_sheets.xlsx")

# Load effect extraction data
est_data <- readxl::read_excel(dat_path, sheet = "estimates_0329") %>%
  # Fix dates
  mutate(across(c("eff_tl", "eff_tr"), ~ char_to_date(.))) %>% 
  # Unpack countries for multi-country studies
  mutate(country_list = map(eff_country, ~ str_split(., pattern = "\\|")[[1]])) %>% 
  # Fix estimates
  mutate(across(c("a_effect", "a_uncert_lb", "a_uncert_ub"), ~ as.numeric(.))) %>% 
  # Fix classes of estimates
  mutate(a_class = case_when(str_detect(a_raw, "correlation") ~ "correlation",
                             a_class == "correlation" & !str_detect(a_class, "correlation") ~ "other",
                             TRUE ~ a_class))  %>%
  # Se ratio to log10 ratio for plots
  mutate(a_effect = case_when(a_class == "ratio" ~ log10(a_effect),
                              TRUE ~ a_effect),
         a_class = case_when(a_class == "ratio" ~ "log10_ratio",
                             TRUE ~ a_class)) %>% 
  # Levels
  mutate(a_direction = factor(a_direction, levels = c("negative", "none", "positive", "not stated"))) %>%
  # Spatial scales
  mutate(
    admin_level = case_when(eff_admin == "adm0" ~ "national", 
                            str_detect(eff_admin, "[1-4]") ~ "sub-national",
                            str_detect(eff_country, "\\|") | str_detect(eff_country, "multi") ~ "regional",
                            TRUE ~ "sub-national"),
    admin_level = factor(admin_level, levels = c("regional", "national", "sub-national", "other"))
  ) %>% 
  add_continent(country_col = "eff_country") 

# Metadata
meta_data <- readxl::read_excel(dat_path, sheet = "metadata_0329") %>% 
  filter(remove_keep == "keep")

# Timelines
est_timings <- est_data %>% 
  group_by(studykey, continent, eff_loc, admin_level) %>% 
  summarise(tl = min(eff_tl, na.rm = T),
            tr = max(eff_tr, na.rm = T)) %>%
  ungroup() %>% 
  mutate(dt = as.numeric(difftime(tr, tl, units = "days"))) 

# Number of papers by exposure class
est_data %>% 
  group_by(studykey, exp_class) %>% 
  slice(1) %>% 
  ungroup() %>% 
  count(exp_class)

# Paper durations
# est_timings %>%
#   ggplot(aes(x = abs(dt), fill = admin_level)) +
#   geom_histogram() +
#   facet_grid(exp_class ~ .) +
#   theme_bw() +
#   scale_x_log10(breaks = c(30, 180, 365, 365*5, 3650, 3650*5),
#                 labels = c("1 month", "6 months", "1 year", "5 years", "10 years", "50 years"),
#                 limits = c(30, max(est_timings$dt)))

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
  scale_fill_viridis_c() +
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
  scale_fill_viridis_c() +
  facet_grid(exp_class ~ admin_level) +
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



# Compute scores ----------------------------------------------------------

scores <- est_data[,colnames(est_data) %in% c("studykey","eff_country","exp_class","exp_temp_scale","exp_spatial_scale","exp_source","out_class_metric","out_class_case_presentation","out_medicattended","out_temp_scale","out_spatial_scale","case_def_raw","exp_out_spatial_scale_match","a_effect","a_adjusted","a_uncert_def")] %>% 
  mutate(
    exp_temp_tf = ifelse(exp_temp_scale == 99,0,1),
    exp_spatial_tf = ifelse(exp_spatial_scale == "not stated",0,1),
    exp_source_tf = ifelse(exp_source == "not stated",0,1),
    out_class_case_presentation_tf = ifelse(out_class_case_presentation == "not stated",0,1),
    out_temp_tf = ifelse(out_temp_scale == 99,0,1),
    out_spatial_tf = ifelse(out_spatial_scale == "not stated",0,1),
    out_case_df_tf = ifelse(case_def_raw == "not stated",0,1),
    a_match_tf = ifelse(exp_out_spatial_scale_match == "unknown",0,1),
    a_class_tf = ifelse(a_effect == "not stated",0,1),
    a_adjusted_tf = ifelse(a_adjusted == "not stated",0,1),
    a_uncert_def_tf = ifelse(a_uncert_def == "not stated",0,1)
    
  ) %>% 
  mutate(
    tot_score = exp_temp_tf + exp_spatial_tf + exp_source_tf + out_class_case_presentation_tf + out_spatial_tf + out_temp_tf + out_case_df_tf + a_match_tf + a_class_tf + a_adjusted_tf + a_uncert_def_tf
  ) 


# Add score information to estimates
est_data <- est_data %>%
  inner_join(tot_data_by_study %>% 
               mutate(full_score = avg_tot_score > 10))

# Directions of estimates -------------------------------------------------

# Thresholds for ploting
a_class_thresh <-  tribble(
  ~a_class, ~y,
  "correlation", 0,
  "log10_ratio", 0,
  "risk difference", 0
)

colors_continent <- function() {
  colors <- c("#FFA378", "#A8B545", "#8C796D")
  names(colors) <- c("Africa", "Asia", "Other")
  colors
}

# Proportion of admin levels by continent
# --> Africa has more ADM0 level effects
p_estimates_admin <- est_data  %>%
  filter(a_direction != "not stated") %>%
  count(admin_level, continent) %>% 
  group_by(admin_level) %>% 
  mutate(frac = n/sum(n)) %>% 
  ggplot(aes(x = continent, y = n, fill = continent)) +
  geom_bar(stat = "identity") +
  # geom_text(aes(label = n),
  #           position = position_fill(.5)) +
  theme_bw() +
  labs(x = "Continent", y = "Number of estimates") +
  facet_grid(admin_level ~ ., switch = "y", scales = "free") +
  scale_fill_manual(values = colors_continent()) +
  theme(
    # panel.border = element_blank(),
    axis.line.x = element_line(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.direction = "vertical")

# Bar plot of directions
p_estimates_direction <- est_data %>%
  filter(a_direction != "not stated") %>% 
  count(admin_level, exp_class, continent, a_direction) %>%
  group_by(admin_level, continent, exp_class) %>% 
  mutate(frac = n/sum(n)) %>% 
  ggplot(aes(x = continent, y = frac, fill = a_direction)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = .5), lty = 2, lwd = .5)  +
  geom_text(aes(label = n),
            position = position_fill(.5), color = "white") +
  facet_grid(admin_level ~ exp_class, switch = "y") +
  theme_bw() +
  scale_fill_manual(values = c("red", "gray", "blue")) +
  labs(x = "Continent", y = "Proportion of estimates", fill = "Estimate direction") +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(),
        legend.position = "right",
        legend.direction = "vertical")

p_timings <- est_timings %>% 
  group_by(studykey) %>% 
  mutate(min_tl = min(tl)) %>% 
  ungroup() %>% 
  mutate(studykey = factor(studykey) %>% 
           forcats::fct_reorder(min_tl)) %>% 
  filter(tl > "1900-01-01") %>% 
  ggplot(aes(xmin = tl, xmax = tr, y = studykey)) +
  geom_linerange(aes(color = continent)) +
  geom_point(aes(color = continent, x = tl), size = .7) +
  geom_point(aes(color = continent, x = tr), size = .7) +
  facet_grid(continent ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  scale_x_date(date_breaks = "5 year", 
               date_labels = "%Y", limits = c("1970-01-01", "2022-01-01") %>% as.Date()) +
  scale_color_manual(values = colors_continent()) +
  labs(y = "study", x = "year")


p_figure_1 <- cowplot::plot_grid(
  p_paper_admin_lev +
    theme(plot.background = element_blank(),
          legend.position = "right",
          plot.margin = unit(c(1, 1, 1, 1), units = "lines")) +
    guides(fill = guide_colorbar("# papers")),
  p_timings +
    theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor.x = element_blank()),
  rel_widths = c(1, .7),
  labels = "auto",
  align = "h",
  axis = "tb",
  nrow = 1
) +
  theme(plot.background = element_blank()) +
  labs(y = "study")

ggsave(p_figure_1, filename = "figures/figure_1.png", 
       width = 14, height = 8, dpi = 400)

p_figure_2 <- cowplot::plot_grid(
  p_estimates_admin +
    theme(plot.margin = unit(c(1, 1, 0, 0), units = "lines"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) ,
  p_estimates_direction +
    theme(plot.margin = unit(c(1, 1, 0, 1), units = "lines"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(fill = guide_legend("Effect direction")),
  rel_widths = c(.3, 1),
  labels = "auto",
  align = "h",
  axis = "tb",
  nrow = 1
)

ggsave(p_figure_2, filename = "figures/figure_2.png", 
       width = 14, height = 6, dpi = 300)

ggsave(p_estimates_direction +
         theme(legend.position = "bottom",
               legend.direction = "horizontal") ,
       filename = "figures/estimate_direction.png", 
       width = 8, height = 6, dpi = 300)

# Boxplot by full score and continent
p_corr <- est_data %>%
  filter(a_class == "correlation") %>%
  filter(a_direction != "not stated") %>% 
  ggplot(aes(y = a_effect, x = continent)) +
  geom_boxplot(width = .4,
               position = position_dodge(preserve = "single", width = .5), 
               outliers = FALSE) +
  geom_jitter(aes(color = a_direction, pch = admin_level), width = .1, height = 0, alpha = .5) +
  geom_hline(data = a_class_thresh,
             aes(yintercept = y), lty = 2, lwd = .4) +
  facet_grid(. ~ exp_class, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("red", "gray", "blue")) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(y = "Correlation coefficient [-]", 
       color = "Effect direction",
       shape = "Admin level")

ggsave(p_corr, filename = "figures/correlation_by_continent.png", 
       width = 12, height = 5, dpi = 300)

# Boxplot by full score and continent
est_data %>%
  filter(!is.na(full_score)) %>% 
  ggplot(aes(y = a_effect, x = continent)) +
  geom_boxplot(aes(color = full_score), width = .4,
               position = position_dodge(preserve = "single", width = .5)) +
  geom_hline(data = a_class_thresh,
             aes(yintercept = y), lty = 2, lwd = .4) +
  facet_grid(a_class ~ exp_class, scales = "free") +
  theme_bw()

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
