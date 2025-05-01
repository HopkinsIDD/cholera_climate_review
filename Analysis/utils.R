# This script contains utility functions for the analysis and plots



# Data processing ---------------------------------------------------------

#' add_continent
#'
#' @param df 
#' @param country_col 
#'
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

#' count_countries_single
#'
#' @param df 
#' @param all_isocodes 
#'
count_countries_single <- function(df,
                                   all_isocodes) {
  map_df(all_isocodes, function(x) {
    count <- map_lgl(df$country_list, ~ x %in% .)
    tibble(country = x,
           count = sum(count))
  })
}


#' count_countries
#'
#' @param df 
#' @param all_isocodes 
#' @param group_col 
#'
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


#' char_to_date
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
char_to_date <- function(x) {
  
  out <- vector(length = length(x))
  for (i in 1:length(x)) {
    if (nchar(as.character(x[i])) == 5) {
      out[i] <- as.character(as.Date(as.numeric(x[i]), origin = "1899-12-30"))
    } else  if (nchar(as.character(x[i])) < 7) {
      out[i] <- as.character(as.Date(x[i]))
    } else {
      out[i] <- as.Date(x[i] %>% 
                          str_replace_all("-00-", "-01-") %>% 
                          str_replace_all("-00", "-01")) %>% 
        as.character()
    }
  }
  out
}

add_lag_class_def <- function(df) {
  df %>% 
    mutate(
      lag_class_def = ifelse(!is.na(lag_class), "lag specificed", "unspecified"),
      lag_class_def = factor(lag_class_def, levels = c("lag specificed", "unspecified"))
    )
}

add_rainy_season_def <- function(df) {
  df %>% 
    mutate(rainy_season_def = case_when(is.na(rainy_season) ~ "unspecified",
                                        rainy_season != "not stated" ~ "specificed", 
                                        T ~ "unspecified"),
           rainy_season_def = factor(rainy_season_def, levels = c("specificed", "unspecified")))
}

add_urban_rural_def <- function(df) {
  df %>% 
    mutate(urban_rural_class_def = case_when(is.na(urban_rural_class) ~ "unspecified",
                                             urban_rural_class != "not stated" ~ "specificed", 
                                             T ~ "unspecified"),
           urban_rural_class_def = factor(urban_rural_class_def, levels = c("specificed", "unspecified")))
}


simplify_direction_levels <- function(df) {
  df %>% 
    mutate(a_direction = factor(a_direction, levels = c("negative", "positive")))
}

define_stat_significance <- function(df) {
  df %>% 
    mutate(
      stat_sign = case_when(
        is.na(a_uncert_significant) | a_uncert_significant == "NA" ~ "not reported",
        a_uncert_significant == "TRUE" ~ "statistically significant",
        T ~ "not statistically significant"
      ),
      stat_sign = factor(stat_sign, 
                         rev(c("statistically significant", "not statistically significant", "not reported")))
    )
}

make_direction_counts <- function(df, 
                                  vars) {
  
  df %>% 
    count_(vars =  c("a_direction", "a_uncert_significant", vars)) %>%
    rowwise() %>% 
    mutate(n = ifelse(a_direction == "negative", -n, n)) %>% 
    ungroup() %>% 
    define_stat_significance()
}


correct_effect <- function(df, 
                           thresh = 0) {
  df %>% 
    mutate(a_effect_raw = a_effect,
           a_effect = case_when(a_direction == "negative" & a_effect > thresh ~ -a_effect,
                                a_direction == "positive" & a_effect < thresh ~ -a_effect,
                                T ~ a_effect),
           a_uncert_lb = case_when(a_direction == "negative" & a_effect_raw > thresh ~ -a_uncert_lb,
                                   a_direction == "positive" & a_effect_raw < thresh ~ -a_uncert_lb,
                                   T ~ a_uncert_lb),
           a_uncert_ub = case_when(a_direction == "negative" & a_effect_raw > thresh ~ -a_uncert_ub,
                                   a_direction == "positive" & a_effect_raw < thresh ~ -a_uncert_ub,
                                   T ~ a_uncert_ub))
}

get_estimate_thresholds <- function() {
  tribble(
    ~a_class, ~y,
    "correlation", 0,
    "log10_ratio", 0,
    "risk difference", 0
  )
}

get_all_covar_patterns <- function() {
  "continent|admin_level|full_score|lag_class|urban_rural_class|rainy_season|inland_coast|hdi_cat|new_climatic_zone"
}

get_all_effects <- function() {
  c("corr", "diff", "ratio", "direction_positive", "direction_negative")
}

get_all_exp_class <- function() {
  c("flood", "drought", "precipitation", "temperature")
}

get_all_out_class <- function() {
  c("occurrence", "incidence")
}

#' @name identify_endemicity
#' @export 
#' @param processed_data extracted and processed data
#' @param study_year_col the column name for the year of the study/estimate 
#' @param study_country_col the column name for the country of the study/estimate 
#' @param who_data WHO annual report data
#' @param who_year_col the column name for the year of the WHO annual reports 
#' @param who_country_col the column name for the country of the WHO annual reports
#' @param who_case_col the colname name for the cases of the WHO annual reports
#' @param prior_year_with_sCh the number of years with WHO reported cases before the study year used to define endemicity  
#' @returns processed_data with a new column called endemicity with three levels: Endemic, Non-Endemic, and Unknown

identify_endemicity <- function(
    processed_data = processed_extracted_data,
    study_year_col = "year_tl",
    study_year_lb = 1000,
    study_country_col = "eff_country",
    who_data = report %>% mutate(year = lubridate::year(time_left)),
    who_year_col = "year",
    who_country_col = "region",
    who_case_col = "suspected_cases",
    prior_year_with_sCh = 3
) {
  
  library(dplyr)
  library(tidyr)
  
  processed_data <- processed_data %>% 
    mutate(!!study_year_col := as.numeric(!!sym(study_year_col)))
  
  processed_data <- processed_data %>%
    rowwise() %>%
    mutate(
      endemicity_status = {
        study_year <- get(study_year_col)
        study_country <- get(study_country_col)
        
        if (study_year < study_year_lb + prior_year_with_sCh) {
          status <- "Unknown"
          
        } else {
          valid_years <- seq(max(study_year_lb, study_year - 100), study_year - 1)
          
          full_years_df <- tibble(
            !!who_year_col := valid_years,
            !!who_country_col := study_country
          )
          
          who_filtered <- full_years_df %>%
            left_join(who_data, by = c(who_year_col, who_country_col))
          
          if (all(is.na(who_filtered[[who_case_col]]))) {
            status <- "Unknown"
          } else {
            # Replace NA with 0 if there is any non-NA year
            who_filtered[[who_case_col]][is.na(who_filtered[[who_case_col]])] <- 0
            reported_years <- sum(who_filtered[[who_case_col]] > 0)
            
            status <- ifelse(reported_years >= prior_year_with_sCh, "Endemic", "Non-Endemic")
          }
        }
        status
      }
    ) %>%
    ungroup()
  
  return(processed_data)
}

# Data summary -----------------------------------------------

make_scores <- function(df) {
  df  %>% 
    mutate(
      exp_temp_tf = ifelse(exp_temp_scale == 99, 0, 1),
      exp_spatial_tf = ifelse(exp_spatial_scale == "not stated", 0, 1),
      exp_source_tf = ifelse(exp_source == "not stated", 0, 1),
      out_class_case_presentation_tf = ifelse(out_class_case_presentation == "not stated", 0, 1),
      out_temp_tf = ifelse(out_temp_scale == 99, 0, 1),
      out_spatial_tf = ifelse(out_spatial_scale == "not stated", 0, 1),
      out_case_df_tf = ifelse(case_def_raw == "not stated", 0, 1),
      a_match_tf = ifelse(exp_out_spatial_scale_match == "unknown", 0, 1),
      a_class_tf = ifelse(a_effect == "not stated", 0, 1),
      a_adjusted_tf = ifelse(a_adjusted == "not stated", 0, 1),
      a_uncert_def_tf = ifelse(a_uncert_def == "not stated", 0, 1),
      tot_score = exp_temp_tf + exp_spatial_tf + exp_source_tf + out_class_case_presentation_tf +
        out_spatial_tf + out_temp_tf + out_case_df_tf + a_match_tf +
        a_class_tf + a_adjusted_tf + a_uncert_def_tf
    )
}

#' @name generate_summary_stats
#' @description
#' A function to summarize basic statistics for a given exposure and outcome combination
#' @export
#' @param data the processed extracted data 
#' @param exposures the exposure over which the statistics are summarized (the default is NULL meaning across all exposure groups)
#' @param outcomes the exposure over which the statistics are summarized (the default is NULL meaning across all outcome groups)
#' @param direction_var the column name for the direction of the effect in the data
#' @param study_id_var the column name for the unique identifier for the study in the data
generate_summary_stats <- function(data,
                                   exposures = NULL,
                                   outcomes = NULL,
                                   direction_var = "a_direction",
                                   study_id_var = "studykey") {
  
  # 0. Calculate the quality score as tot_score
  data <- data %>%
    make_scores()
  
  # 1. Location-time product per study
  location_time <- data %>% mutate(year = substr(eff_tl, 1, 4)) %>% group_by(!!sym(study_id_var)) %>%
    summarise(unique_location_time_comb = n_distinct(eff_loc) * n_distinct(year)) %>%
    ungroup()
  
  # 2. Study counts by exp_class
  study_counts <- data %>% distinct(!!sym(study_id_var), exp_class) %>%
    count(exp_class, name = "n_studies")
  
  # 3. Quality summary by exposure class
  quality_vars <- c("exp_temp_scale", "exp_spatial_scale", "exp_source",
                    "out_class_case_presentation", "out_temp_scale", "out_spatial_scale",
                    "case_def_raw", "exp_out_spatial_scale_match", "a_effect",
                    "a_adjusted", "a_uncert_def")
  quality_summary_by_exp <- data %>% distinct(!!sym(study_id_var), exp_class, tot_score) %>%
    group_by(!!sym(study_id_var), exp_class) %>% 
    summarize(avg_tot_score_est = mean(tot_score))%>% 
    ungroup() %>% 
    count(exp_class, avg_tot_score_est, name = "n_studies")
  
  # 4. Quality summary by outcome class
  quality_vars <- c("exp_temp_scale", "exp_spatial_scale", "exp_source",
                    "out_class_case_presentation", "out_temp_scale", "out_spatial_scale",
                    "case_def_raw", "exp_out_spatial_scale_match", "a_effect",
                    "a_adjusted", "a_uncert_def")
  quality_summary_by_out <- data %>% distinct(!!sym(study_id_var), out_class_metric, tot_score) %>%
    group_by(!!sym(study_id_var), out_class_metric) %>% 
    summarize(avg_tot_score_est = mean(tot_score))%>% 
    ungroup() %>% 
    count(out_class_metric, avg_tot_score_est, name = "n_studies")
  
  # 5. Quality summary by continent class
  quality_vars <- c("exp_temp_scale", "exp_spatial_scale", "exp_source",
                    "out_class_case_presentation", "out_temp_scale", "out_spatial_scale",
                    "case_def_raw", "exp_out_spatial_scale_match", "a_effect",
                    "a_adjusted", "a_uncert_def")
  quality_summary_by_cont <- data %>% distinct(!!sym(study_id_var), continent, tot_score) %>%
    group_by(!!sym(study_id_var), continent) %>% 
    summarize(avg_tot_score_est = mean(tot_score))%>% 
    ungroup() %>% 
    count(continent, avg_tot_score_est, name = "n_studies")
  
  # 6. Direction of effect by quality score
  direction_by_quality <- data %>% count(tot_score, !!sym(direction_var)) %>%
    pivot_wider(names_from = !!sym(direction_var), values_from = n, values_fill = 0)
  
  # If exposure / outcome is specified in the parameter:
  filtered_data <- data
  if (!is.null(exposures)) {
    filtered_data <- filtered_data %>% filter(exp_class %in% exposures)
  }
  if (!is.null(outcomes)) {
    filtered_data <- filtered_data %>% filter(out_class_metric %in% outcomes)
  }
  
  # 7. Vote counting
  vote_combos <- expand_grid(
    exposure = if (is.null(exposures)) unique(data$exp_class) else exposures,
    outcome = if (is.null(outcomes)) unique(data$out_class_metric) else outcomes
  )
  
  # 8. Summarize study-level directions
  study_level_dir <- filtered_data %>%
    group_by(!!sym(study_id_var)) %>%
    summarise_directions(thresh = 0.7) %>%
    select(!!sym(study_id_var), overall_dir)
  
  # 9. Vote counting summary using study-level directions
  vote_combos <- expand_grid(
    exposure = if (is.null(exposures)) unique(data$exp_class) else exposures,
    outcome = if (is.null(outcomes)) unique(data$out_class_metric) else outcomes
  )
  
  vote_counting_summary <- vote_combos %>%
    rowwise() %>%
    mutate(result = list({
      subdata <- filtered_data %>%
        filter(exp_class == exposure, out_class_metric == outcome)
      
      if (nrow(subdata) == 0) {
        tibble(
          direction = NA_character_,
          p_value = NA_real_,
          total_higher_votes = 0,
          total_lower_votes = 0,
          total_uncertain_votes = 0
        )
      } else {
        study_dir_sub <- study_level_dir %>%
          filter(!!sym(study_id_var) %in% subdata[[study_id_var]])
        
        vote_counts <- study_dir_sub %>%
          count(overall_dir) %>%
          pivot_wider(names_from = overall_dir, values_from = n, values_fill = 0)
        
        n_pos <- vote_counts[["higher"]] %||% 0
        n_neg <- vote_counts[["lower"]] %||% 0
        n_uncertain <- vote_counts[["uncertain"]] %||% 0
        n_dir <- n_pos + n_neg
        
        ref_direction <- case_when(
          n_pos > n_neg ~ "higher",
          n_neg > n_pos ~ "lower",
          n_pos == n_neg ~ "higher"  # break ties in favor of higher
        )
        
        success_count <- if (ref_direction == "higher") n_pos else n_neg # Use the direction that's most frequent as the reference
        failure_count <- if (ref_direction == "higher") n_neg else n_pos # Use the direction that's most frequent as the reference
        
        p_val <- if (n_dir > 0) {
          binom.test(success_count, n_dir, p = 0.5, alternative = "two.sided")$p.value # Two sided
        } else {
          NA_real_
        }
        
        p_val_trinom <- if (n_dir > 0) {
          trinom_test(success_count, failure_count, n_uncertain)
        } else {
          NA_real_
        }
        
        tibble(
          direction = if (n_dir == 0) "no direction" else ref_direction, #only higher or lower can be direction
          p_value = p_val,
          p_value_trinom = p_val_trinom,
          total_higher_votes = n_pos,
          total_lower_votes = n_neg,
          total_uncertain_votes = n_uncertain
        )
      }
    })) %>%
    unnest(result)
  
  # 8. Summary of studies by exposure-outcome combinations
  exposure_outcome_summary <- filtered_data %>%
    distinct(!!sym(study_id_var), exp_class, out_class_metric) %>%
    count(exp_class, out_class_metric, name = "n_studies")
  
  # 9. Final output
  list(
    location_time = location_time,
    study_counts = study_counts,
    quality_summary_by_exp = quality_summary_by_exp,
    quality_summary_by_out = quality_summary_by_out, 
    quality_summary_by_cont = quality_summary_by_cont,
    direction_by_quality = direction_by_quality,
    vote_counting_summary = vote_counting_summary,
    exposure_outcome_summary = exposure_outcome_summary
  )
}

trinom_test <- function(pos, neg, ties) {
  # From stikpetR::ts_trinomial_os based on (Bian et al., 2011, p. 6)
  n <- pos + neg + ties
  p0 <- ties/n
  p1 <- (1 - p0)/2
  k <- abs(pos - neg)
  sig <- 0
  for (z in k:n) {
    for (i in 0:floor((n - z)/2)) {
      sig = sig + dmultinom(c(i, i + z, n - i - (i + z)), 
                            prob = c(p1, p1, p0))
    }
  }
  return(pmin(1, 2*sig))
}
# Direction of effect plot ------------------------------------------------

#' summaries_directions
#' Summarise directions for direction of effect plot following Boon and Thomson 2021 https://doi.org/10.1002/jrsm.1458
#' 
#' @param df 
#' @param thresh 
#'
summarise_directions <- function(df, 
                                 thresh = 0.7) {
  df %>% 
    # Apply the algorithm in Box 1 of Boon and Thomson
    summarise(N = n(),    # Number of estimates
              n_pos = sum(a_direction == "positive"),
              n_neg = sum(a_direction == "negative"),
              n_dir = n_pos + n_neg,
              overall_dir = case_when(n_pos >= thresh * n_dir ~ "positive",
                                      n_neg >= thresh * n_dir ~ "negative",
                                      TRUE ~ "uncertain") %>% 
                factor(levels = c("positive", "negative", "uncertain"),
                       labels = c("higher", "lower", "uncertain")))
}

vert_lab <- function(label){
  paste0("<vertical-text>", label, "</vertical-text>")
}


add_doe_styling <- function(gtx) {
  fill_colors <- c("low" = "#FFACA6", "med" = "#FFC28C", "high" = "#D6FFC2")
  
  gtx %>% 
    # Add colors for score quality
    tab_style(
      style = list(cell_fill(color = fill_colors["low"])),
      locations = cells_body(rows = score_cat == "low")
    ) %>% 
    tab_style(
      style = list(cell_fill(color = fill_colors["med"])),
      locations = cells_body(rows = score_cat == "medium")
    ) %>% 
    tab_style(
      style = list(cell_fill(color = fill_colors["high"])),
      locations = cells_body(rows = score_cat == "high")
    ) %>% 
    # Add bottom line for change in directions
    tab_style(
      style = cell_borders(sides = c("bottom"),  weight = px(0.5)),
      locations = cells_body(rows = bottom_line == TRUE)
    ) %>% 
    cols_hide(c("score_cat", "bottom_line", "overall_dir"))  %>% 
    tab_stub_indent(
      rows = everything(),
      indent = 3
    ) %>% 
    tab_options(
      table.font.size = 10,
      row_group.font.weight = "bold",
      data_row.padding = px(1)
    ) %>% 
    opt_css(
      css = "vertical-text {
        writing-mode: vertical-lr;
        transform: rotate(180deg);
      }"
    ) %>% 
    cols_label(
      # Covariates
      continent = md("Continent"),
      admin_level = md("Admin level"),
      out_class_case_presentation = md("Confirmation"),
      # Exposures
      precipitation = md(vert_lab("precipitation")),
      temperature = md(vert_lab("temperature")),
      flood = md(vert_lab("flood")),
      drought = md(vert_lab("drought"))
    ) %>%
    tab_spanner(
      label = "Exposures",
      columns = c(precipitation, temperature, flood, drought)
    ) %>%
    tab_spanner(
      label = "Study characteristics",
      columns = c(continent, admin_level, out_class_case_presentation)
    ) %>% 
    sub_missing(
      missing_text = "-"
    )
}




add_doe_styling2 <- function(gtx) {
  fill_colors <- c("low" = "#FFACA6", "med" = "#FFC28C", "high" = "#D6FFC2")
  
  gtx %>% 
    # Add colors for score quality
    tab_style(
      style = list(cell_fill(color = fill_colors["low"])),
      locations = cells_body(rows = score_cat == "low")
    ) %>% 
    tab_style(
      style = list(cell_fill(color = fill_colors["med"])),
      locations = cells_body(rows = score_cat == "medium")
    ) %>% 
    tab_style(
      style = list(cell_fill(color = fill_colors["high"])),
      locations = cells_body(rows = score_cat == "high")
    ) %>% 
    # Add bottom line for change in directions
    tab_style(
      style = cell_borders(sides = c("bottom"),  weight = px(0.5)),
      locations = cells_body(rows = bottom_line == TRUE)
    ) %>% 
    cols_hide(c("score_cat", "bottom_line", "overall_dir"))  %>% 
    tab_stub_indent(
      rows = everything(),
      indent = 3
    ) %>% 
    tab_options(
      table.font.size = 10,
      row_group.font.weight = "bold",
      data_row.padding = px(1)
    ) %>% 
    opt_css(
      css = "vertical-text {
        writing-mode: vertical-lr;
        transform: rotate(180deg);
      }"
    ) %>% 
    cols_label(
      # Covariates
      continent = md("Continent"),
      admin_level = md("Admin level"),
      out_class_case_presentation = md("Confirmation"),
      # Outcomes
      incidence = md(vert_lab("incidence")),
      occurrence = md(vert_lab("occurrence"))
    ) %>%
    tab_spanner(
      label = "Outcomes",
      columns = c(incidence, occurrence)
    ) %>%
    tab_spanner(
      label = "Study characteristics",
      columns = c(continent, admin_level, out_class_case_presentation)
    ) %>% 
    sub_missing(
      missing_text = "-"
    )
}
# Plots -------------------------------------------------------------------


#' colors_continent
#' Colors for continents for plots
#' @return
#'
colors_continent <- function() {
  # colors <- c("#FFA378", "#A8B545", "#8C796D")
  colors <- c("#EB9711", "#8C07DE", "#474357")
  
  names(colors) <- c("Africa", "Asia", "Other")
  colors
}



make_facet_formula <- function(x_var) {
  as.formula(str_glue("exp_class ~ {x_var}_def")) 
}

#' plot_direction_by_var
#'
#' @param df 
#' @param x_var 
#' @param facet_formula 
#' @param facet_space 
#' @param facet_scales 
#'
plot_direction_by_var <- function(df,
                                  x_var,
                                  do_facet = TRUE,
                                  facet_formula = NULL,
                                  facet_space = "free",
                                  facet_scales = "free_x",
                                  facet_drop = FALSE) {
  
  p <- df %>% 
    ggplot(aes(y = n, x = !!rlang::sym(x_var), fill = a_direction, alpha = stat_sign)) +
    geom_bar(stat = "identity", width = .5) +
    geom_hline(aes(yintercept = 0), lty = 2, lwd = .5)
  
  
  if (do_facet) {
    if (is.null(facet_formula)) {
      facet_formula <- make_facet_formula(x_var)
    }
    p <- p + facet_grid(facet_formula, space = facet_space, scales = facet_scales, 
                        switch = "y", drop = facet_drop)
  }
  
  p +
    theme_bw() +
    scale_fill_manual(values = rev(c("red", "blue")), drop = F) +
    scale_alpha_manual(values = c(.2, .4, 1), drop = F) +
    labs(x = "Lag", y = "Number of estimates", fill = "Estimate direction",
         alpha = "Statistical significance") +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(),
          legend.position = "right",
          legend.direction = "vertical",
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}  


#' plot_direction_by_var
#'
#' @param df 
#' @param x_var 
#' @param facet_formula 
#' @param facet_space 
#' @param facet_scales 
#'
plot_effect_by_var <- function(df,
                               x_var,
                               do_facet = TRUE,
                               facet_formula = NULL,
                               facet_space = "free",
                               facet_scales = "free_x",
                               facet_drop = F,
                               yintercept_dat = tibble(yintercept = 0),
                               ylim = c(-1, 1)) {
  
  p <- df %>% 
    ggplot(aes(x = !!rlang::sym(x_var), y = a_effect)) +
    geom_boxplot(aes(group = !!rlang::sym(x_var)), outliers = F,
                 width = .3) +
    geom_jitter(aes(alpha = stat_sign, color = a_direction, pch = out_class_metric), 
                width = .2, height = 0) +
    geom_hline(data = yintercept_dat,
               aes(yintercept = yintercept), lwd = .5, lty = 2) +
    labs(x = "Covariate", y = "Effect size", color = "Effect direction",
         alpha = "Statistical significance", shape = "Outcome type")
  
  if (do_facet) {
    if (is.null(facet_formula)) {
      facet_formula <- make_facet_formula(x_var)
    }
    p <- p + facet_grid(facet_formula, space = facet_space, scales = facet_scales, 
                        switch = "y", drop = facet_drop)
  }
  
  if (!is.null(ylim)) {
    p <- p +
      scale_y_continuous(limits = ylim)
  }
  
  p +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_manual(values = rev(c("red", "blue")), drop = F) +
    scale_alpha_manual(values = c(.2, .4, 1), drop = F)
}  


map_plot <- function(vars = c("lag_class", "rainy_season", "urban_rural_class"),
                     plot_data,
                     y_lab = NULL,
                     plot_function = plot_effect_by_var,
                     plot_function_opt = list(do_facet = T,
                                              facet_drop = T)) {
  
  map(
    vars, 
    function(x) {
      
      if (x == "urban_rural_class") {
        df <- filter(plot_data, admin_level == "sub-national")
      } else {
        df <- plot_data
      }
      
      p <- do.call(plot_function, 
                   c(list(df = df,
                          x_var = x), 
                     plot_function_opt)) +
        ggtitle(covar_dict[x])
      
      if (x != last(vars)) {
        p <- p + guides(color = "none", alpha = "none", shape = "none")
      }
      
      if (x != first(vars)) {
        p <- p + 
          theme(axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                strip.text.y = element_blank(),
                strip.background.y = element_blank(),
                plot.margin = unit(c(1, 0, 1, 1), units = "line"))
      } else {
        p <- p +
          theme(plot.margin = unit(c(1, 0, 1, 1), units = "line")) +
          labs(y = y_lab)
      }
      
      p
    })
}

plot_spatial_distribution_map <- function(plot_data,
                                          country_col,
                                          group_col = NULL,
                                          group_levels = NULL) {
  
  library(rgeoboundaries)
  library(sf)
  library(tidyverse)
  
  if (!is.null(group_col)) {
    country_counts <- plot_data %>% 
      distinct(studykey, !!sym(country_col), !!sym(group_col)) %>% 
      separate_rows(!!sym(country_col), sep = "\\|") %>% 
      distinct(studykey, !!sym(country_col), !!sym(group_col)) %>%
      count(!!sym(country_col), !!sym(group_col), name = "n_studies")
  } else {
    country_counts <- plot_data %>%
      distinct(studykey, !!sym(country_col)) %>%
      separate_rows(!!sym(country_col), sep = "\\|") %>%
      count(!!sym(country_col), name = "n_studies")
  }
  
  world_sf <- geoboundaries("all", adm_lvl = "adm0")
  
  # Add India's ADM0 from ADM1 union
  india_adm1 <- geoboundaries("IND", adm_lvl = "adm1")
  india_union_geom <- st_union(india_adm1)
  india_union_sf <- st_sf(
    shapeName = "India",
    shapeISO = "IND",
    geometry = india_union_geom,
    crs = st_crs(india_adm1)
  )
  world_sf_combined <- bind_rows(world_sf %>% filter(shapeISO != "IND"), india_union_sf)
  
  # Highlight affected countries
  most_affected_countries <- c("YEM", "HTI", "AFG", "COD", "NGA", "SOM", "MOZ", "ZWE", "ETH", "ZAF")
  shps_affected_countries <- rgeoboundaries::geoboundaries(
    country = most_affected_countries, adm_lvl = 'adm0'
  ) %>% 
    select(shapeISO, shapeName, geometry)
  
  if (!is.null(group_col)) {
    if (is.null(group_levels)) {
      group_levels <- sort(unique(plot_data[[group_col]]))
    }
    data_shp <- world_sf_combined %>%
      left_join(country_counts, by = c("shapeISO" = country_col)) %>%
      filter(!is.na(n_studies)) %>%
      mutate(!!sym(group_col) := factor(!!sym(group_col), levels = group_levels))
  } else {
    data_shp <- world_sf_combined %>%
      left_join(country_counts, by = c("shapeISO" = country_col)) %>%
      filter(!is.na(n_studies))
  }
  
  ggplot() +
    geom_sf(data = world_sf_combined, fill = "gray95", color = "gray95", size = 0.1) +
    geom_sf(data = data_shp, aes(fill = n_studies), color = "gray95", size = 0.1) +
    geom_sf(data = shps_affected_countries %>%
              mutate(group = "Cholera Most Affected\nCountries in 2000-2023"),
            aes(color  = group), fill=NA,linewidth = 0.1, alpha = 0.2) +
    scale_fill_viridis_c(direction = -1, name = "Number of Studies", na.value = "white") +
    scale_color_manual(name = "", values = c("Cholera Most Affected\nCountries in 2000-2023" = "red")) +
    {
      if (!is.null(group_col)) facet_grid(rows = vars(!!sym(group_col)), switch = "y") else NULL
    } +
    coord_sf(xlim = c(st_bbox(data_shp)["xmin"], st_bbox(data_shp)["xmax"]),
             ylim = c(st_bbox(data_shp)["ymin"], st_bbox(data_shp)["ymax"]),
             expand = FALSE) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text.y.left = element_text(size = 10, face = "bold"),
      strip.placement = "outside",
      plot.title = element_text(face = "bold", hjust = 0)
    )
}


# Misc --------------------------------------------------------------------

z_transform <- function(x) {
  0.5 * log((1+x)/(1-x))
}

inv_z_transform <- function(x) {
  (exp(2*x) - 1)/(exp(2*x) + 1)
}


# File names --------------------------------------------------------------
#' make_filename
#'
#' @param basename 
#' @param outdir 
#' @param opt 
#' @param parse_opt 
#' @param type 
#' @param verbose 
#'
make_filename <- function(basename, 
                          outdir = here("generated_data/"), 
                          opt = NULL, 
                          parse_opt = NULL,
                          type = "rds",
                          verbose = TRUE) {
  
  if (!is.null(parse_opt)) {
    parsed_opt <- str_c("_", parse_opt(opt))
  } else {
    parsed_opt <- ""
  }
  
  file_name <- str_glue("{outdir}{basename}{parsed_opt}.{type}")
  
  if (verbose) {
    cat("-- file name:", file_name, "\n")
  }
  
  file_name
}

parse_stan_opt <- function(opt) {
  res <- str_glue("effect-{opt$effect}_exp-{opt$exp_class}_out-{opt$out_class}")
  
  if (is.null(res)) {
    stop("Failed parsing options.")
  }
  
  res
}

#' make_inc_period_stan_data_filename
#' filename for stan data for incubation period
#' 
#' @param outdir 
#' @param opt 
#'
make_stan_data_filename <- function(outdir = here("generated_data/"),
                                    opt) {
  make_filename(basename = "stan_data",
                outdir = outdir,
                opt = opt,
                parse_opt = parse_stan_opt,
                type = "rds")
}


#' make_stan_fit_filename
#' filename for stan data for incubation period
#' 
#' @param outdir 
#' @param opt 
#'
make_stan_fit_filename <- function(outdir = here("generated_data/"),
                                   opt) {
  
  make_filename(basename = "stan_fit",
                outdir = outdir,
                opt = opt,
                parse_opt = parse_stan_opt,
                type = "rds")
}


# Post-stratification -----------------------------------------------------

#' getPostStratMatrix
#'
#' @param variable 
#' @param what 
#' @param eq 
#' @param u_covar_df 
#' @param u_covar_mat 
#' @param lab_data 
#'
#' @return
#' @export
#'
#' @examples
getPostStratMatrix <- function(variable,
                               eq,
                               u_covar_df,
                               u_covar_mat,
                               est_data) {
  
  
  # Dimension along which to post-stratify
  
  if (variable == "all") {
    u_poststrat_strata <- "all"
    U_poststrat <- 1
    poststrat_col <- NULL
  } else {
    u_poststrat_strata <- unique(u_covar_df[[variable]]) %>% sort()
    U_poststrat <- length(u_poststrat_strata)
    poststrat_col <- variable
  }
  
  if (is.na(variable)) {
    cat("-- No stratification for", variable, "in", "(not found in data)\n")
    return(matrix(1, nrow = 1, ncol = 1))
  }
  
  vars <- eq %>% 
    as.character() %>% 
    .[[2]] %>% 
    str_split("\\+") %>% .[[1]] %>% 
    str_subset(variable, negate = T) %>% 
    str_subset("1", negate = T) %>% 
    str_trim()
  
  if (length(vars) == 0 & is.null(u_poststrat_strata)) {
    cat("-- No stratification for", what, "\n")
    return(matrix(1, nrow = 1, ncol = 1))
  }
  
  # Compute post-stratified estimates for each level.
  # For categorical variables these are computed directly from the data
  # Compute frequencies
  poststrat_freq <- est_data %>% 
    group_by_at(c(poststrat_col, vars)) %>% 
    tally() %>% 
    na.omit() %>% 
    group_by_at(poststrat_col) %>% 
    mutate(frac = n/sum(n)) %>% 
    ungroup()
  
  # Make covariate matrix corresponding to frequency table
  covar_mat_poststrat_freq <- model.matrix(data = poststrat_freq, eq)
  
  poststrat_mat <- matrix(0, nrow = U_poststrat, ncol = nrow(u_covar_mat))
  counts <- rep(0, length(u_poststrat_strata))
  
  # Loop over poststrat_freq rows to match frequency tables
  for (i in 1:nrow(poststrat_freq)) {
    if (is.null(poststrat_col)) {
      ii <- 1
    } else {
      # What strata is the post-stratification for?
      ii <- which(u_poststrat_strata == poststrat_freq[[poststrat_col]][i])
    }
    
    # Get index of row in u_covar_mat that matches
    j <- map_lgl(1:nrow(u_covar_mat), function(k){
      identical(covar_mat_poststrat_freq[i,],
                u_covar_mat[k, ])
    }) %>% 
      which()
    
    if (length(j) > 0) {
      counts[ii] <- counts[ii] + 1
      # Set the frequency
      poststrat_mat[ii, j] <- poststrat_freq$frac[i]
    }
  }
  
  # Add intercept if none present
  for (i in 1:nrow(poststrat_mat)) {
    if (sum(poststrat_mat[i, ]) == 0) {
      poststrat_mat[i, 1] <- 1
    }
  }
  
  # Check that rows sum to 1
  testthat::expect_equal(rowSums(poststrat_mat), rep(1, nrow(poststrat_mat)))
  rownames(poststrat_mat) <- u_poststrat_strata
  poststrat_mat
}

# Post-processing ---------------------------------------------------------


get_effect_sizes <- function(res_files) {
  map_df(res_files, function(x) {
    stan_fit <- readRDS(x)
    effect <- str_extract(x, get_all_effects() %>% str_c(collapse = "|"))
    exp_class <- str_extract(x, get_all_exp_class() %>% str_c(collapse = "|"))
    out_class <- str_extract(x, get_all_out_class() %>% str_c(collapse = "|"))
    
    stan_data <- readRDS(str_replace(x, "stan_fit", "stan_data"))
    strata_levels <- unlist(map(stan_data$all_covars, ~ levels(stan_data$data[[.]])))
    
    stan_fit$summary("beta") %>%
      mutate(var = colnames(stan_data$X),
             what = str_extract(var, get_all_covar_patterns()),
             strata = str_remove(var, what),
             strata = factor(strata, levels = strata_levels)) %>%
      filter(!is.na(what)) %>% 
      mutate(what = map_chr(what, function(w) {
        ref <- levels(stan_data$data[[w]])[1]
        str_c(w, "\nref:", ref)
      })) %>% 
      select(what, strata, mean, q5, q95) %>% 
      mutate(effect = effect,
             exp_class = exp_class,
             out_class = out_class)
  })
}


get_poststrat <- function(res_files) {
  map_df(res_files, function(x) {
    stan_fit <- readRDS(x)
    effect <- str_extract(x, get_all_effects() %>% str_c(collapse = "|"))
    exp_class <- str_extract(x, get_all_exp_class() %>% str_c(collapse = "|"))
    out_class <- str_extract(x, get_all_out_class() %>% str_c(collapse = "|"))
    
    stan_data <- readRDS(str_replace(x, "stan_fit", "stan_data"))
    strata_levels <- unlist(map(stan_data$all_covars, ~ levels(stan_data$data[[.]])))
    
    stan_fit$summary("poststrat") %>%
      mutate(var = stan_data$postrat_names,
             what = str_split(var, "-") %>% map_chr(~.[1]),
             strata = str_split(var, "-") %>% map_chr(~.[2]),
             strata = factor(strata, levels = stan_data$postrat_names %>% str_split("-") %>% map_chr(~.[2]))) %>%
      filter(!is.na(what)) %>% 
      select(what, strata, mean, q5, q95) %>% 
      mutate(effect = effect,
             exp_class = exp_class,
             out_class = out_class)
  })
}


# Table rendering ---------------------------------------------------------

#' print_vote_table
#' Print the table of votes
#' 
#' @param stats 
#' @param this_exposure 
#'
#' @returns a gt object
#' 
print_vote_table <- function(stats, 
                             this_exposure) {
  
  vc <- stats$vote_counting_summary %>% 
    filter(exposure == this_exposure)
  
  vc %>%
    select(-exposure) %>%  # optional: drop exposure column if redundant
    rename_with(~ str_replace_all(., "_", " ") %>% str_to_title()) %>%  # pretty column names
    gt() %>%
    tab_header(title = str_glue("Vote Counting Summary for {str_to_title(this_exposure)}")) %>%
    fmt_number(
      columns = contains("P Value"),
      decimals = 4
    ) %>% 
    fmt_number(
      columns = contains("Total"),
      decimals = 0
    )
}

#' print_vote_table_continents
#' Makes a gt_group object out of a list of vote counting tables produced with print_vote_table
#' 
#' @param stats_by_continent 
#' @param this_exposure 
#'
#' @returns gt_group object
#'
print_vote_table_continents <- function(stats_by_continent, 
                                        this_exposure) {
  
  map(unique(data$continent), function(continent) {
    print_vote_table(stats = stats_by_continent[[continent]], 
                     this_exposure = this_exposure) %>%
      tab_header(title = str_glue("Vote Counting Summary for {str_to_title(this_exposure)} and continent:{continent}"))
    
  }) %>% 
    gt_group(.list = .)
  
}

#' print_vote_table_scores
#' Makes a gt_group object out of a list of vote counting tables produced with print_vote_table
#' 
#' @param stats_by_score 
#' @param this_exposure 
#'
#' @returns gt_group object
#'
print_vote_table_scores <- function(stats_by_score, 
                                    this_exposure) {
  
  map(names(stats_by_score), function(score) {
    print_vote_table(stats = stats_by_score[[score]], 
                     this_exposure = this_exposure) %>%
      tab_header(title = str_glue("Vote Counting Summary for {str_to_title(this_exposure)} and score:{score}"))
  }) %>% 
    gt_group(.list = .)
  
}
