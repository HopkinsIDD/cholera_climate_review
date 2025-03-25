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
           a_effect = case_when(a_direction == "negative" & a_effect > thresh ~ - a_effect,
                                a_direction == "positive" & a_effect < thresh ~ -a_effect,
                                T ~ a_effect))
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
  "continent|admin_level|full_score|lag_class|urban_rural_class|rainy_season"
}

get_all_effects <- function() {
  c("corr", "diff", "ratio", "direction_positive", "direction_negative")
}

get_all_exp_class <- function() {
  c("flood", "drought", "precipitation", "temperature")
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


# Misc --------------------------------------------------------------------

z_transform <- function(x) {
  0.5 * log((1+x)/(1-x))
}

inv_z_transform <- function(x) {
  exp(2*x - 1)/exp(2*x + 1)
}


# File names --------------------------------------------------------------

make_output_filename <- function(effect, exposure) {
  str_glue("generated_data/stan_output_{effect}_{exposure}.rds")
}

make_stan_data_filename <- function(effect, exposure) {
  str_glue("generated_data/stan_data_{effect}_{exposure}.rds")
}



# Post-processing ---------------------------------------------------------


get_effect_sizes <- function(res_files) {
  map_df(res_files, function(x) {
    stan_fit <- readRDS(x)
    effect <- str_extract(x, get_all_effects() %>% str_c(collapse = "|"))
    exp_class <- str_extract(x, get_all_exp_class() %>% str_c(collapse = "|"))
    stan_data <- readRDS(make_stan_data_filename(effect = effect, exposure = exp_class))
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
             exp_class = exp_class)
  })
}
