# Direction of the weather-cholera association by exposure-outcome lag.
#
# Descriptive part: coverage of lag information, direction of estimates and
# study-level votes by lag class, and within-study lag profiles.
# Model part: Bernoulli-logit meta-regression of P(positive direction) on lag
# class and temporal resolution of the outcome, with a study random effect
# (same Stan model as the main direction meta-regression).
#
# Run from the repo root after 01_process_extracted_data.R:
#   Rscript Analysis/05_lag_analysis.R            # everything
#   Rscript Analysis/05_lag_analysis.R -m FALSE   # descriptive part only

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(here)
library(optparse)

source("Analysis/utils.R")

option_list <- list(
  make_option(c("-m", "--do_model"),
              default = TRUE, action = "store", type = "logical", help = "Fit the lag models"),
  make_option(c("-r", "--redo"),
              default = TRUE, action = "store", type = "logical", help = "Redo fits")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Lag classes used in the models (1 and 2 weeks are merged: too few estimates)
lag_levels_model <- c("no lag", "1-2 weeks", "1 month", "2 months", "3 months", "> 3 months")
res_levels <- c("0-14 days", "15-30 days", "> 30 days")
# Prior sd of the model coefficients: 0.5 as in the main meta-regression, and a
# wider one that lets small lag classes depart from the overall mean
prior_sds <- c(0.5, 1.5)

# Load data ---------------------------------------------------------------

est_data <- readRDS("generated_data/processed_estimate_data.rds") %>%
  filter(a_direction != "not stated") %>%
  mutate(
    a_direction = fct_drop(a_direction),
    # Is there a usable lag (exposure before or concurrent with outcome)?
    lag_status = case_when(is.na(lag_days) ~ "not stated",
                           lag_days < 0 ~ "negative lag",
                           TRUE ~ "lag stated") %>%
      factor(levels = c("lag stated", "negative lag", "not stated")),
    lag_class_model = fct_collapse(lag_class, "1-2 weeks" = c("1 week", "2 weeks")) %>%
      fct_relevel(lag_levels_model),
    # Temporal resolution of the outcome data
    out_res = case_when(out_temp_scale %in% c("31 days-1 year", "> 1 year") ~ "> 30 days",
                        out_temp_scale == "missing" ~ NA_character_,
                        TRUE ~ as.character(out_temp_scale)) %>%
      factor(levels = res_levels),
    out_res_model = if_else(out_res == "0-14 days", "<= 14 days", "> 14 days") %>%
      factor(levels = c("> 14 days", "<= 14 days"))
  )

# Estimates with a lag and a direction (positive or negative)
lag_data <- est_data %>%
  filter(lag_status == "lag stated")

dir_data <- lag_data %>%
  filter(a_direction != "none")

# 1. Coverage -------------------------------------------------------------

coverage_status <- est_data %>%
  group_by(exp_class, out_class_metric, lag_status) %>%
  summarise(n_estimates = n(),
            n_studies = n_distinct(studykey),
            .groups = "drop")

coverage_lag <- lag_data %>%
  group_by(exp_class, out_class_metric, lag_class) %>%
  summarise(n_estimates = n(),
            n_studies = n_distinct(studykey),
            .groups = "drop")

coverage_res <- lag_data %>%
  count(exp_class, out_class_metric, lag_class, out_res, name = "n_estimates")

# Number of lag classes reported by each study (within exposure and outcome)
n_lags_study <- lag_data %>%
  group_by(studykey, exp_class, out_class_metric) %>%
  summarise(n_lag_classes = n_distinct(lag_class),
            n_lag_days = n_distinct(lag_days),
            n_estimates = n(),
            .groups = "drop")

cat("\n-- Lag information by exposure and outcome\n")
print(as.data.frame(coverage_status))

# 2. Direction of estimates by lag ----------------------------------------

#' Counts and share of positive estimates by a set of variables
#' The interval is a Wilson interval on estimates, ignoring clustering by study
summarise_lag_directions <- function(df, vars) {
  df %>%
    group_by(across(all_of(vars))) %>%
    summarise(n_pos = sum(a_direction == "positive"),
              n_neg = sum(a_direction == "negative"),
              n_none = sum(a_direction == "none"),
              n_studies = n_distinct(studykey),
              .groups = "drop") %>%
    mutate(n_dir = n_pos + n_neg,
           frac_pos = n_pos/n_dir,
           ci = map2(n_pos, n_dir, function(x, n) {
             if (n == 0) return(c(NA_real_, NA_real_))
             suppressWarnings(as.numeric(prop.test(x, n, correct = FALSE)$conf.int))
           }),
           frac_pos_lo = map_dbl(ci, 1),
           frac_pos_hi = map_dbl(ci, 2)) %>%
    select(-ci)
}

dir_by_lag <- lag_data %>%
  summarise_lag_directions(c("exp_class", "out_class_metric", "lag_class"))

dir_by_lag_res <- lag_data %>%
  filter(!is.na(out_res)) %>%
  summarise_lag_directions(c("exp_class", "out_class_metric", "out_res", "lag_class"))

cat("\n-- Direction of estimates by lag (incidence)\n")
print(dir_by_lag %>% filter(out_class_metric == "incidence") %>% as.data.frame(), digits = 2)

# 3. Study-level votes by lag ---------------------------------------------

#' Vote counting by lag class: one vote per study, exposure, outcome and lag class
#' (Boon and Thomson rule in summarise_directions), with a binomial sign test
count_lag_votes <- function(df, vars) {
  df %>%
    group_by(across(all_of(c("studykey", vars)))) %>%
    summarise_directions() %>%
    ungroup() %>%
    count(across(all_of(vars)), overall_dir, .drop = FALSE) %>%
    pivot_wider(names_from = overall_dir, values_from = n, values_fill = 0) %>%
    filter(higher + lower + uncertain > 0) %>%
    rowwise() %>%
    mutate(p_binom = if (higher + lower > 0) binom.test(higher, higher + lower)$p.value else NA_real_) %>%
    ungroup()
}

votes_by_lag <- lag_data %>%
  count_lag_votes(c("exp_class", "out_class_metric", "lag_class"))

votes_by_lag_res <- lag_data %>%
  filter(!is.na(out_res)) %>%
  count_lag_votes(c("exp_class", "out_class_metric", "out_res", "lag_class"))

cat("\n-- Study-level votes by lag (incidence)\n")
print(votes_by_lag %>% filter(out_class_metric == "incidence") %>% as.data.frame(), digits = 2)

# 4. Within-study lag profiles --------------------------------------------

# Studies reporting two or more lag classes for the same exposure and outcome
profiles <- dir_data %>%
  group_by(studykey, exp_class, out_class_metric) %>%
  filter(n_distinct(lag_class) > 1) %>%
  group_by(studykey, exp_class, out_class_metric, lag_class) %>%
  summarise(n_dir = n(),
            frac_pos = mean(a_direction == "positive"),
            .groups = "drop")

# Does the direction change with lag within a study?
profile_summary <- profiles %>%
  group_by(studykey, exp_class, out_class_metric) %>%
  summarise(n_lag_classes = n(),
            profile = case_when(all(frac_pos >= .7) ~ "positive at all lags",
                                all(frac_pos <= .3) ~ "negative at all lags",
                                TRUE ~ "direction varies with lag"),
            .groups = "drop")

cat("\n-- Within-study lag profiles\n")
print(profile_summary %>% count(exp_class, out_class_metric, profile) %>% as.data.frame())

# 5. Figures --------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

# Number of estimates by direction and lag (same style as effect_directions_lags.png)
p_counts <- dir_data %>%
  filter(out_class_metric == "incidence", !is.na(out_res)) %>%
  make_direction_counts(vars = c("exp_class", "lag_class", "out_res")) %>%
  simplify_direction_levels() %>%
  plot_direction_by_var(x_var = "lag_class",
                        facet_formula = exp_class ~ out_res,
                        facet_space = "fixed",
                        facet_scales = "free_y") +
  labs(x = "Lag between exposure and outcome",
       title = "Direction of estimates by lag and temporal resolution of the outcome (incidence)")

ggsave(p_counts, filename = "figures/lag_direction_counts.png",
       width = 11, height = 8, dpi = 300)

# Share of positive estimates by lag
p_frac <- dir_by_lag %>%
  filter(out_class_metric == "incidence", n_dir > 0) %>%
  ggplot(aes(x = lag_class, y = frac_pos)) +
  geom_hline(yintercept = .5, lty = 2, lwd = .3) +
  geom_pointrange(aes(ymin = frac_pos_lo, ymax = frac_pos_hi)) +
  geom_text(aes(y = 1.08, label = str_glue("{n_dir} ({n_studies})")), size = 2.8) +
  facet_wrap(~ exp_class, nrow = 1) +
  scale_y_continuous(breaks = seq(0, 1, by = .25), limits = c(0, 1.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.minor = element_blank()) +
  labs(x = "Lag between exposure and outcome",
       y = "Share of positive estimates",
       title = "Share of positive estimates by lag (incidence)",
       caption = "Positive and negative estimates only. Bars: 95% Wilson intervals ignoring clustering by study.\nLabels: number of estimates (number of studies).")

ggsave(p_frac, filename = "figures/lag_share_positive.png",
       width = 11, height = 5, dpi = 300)

# Within-study lag profiles
p_profiles <- profiles %>%
  filter(out_class_metric == "incidence") %>%
  mutate(studykey = fct_rev(factor(studykey))) %>%
  ggplot(aes(x = lag_class, y = studykey, fill = frac_pos)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_dir), size = 2.5) +
  facet_grid(exp_class ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = .5,
                       limits = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank()) +
  labs(x = "Lag between exposure and outcome", y = NULL,
       fill = "Share of\npositive\nestimates",
       title = "Lag profiles of studies reporting two or more lag classes (incidence)",
       caption = "Numbers: positive and negative estimates in the cell.")

ggsave(p_profiles, filename = "figures/lag_within_study_profiles.png",
       width = 8, height = 10, dpi = 300)

# 6. Lag models -----------------------------------------------------------

#' Build the Stan data for the lag model of one exposure
#' Same structure as 02_prepare_meta_analysis.R, without missing covariates
make_lag_stan_data <- function(data, eq) {

  covars <- all.vars(eq)
  X <- model.matrix(eq, data = data)

  # All covariate combinations, for post-stratification
  all_u_covar <- map(covars, ~ levels(data[[.]])) %>%
    set_names(covars) %>%
    expand.grid() %>%
    as_tibble() %>%
    mutate(across(everything(), ~ factor(., levels = levels(data[[cur_column()]]))))

  u_covar_mat <- model.matrix(data = all_u_covar, eq) %>%
    unique()

  poststrat_matrices <- map(c("all", covars), function(x) {
    getPostStratMatrix(variable = x,
                       u_covar_mat = u_covar_mat,
                       u_covar_df = all_u_covar,
                       eq = eq,
                       est_data = data)
  })

  postrat_names <- map2(c("all", covars), poststrat_matrices,
                        ~ str_c(.x, rownames(.y), sep = "-")) %>%
    unlist()

  poststrat_mat_comb <- do.call(rbind, poststrat_matrices)

  u_studies <- data$studykey %>% unique() %>% sort()
  signif <- replace_na(as.character(data$a_uncert_significant) == "TRUE", FALSE)
  u_signif <- sort(unique(signif))

  stan_data <- list(
    y = as.integer(data$a_direction == "positive"),
    N = nrow(X),
    J = ncol(X),
    X = X,
    K = length(u_studies),
    L = length(u_signif),
    N_full = nrow(X),
    N_missing = 0,
    M = 0,
    ind_full = seq_len(nrow(X)),
    ind_missing = integer(0),
    X_missing = matrix(0, nrow = 0, ncol = ncol(X)),
    missing_covar_starts = integer(0),
    missing_covar_ends = integer(0),
    p_missing_covar = numeric(0),
    map_to_study = match(data$studykey, u_studies),
    map_to_weight = match(signif, u_signif),
    U_poststrat = nrow(poststrat_mat_comb),
    poststrat_mat_comb = poststrat_mat_comb,
    U_covar = nrow(u_covar_mat),
    u_covar_mat = u_covar_mat
  )

  list(data = data, X = X, stan_data = stan_data, postrat_names = postrat_names)
}

#' Fit the lag model and summarise odds ratios and post-stratified probabilities
run_lag_model <- function(data, eq, stan_model, fit_file, sd_beta = .5, redo = TRUE) {

  # Drop empty levels so that every coefficient is informed by data
  data <- data %>%
    mutate(across(all_of(all.vars(eq)), fct_drop))

  sd <- make_lag_stan_data(data = data, eq = eq)
  sd$stan_data$sd_beta <- sd_beta

  if (file.exists(fit_file) & !redo) {
    stan_fit <- readRDS(fit_file)
  } else {
    stan_fit <- stan_model$sample(
      data = sd$stan_data,
      iter_warmup = 250,
      iter_sampling = 500,
      parallel_chains = 4,
      init = .1,
      max_treedepth = 10,
      seed = 1234,
      refresh = 0
    )
    stan_fit$save_object(fit_file)
  }

  q95 <- function(x) quantile(x, c(.025, .975), names = FALSE) %>% set_names(c("q2.5", "q97.5"))

  odds_ratios <- stan_fit$summary("eff_size", mean, q95) %>%
    mutate(var = colnames(sd$X)) %>%
    filter(var != "(Intercept)") %>%
    mutate(what = str_extract(var, "lag_class_model|out_res_model"),
           strata = str_remove(var, what)) %>%
    select(what, strata, mean, q2.5, q97.5)

  poststrat <- stan_fit$summary("poststrat", mean, q95) %>%
    mutate(var = sd$postrat_names,
           what = str_extract(var, "^[^-]+"),
           strata = str_remove(var, "^[^-]+-")) %>%
    select(what, strata, mean, q2.5, q97.5)

  # Observed share of positive estimates, to compare with the model
  observed <- data %>%
    group_by(lag_class_model) %>%
    summarise(n = n(),
              n_studies = n_distinct(studykey),
              obs_frac_pos = mean(a_direction == "positive"),
              .groups = "drop")

  diagnostics <- tibble(
    n = nrow(data),
    n_studies = n_distinct(data$studykey),
    n_divergent = sum(stan_fit$diagnostic_summary(quiet = TRUE)$num_divergent),
    max_rhat = max(stan_fit$summary(c("beta", "sigma_re", "sigma_epsilon"))$rhat)
  )

  list(odds_ratios = odds_ratios,
       poststrat = poststrat,
       observed = observed,
       diagnostics = diagnostics)
}

models <- NULL

if (opt$do_model) {

  library(cmdstanr)

  # Main direction model, with the sd of the prior on coefficients as data
  stan_model <- cmdstan_model("Analysis/stan/direction_meta_analysis_lag.stan")
  eq <- ~ lag_class_model + out_res_model

  # Same inclusion criteria as the main direction meta-regression
  model_data <- dir_data %>%
    filter(out_class_metric == "incidence",
           admin_level %in% c("national", "sub-national"),
           !is.na(out_res_model))

  # Studies with two or more lag classes: within-study contrast only
  multi_lag_data <- model_data %>%
    group_by(studykey, exp_class) %>%
    filter(n_distinct(lag_class_model) > 1) %>%
    ungroup()

  runs <- expand_grid(exp_class = c("precipitation", "temperature"),
                      subset = c("all", "multilag"),
                      sd_beta = prior_sds)

  models <- pmap(runs, function(exp_class, subset, sd_beta) {
    cat("\n-- Lag model for", exp_class, "/", subset, "/ prior sd", sd_beta, "\n")

    this_data <- if (subset == "all") model_data else multi_lag_data
    this_data <- filter(this_data, exp_class == !!exp_class)

    res <- run_lag_model(
      data = this_data,
      eq = eq,
      stan_model = stan_model,
      fit_file = str_glue("generated_data/lag_stan_fit_exp-{exp_class}_subset-{subset}_sdbeta-{sd_beta}.rds"),
      sd_beta = sd_beta,
      redo = opt$redo)

    map(res, ~ mutate(., exp_class = exp_class, subset = subset, sd_beta = sd_beta, .before = 1))
  }) %>%
    purrr::transpose() %>%
    map(bind_rows)

  cat("\n-- Diagnostics\n")
  print(as.data.frame(models$diagnostics), digits = 4)
  cat("\n-- Odds ratios of a positive direction (reference: no lag, outcome resolution > 14 days)\n")
  print(as.data.frame(models$odds_ratios), digits = 3)

  # Post-stratified probability of a positive direction by lag
  p_model <- models$poststrat %>%
    filter(what == "lag_class_model") %>%
    rename(lag_class_model = strata) %>%
    left_join(models$observed, by = c("exp_class", "subset", "sd_beta", "lag_class_model")) %>%
    mutate(lag_class_model = factor(lag_class_model, levels = lag_levels_model),
           subset = factor(subset, levels = c("all", "multilag"),
                           labels = c("All studies", "Studies with 2+ lag classes")),
           prior = factor(str_glue("N(0, {sd_beta})"), levels = str_glue("N(0, {prior_sds})"))) %>%
    ggplot(aes(x = lag_class_model, y = mean)) +
    geom_hline(yintercept = .5, lty = 2, lwd = .3) +
    geom_pointrange(aes(ymin = q2.5, ymax = q97.5, color = prior), 
                    position = position_dodge(width = .5), size = .3) +
    scale_color_manual(values = c("black", "steelblue")) +
    geom_point(aes(y = obs_frac_pos), shape = 4, color = "red", size = 2) +
    geom_text(aes(y = 1.08, label = str_glue("{n} ({n_studies})")), size = 2.8) +
    facet_grid(subset ~ exp_class) +
    scale_y_continuous(breaks = seq(0, 1, by = .25), limits = c(0, 1.1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank()) +
    labs(x = "Lag between exposure and outcome",
         y = "Probability of a positive association",
         color = "Prior on\ncoefficients",
         title = "Modelled probability of a positive association by lag (incidence)",
         caption = "Points and bars: posterior mean and 95% CrI, post-stratified over the temporal resolution of the outcome.\nRed crosses: observed share of positive estimates. Labels: number of estimates (number of studies).")

  ggsave(p_model, filename = "figures/lag_model_prob_positive.png",
         width = 9, height = 7, dpi = 300)
}

# Save --------------------------------------------------------------------

saveRDS(
  list(coverage_status = coverage_status,
       coverage_lag = coverage_lag,
       coverage_res = coverage_res,
       n_lags_study = n_lags_study,
       dir_by_lag = dir_by_lag,
       dir_by_lag_res = dir_by_lag_res,
       votes_by_lag = votes_by_lag,
       votes_by_lag_res = votes_by_lag_res,
       profiles = profiles,
       profile_summary = profile_summary,
       models = models),
  file = "generated_data/lag_analysis_results.rds")
