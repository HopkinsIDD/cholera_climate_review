# Compare results of two pipeline runs (snapshot vs current generated_data).
# Run from the repo root: Rscript Analysis/scratchpads/s04_compare_runs.R
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(cmdstanr)
})
source("Analysis/utils.R")
old_dir <- "results_snapshots/before_lag_fix_2026-11-18/generated_data/"
new_dir <- "generated_data/"
out <- list()

# Old direction rule (before the n_dir == 0 fix)
old_summarise_directions <- function(df, thresh = 0.7) {
  df %>% summarise(N = n(), n_pos = sum(a_direction == "positive"),
                   n_neg = sum(a_direction == "negative"), n_dir = n_pos + n_neg,
                   overall_dir = case_when(n_pos >= thresh * n_dir ~ "higher",
                                           n_neg >= thresh * n_dir ~ "lower",
                                           TRUE ~ "uncertain"), .groups = "drop")
}

old_est <- readRDS(str_c(old_dir, "processed_estimate_data.rds"))
new_est <- readRDS(str_c(new_dir, "processed_estimate_data.rds"))

# 1. Data -----------------------------------------------------------------
cat("\n== 1. Processed data ==\n")
cat("rows old/new:", nrow(old_est), nrow(new_est), "\n")
cat("a_class old:\n"); print(table(old_est$a_class)); cat("a_class new:\n"); print(table(new_est$a_class))
cat("pct_years_cholera_cat old/new:\n")
print(rbind(old = table(old_est$pct_years_cholera_cat), new = table(new_est$pct_years_cholera_cat)))
load(str_c(old_dir, "processed_effect_size_data_by_type.rdata")); o_es <- c(corr = nrow(corr_data), ratio = nrow(log10_ratio_data), diff = nrow(risk_diff_data))
load(str_c(new_dir, "processed_effect_size_data_by_type.rdata")); n_es <- c(corr = nrow(corr_data), ratio = nrow(log10_ratio_data), diff = nrow(risk_diff_data))
cat("effect-size datasets old/new:\n"); print(rbind(old = o_es, new = n_es))

# 2. Quality scores ---------------------------------------------------------
cat("\n== 2. Quality scores (study = min over estimates) ==\n")
sc_old <- get_score_categories(str_c(old_dir, "paper_scores.rds"))
sc_new <- get_score_categories(str_c(new_dir, "paper_scores.rds"))
cat("tot_score old:\n"); print(table(sc_old$tot_score)); cat("tot_score new:\n"); print(table(sc_new$tot_score))
sc <- full_join(sc_old, sc_new, by = "studykey", suffix = c("_old", "_new"))
cat("score_cat old x new:\n"); print(table(old = sc$score_cat_old, new = sc$score_cat_new, useNA = "ifany"))
out$scores <- sc

# 3. Direction of effect (Figure 3) ----------------------------------------
cat("\n== 3. Study-level direction (Figure 3 groups) ==\n")
grp <- c("studykey", "out_class_case_presentation", "out_class_metric", "exp_class", "continent", "admin_level")
d_old <- old_est %>% group_by(across(all_of(grp))) %>% old_summarise_directions()
d_new <- new_est %>% group_by(across(all_of(grp))) %>% summarise_directions() %>% ungroup() %>% mutate(overall_dir = as.character(overall_dir))
d <- full_join(d_old %>% select(all_of(grp), dir_old = overall_dir), d_new %>% select(all_of(grp), dir_new = overall_dir, n_dir), by = grp)
print(table(old = d$dir_old, new = d$dir_new, useNA = "ifany"))
cat("changed groups:\n"); print(d %>% filter(dir_old != dir_new) %>% select(studykey, exp_class, out_class_metric, continent, dir_old, dir_new) %>% as.data.frame())

# 4. Vote counting ----------------------------------------------------------
cat("\n== 4. Vote counting (study x exposure x outcome) ==\n")
vote <- function(est, fun) {
  est %>% group_by(studykey, exp_class, out_class_metric) %>% fun() %>% ungroup() %>%
    mutate(overall_dir = as.character(overall_dir)) %>%
    count(exp_class, out_class_metric, overall_dir) %>%
    pivot_wider(names_from = overall_dir, values_from = n, values_fill = 0)
}
v <- full_join(vote(old_est, old_summarise_directions), vote(new_est, summarise_directions),
               by = c("exp_class", "out_class_metric"), suffix = c("_old", "_new")) %>%
  filter(out_class_metric != "other") %>%
  rowwise() %>%
  mutate(p_old = if (higher_old + lower_old > 0) binom.test(max(higher_old, lower_old), higher_old + lower_old)$p.value else NA,
         p_new = if (higher_new + lower_new > 0) binom.test(max(higher_new, lower_new), higher_new + lower_new)$p.value else NA) %>%
  ungroup()
print(v %>% select(exp_class, out_class_metric, higher_old, higher_new, lower_old, lower_new, uncertain_old, uncertain_new, p_old, p_new) %>% as.data.frame(), digits = 3)

# 5. Meta-regression (Figure 4 + stratified supplement) --------------------
cat("\n== 5. Meta-regression odds ratios ==\n")
gq_files <- function(dir) dir(dir, pattern = "stan_genquant", full.names = TRUE) %>% str_subset("strat") %>% str_subset("direction")
es_old <- get_effect_sizes(gq_files(old_dir))
es_new <- get_effect_sizes(gq_files(new_dir))
key <- c("exp_class", "out_class", "stratification", "level", "what", "strata")
es <- full_join(es_old %>% mutate(strata = as.character(strata)) %>% select(all_of(key), mean, q2.5, q97.5),
                es_new %>% mutate(strata = as.character(strata)) %>% select(all_of(key), mean, q2.5, q97.5),
                by = key, suffix = c("_old", "_new")) %>%
  mutate(sig_old = q2.5_old > 1 | q97.5_old < 1, sig_new = q2.5_new > 1 | q97.5_new < 1,
         ratio = mean_new / mean_old)
cat("fits old/new:", length(gq_files(old_dir)), length(gq_files(new_dir)), "\n")
cat("coefficients only in old:", sum(is.na(es$mean_new)), " only in new:", sum(is.na(es$mean_old)), "\n")
cat("strata only in old:\n"); print(es %>% filter(is.na(mean_new)) %>% distinct(exp_class, out_class, stratification, level) %>% as.data.frame())
cat("strata only in new:\n"); print(es %>% filter(is.na(mean_old)) %>% distinct(exp_class, out_class, stratification, level) %>% as.data.frame())
both <- es %>% filter(!is.na(mean_old), !is.na(mean_new))
cat("\nmatched coefficients:", nrow(both), "; median |log OR change|:", round(median(abs(log(both$ratio))), 3), "\n")
cat("significance (95% CI excludes 1) old x new:\n"); print(table(old = both$sig_old, new = both$sig_new))
cat("\nFigure 4 (unstratified) coefficients whose significance changed:\n")
print(both %>% filter(stratification == "none", sig_old != sig_new) %>%
        mutate(what = str_remove(what, "\nref:.*") %>% str_replace("\n", " ")) %>%
        select(exp_class, out_class, what, strata, mean_old, q2.5_old, q97.5_old, mean_new, q2.5_new, q97.5_new) %>% as.data.frame(), digits = 3)
cat("\nLargest changes, Figure 4 coefficients (|log ratio| top 10):\n")
print(both %>% filter(stratification == "none") %>% arrange(desc(abs(log(ratio)))) %>% head(10) %>%
        mutate(what = str_remove(what, "\nref:.*") %>% str_replace("\n", " ")) %>%
        select(exp_class, out_class, what, strata, mean_old, mean_new, q2.5_new, q97.5_new) %>% as.data.frame(), digits = 3)
cat("\nStratified fits whose significance changed (count by stratification):\n")
print(both %>% filter(stratification != "none", sig_old != sig_new) %>% count(stratification, level, exp_class, out_class) %>% as.data.frame())
out$es <- es

# 6. I^2 --------------------------------------------------------------------
cat("\n== 6. I^2 (unstratified fits) ==\n")
i2 <- function(dir) get_i_squared(dir(dir, pattern = "stan_genquant", full.names = TRUE) %>% str_subset("direction") %>% str_subset("none-all"))
print(full_join(i2(old_dir) %>% select(exp_class, out_class, mean_old = mean),
                i2(new_dir) %>% select(exp_class, out_class, mean_new = mean), by = c("exp_class", "out_class")) %>% as.data.frame(), digits = 3)

saveRDS(out, "logs_rerun/comparison_lagfix.rds")
