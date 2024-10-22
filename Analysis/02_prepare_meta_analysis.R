# This script prepares the data for meta-analysis

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(here)
library(optparse)

source("Analysis/utils.R")

# User-supplied options
option_list <- list(
  make_option(c("-e", "--effect"), 
              default = "corr", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-c", "--exp_class"), 
              default = "precipitation", action ="store", type = "character", help = "Make intermediate plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Load data ---------------------------------------------------------------

# Load the effect sizes by type
load("generated_data/processed_effect_size_data_by_type.rdata")

full_scores <- readRDS("generated_data/sum_scores_long.rds") %>% 
  filter(what == "overall", value)

if (opt$effect == "corr") {
  data <- corr_data 
} else if (opt$effect == "ratio") {
  data <- log10_ratio_data
  y <- data$a_effect
} else if (opt$effect == "diff") {
  data <- risk_diff_data
} else {
  stop("Effect ", opt$effect, " not found.")
}

# Filter for exposure of interest
data <- data %>% 
  filter(exp_class == opt$exp_class)

# Transform data if necessary
if (opt$effect == "corr") {
  y <- z_transform(data$a_effect)
} else {
  y <- data$a_effect
}

# Fill in urban_rural_class
data <- data %>% 
  mutate(
    full_score = studykey %in% full_scores$studykey,
    full_score = factor(full_score),
    urban_rural_class = case_when(
      is.na(urban_rural_class) & admin_level != "sub-national" ~ "mixed",
      T ~ urban_rural_class
    ),
    urban_rural_class = factor(urban_rural_class, levels = c("rural", "mixed", "urban")))

if (nrow(data) == 0) {
  # https://stackoverflow.com/questions/17837289/break-exit-script
  exit <- function() { invokeRestart("abort") }    
  print("No data to run on.")
  exit()
}


# Covariate matrix --------------------------------------------------------

options(na.action='na.pass')
if (length(unique(data$full_score)) > 1) {
  X <- model.matrix(~ continent + admin_level + full_score + lag_class + rainy_season + urban_rural_class, 
                    data = data)
} else {
  X <- model.matrix(~ continent + admin_level + lag_class + rainy_season + urban_rural_class, 
                    data = data)
}

# Prepare information for marginalizing
is_ind_missing <- apply(X, 1, function(x) any(is.na(x)))
ind_missing <- which(is_ind_missing)
ind_full <- which(!is_ind_missing)
missing_covar_levels <- matrix(0, nrow  = length(ind_missing), ncol = ncol(X))

counter <- 1
all_missing_covar <- c()

all_covars_pattern <- get_all_covar_patterns()
all_covars <- str_extract_all(colnames(X), all_covars_pattern) %>% unlist() %>% unique()

for (i in ind_missing) {
  # Define later
  cols <- colnames(X)[is.na(X[i, ])]
  
  missing_covar <- str_extract_all(cols, all_covars_pattern) %>% 
    unlist() %>% 
    unique()
  
  all_missing_covar <- c(all_missing_covar, missing_covar) %>% unique()
  
  for (covar in missing_covar) {
    # Get missing covariate levels
    covar_levels <- levels(data[[covar]])
    # Get columns where covariate is missing
    col_missing <- which(colnames(X) %in%  str_c(covar, covar_levels[-1]))
    X[i, col_missing] <- 0
    missing_covar_levels[counter, col_missing] <- col_missing
  }
  counter <- counter + 1
}

# Prior is equal probability between classes. Can make this different later.
p_missing_covar <- rep(0, ncol(X))

for (covar in all_missing_covar) {
  covar_levels <- levels(data[[covar]])
  col_missing <- which(colnames(X) %in%  str_c(covar, covar_levels[-1]))
  
  # Empirical counts
  prop <- data %>% 
    filter(!is.na(!!rlang::sym(covar))) %>% 
    count(!!rlang::sym(covar)) %>% 
    mutate(prop = n/sum(n)) %>% 
    pull(prop)
  
  if (length(prop) != (length(covar_levels) - 1)) {
    p_missing_covar[col_missing] <- 1/(length(covar_levels))
  } else {
    p_missing_covar[col_missing] <- prop[-1]
  }
}



# Maps --------------------------------------------------------------------

u_studies <- data$studykey %>% unique() %>% sort()
map_to_study <- map_dbl(data$studykey, ~ which(u_studies == .))

u_uncertainty <- unique(data$a_uncert_significant) %>% sort()
map_to_weight <- map_dbl(data$a_uncert_significant, ~ which(u_uncertainty == .))


# Stan data ---------------------------------------------------------------

stan_data <- list(
  y = y,
  N = length(y),
  J = ncol(X),
  X = X,
  K = length(u_studies),
  L = length(u_uncertainty),
  N_full = length(ind_full),
  N_missing = length(ind_missing),
  ind_full = ind_full,
  ind_missing = ind_missing,
  missing_covar_levels = missing_covar_levels,
  p_missing_covar = p_missing_covar,
  map_to_study = map_to_study,
  map_to_weight = map_to_weight
)

# Save
saveRDS(list(data = data,
             all_covars = all_covars,
             X = X, 
             stan_data = stan_data),
        file = make_stan_data_filename(effect = opt$effect, exposure = opt$exp_class))
