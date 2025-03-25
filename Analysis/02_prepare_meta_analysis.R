# This script prepares the data for meta-analysis

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(here)
library(optparse)

source("Analysis/utils.R")

# User-supplied options
option_list <- list(
  make_option(c("-e", "--effect"), 
              default = "direction_positive", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-c", "--exp_class"), 
              default = "temperature", action ="store", type = "character", help = "Make intermediate plots."),
  make_option(c("-o", "--out_class"), 
              default = "occurrence", action ="store", type = "character", help = "Make intermediate plots.")
)

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Load data ---------------------------------------------------------------

# Load the effect sizes by type
load("generated_data/processed_effect_size_data_by_type.rdata")
est_data <- readRDS("generated_data/processed_estimate_data.rds")

full_scores <- readRDS("generated_data/sum_scores_long.rds") %>% 
  filter(what == "overall", value)

if (opt$effect == "corr") {
  data <- corr_data 
} else if (opt$effect == "ratio") {
  data <- log10_ratio_data
} else if (opt$effect == "diff") {
  data <- risk_diff_data
} else if (str_detect(opt$effect, "direction")) {
  data <- readRDS(file = "generated_data/processed_estimate_data.rds") %>% 
    filter(!(a_direction %in%  c("not stated", "none")))
} else {
  stop("Effect ", opt$effect, " not found.")
}

# Filter for exposure of interest
data <- data %>% 
  filter(exp_class == opt$exp_class,
         out_class_metric == opt$out_class)

# Fill in urban_rural_class
data <- data %>% 
  mutate(
    full_score = studykey %in% full_scores$studykey,
    full_score = factor(full_score),
    urban_rural_class = case_when(
      is.na(urban_rural_class) & admin_level != "sub-national" ~ "mixed",
      T ~ urban_rural_class
    ),
    urban_rural_class = factor(urban_rural_class, levels = c("rural", "mixed", "urban"))) %>% 
  filter(admin_level %in% c("national", "sub-national"))

if (nrow(data) == 0) {
  # https://stackoverflow.com/questions/17837289/break-exit-script
  exit <- function() { invokeRestart("abort") }    
  print("No data to run on.")
  exit()
}

# Transform data if necessary
if (opt$effect == "corr") {
  y <- z_transform(data$a_effect)
} else if (opt$effect == "direction_positive") {
  y <- data$a_direction == "positive"
} else if (opt$effect == "direction_negative") {
  y <- data$a_direction == "negative"
} else {
  y <- data$a_effect
}

# Covariate matrix --------------------------------------------------------

options(na.action='na.pass')

# Regression equation
eq <- ~ continent + lag_class + rainy_season + urban_rural_class + new_climatic_zone + hdi_cat + inland_coast

# Parse the list of covariates used to save
all_covars <- as.character(eq)[2] %>% 
  str_remove('~') %>% 
  str_split("\\+") %>% 
  .[[1]] %>% 
  str_trim()

# Make the model matrix
X <- model.matrix(eq, data = data)

# Prepare information for marginalizing
is_ind_missing <- apply(X, 1, function(x) any(is.na(x)))
ind_missing <- which(is_ind_missing)
ind_full <- which(!is_ind_missing)

get_missing_covars <- function(X, i) {
  cols <- colnames(X)[is.na(X[i, ])]
  
  str_extract_all(cols, get_all_covar_patterns()) %>% 
    unlist() %>% 
    unique()
}

make_covar_combs <- function(covar, data) {
  map(covar, ~ levels(data[[.]])) %>% 
    expand.grid() %>% 
    as_tibble() %>% 
    magrittr::set_colnames(covar)
}

# Attempt to do everything in one go
X_missing_covar_combs <- map(ind_missing, function(i) {
  miss <- get_missing_covars(X, i)
  X_combs <- make_covar_combs(covar = miss, data = data)
  
  # Complete with data that is not missing
  fake_data <- data %>% 
    slice(i) %>% 
    select(-one_of(miss)) %>% 
    bind_cols(X_combs)
  
  model.matrix(eq, data = fake_data)
})

# Define the starts and ends of fake matrix
X_comb_mat_sizes <- map_dbl(X_missing_covar_combs, nrow)
missing_covar_starts <- c(1, cumsum(X_comb_mat_sizes[-length(X_comb_mat_sizes)]) + 1)
missing_covar_ends <- cumsum(X_comb_mat_sizes)
X_missing <- reduce(X_missing_covar_combs, rbind)

# Prior is equal probability between classes. Can make this different later.
prob_missing_covar_combs <- map(ind_missing, function(i) {
  
  miss <- get_missing_covars(X, i)
  X_combs <- make_covar_combs(covar = miss, data = data)
  
  # Keep only observed combinations
  dat <- est_data
  for (x in 1:length(miss)) {
    dat <- filter(dat, !is.na(!!rlang::sym(miss[x])))
  }
  
  prop <- X_combs %>% 
    left_join(dat %>% 
                count_(vars = miss) %>% 
                mutate(prop = n/sum(n))) %>% 
    replace_na(list(prop = 1e-2)) %>% 
    mutate(prop = prop/sum(prop))
  
  prop
})

p_missing_covar <- prob_missing_covar_combs %>% 
  map(~ .$prop) %>% 
  unlist()

X[is.na(X)] <- 0


# Unique combinations -----------------------------------------------------
# continent + lag_class + rainy_season + urban_rural_class + new_climatic_zone + hdi_cat + inland_coast
covars <- as.character(eq) %>% str_split("\\+") %>% .[[2]] %>% str_trim()

# Define the matrix of all possible covariate combinations
all_u_covar <- map(covars,
                   function(x) {
                     unique(est_data[[x]]) %>% 
                       .[!is.na(.)]
                   }) %>% 
  set_names(covars) %>% 
  expand.grid() %>% 
  as_tibble()

# Make matrix of all covariate combinations
u_covar_mat <- model.matrix(data = all_u_covar, eq) %>% 
  unique()

poststrat_matrices <- map(c("all", covars), function(x){
  getPostStratMatrix(variable = x,
                     u_covar_mat = u_covar_mat,
                     u_covar_df = all_u_covar,
                     eq = eq,
                     est_data = est_data)
})

U_poststrats <- map_dbl(poststrat_matrices, nrow)
postrat_names <- map(1:(length(covars)+1), function(x) {str_c(c("all", covars)[x], rownames(poststrat_matrices[[x]]), sep = "-")}) %>% 
  unlist()

# Pad to maximum
postrat_end <- cumsum(U_poststrats)
postrat_start <- postrat_end - U_poststrats + 1

poststrat_mat_comb <- do.call(rbind, poststrat_matrices)
U_poststrat <- nrow(poststrat_mat_comb)
U_covar <- nrow(u_covar_mat)

# Maps --------------------------------------------------------------------

u_studies <- data$studykey %>% unique() %>% sort()
map_to_study <- map_dbl(data$studykey, ~ which(u_studies == .))

data$a_uncert_significant[is.na(data$a_uncert_significant)] <- FALSE
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
  M = nrow(X_missing),
  ind_full = ind_full,
  ind_missing = ind_missing,
  X_missing = X_missing,
  missing_covar_starts = missing_covar_starts,
  missing_covar_ends = missing_covar_ends,
  p_missing_covar = p_missing_covar,
  map_to_study = map_to_study,
  map_to_weight = map_to_weight,
  U_poststrat = U_poststrat,
  poststrat_mat_comb = poststrat_mat_comb,
  U_covar = U_covar,
  u_covar_mat = u_covar_mat
)

# Save
saveRDS(list(data = data,
             X = X, 
             X = X_missing,
             all_covars = all_covars,
             stan_data = stan_data,
             postrat_names = postrat_names),
        file = make_stan_data_filename(opt = opt))
