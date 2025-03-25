# This script is to validate the extracted data ensuring the accuracy of the data.

# Load packages ----
library(tidyverse)

source("Analysis/utils.R")

# Functions ----
## 1. To extract the publication year from the studykey or eff_tr column ---
extract_year <- function(time_col) {
  # Find the first numeric character position
  match_position <- regexpr("[0-9]", time_col)
  
  # Extract starting from the first numeric character
  numeric_part <- substr(time_col, match_position[1], nchar(time_col))
  
  # Extract the first four digits
  year_match <- regmatches(numeric_part, regexpr("\\d{4}", numeric_part))
  
  return(as.numeric(year_match))
}

## 2. To check all study periods are earlier than the publication date ----
check_study_period <- function(data, study_period_tr = "eff_tr", studykey = "studykey") {
  
  data[[studykey]] <- as.character(data[[studykey]])
  
  data[[study_period_tr]] <- lapply(data[[study_period_tr]], function(x) {
    if (nchar(as.character(x)) == 5) {
      return(as.Date(as.numeric(x), origin = "1899-12-30"))
    } else {
      return(x)
    }
  })
  
  # Apply extraction function to each row
  pub_years <- sapply(data[[studykey]], extract_year)
  eff_tr_years <- sapply(data[[study_period_tr]],extract_year)
  
  # Identify invalid rows where eff_tr is NOT earlier than pub_year
  invalid_rows <- which(!is.na(pub_years) & eff_tr_years > pub_years)
  
  # Print the studykey of invalid cases
  if (length(invalid_rows) > 0) {
    cat("The following rows have invalid study periods:\n -",  str_c(sort(unique(names(invalid_rows))), collapse = "\n - "))
  } else {
    cat("All study periods are valid!\n")
  }
  
}    

# 3. To check if the mean falls in the range of the 95% confidence interval ----
check_est <- function(data,
                      est_col ="a_effect", 
                      ci_lb = "a_uncert_lb", 
                      ci_ub = "a_uncert_ub"){
  
  # Extract estimates and their lower/upper bounds of the confidence interval
  ests <- data[[est_col]]
  ci_lbs <- data[[ci_lb]]
  ci_ubs <- data[[ci_ub]]
  
  # Identify invalid rows where est is not in the range of ci_lbs and ci_ubs
  invalid_rows <- which(!is.na(ests) & !is.na(ci_lbs) &!is.na(ci_ubs) & (ests < ci_lbs | ests > ci_ubs ))
  
  # Print the studykey of invalid cases
  if (length(invalid_rows) > 0) {
    cat("The following rows have invalid estimates or 95% CI:\n",  unique(names(invalid_rows)))
  } else {
    cat("All estimates are valid!\n")
  }
}

# Load extracted data ----
data <- read.csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv")
data$studykey <- iconv(data$studykey, from = "latin1", to = "UTF-8")

# Check data ----
check_study_period(data=data)
check_est(data=data%>%dplyr::filter(!is.na(a_uncert_lb) & !is.na(a_uncert_ub)))


proc_data <- readRDS("generated_data/processed_estimate_data.rds")

check_study_period(data = proc_data)

check_est(data = proc_data %>%
            dplyr::filter(!is.na(a_uncert_lb) & !is.na(a_uncert_ub)) %>%
            correct_effect())
