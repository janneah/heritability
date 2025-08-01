#' Simulate Twin Metadata and Diagnoses
#'
#' @description
#' This script simulates metadata for n twin pairs (with realistic sex configuration),
#' and assigns diagnoses (ICD-8 or ICD-10) to a proportion of individuals
#' based on a specified population prevalence (K).
#'
#' @output
#' - `simulated_twin_df.rds`: Metadata for all simulated twins
#' - `simulated_dx.rds`: Diagnosis table for individuals with specified diagnosis
#' 
# -------------------------------------------------------------------
# Libraries
# -------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
library(tibble)

# -------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------

n_pairs <- 45000
K <- 0.01                # Population prevalence
seed <- 42

# ICD-10 and ICD-8 codes to be sampled (epilepsy in this case)
icd10 <- c("G40", "G400", "G401", "G403", "G409")
icd8  <- c("34511", "34530", "34519")

# -------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------

simulate_twins <- function(n_pairs = 3000, prop_os = 1/3, seed = 42) {
  set.seed(seed)
  
  n_os <- round(n_pairs * prop_os)
  n_ss <- n_pairs - n_os
  total_pairs <- n_os + n_ss
  
  # Sex configuration
  sex_config <- c(rep("SS", n_ss), rep("OS", n_os))
  
  # Assign sex per pair
  sex_pair <- map(sex_config, function(cfg) {
    if (cfg == "OS") sample(c(0, 1))  # one male, one female
    else {
      sex <- sample(c(0, 1), 1, prob = c(0.49, 0.51))  # ~51% male
      c(sex, sex)
    }
  })
  sex_vec <- do.call(rbind, lapply(sex_pair, function(x) matrix(x, nrow = 2)))
  
  # Assign same age per pair
  age_vec <- sample(18:40, total_pairs, replace = TRUE)
  
  # Static follow-up date
  follow_up_date <- as.Date("2023-12-31")
  
  # Metadata frame
  meta_data <- tibble(
    sib_id = rep(1:total_pairs, each = 2),
    id = 1:(2 * total_pairs),
    sex_config = rep(sex_config, each = 2),
    sex = as.vector(sex_vec),
    age = rep(age_vec, each = 2),
    follow_up_end = follow_up_date
  )
  
  # Split SS and OS pairs
  meta_ss <- meta_data %>% filter(sex_config == "SS")
  meta_os <- meta_data %>% filter(sex_config == "OS")
  
  # Assign extra_ss_id to SS (sequential per pair: 1,1 → 2,2 → ...)
  meta_ss <- meta_ss %>%
    mutate(pair_index = rep(1:(n() / 2), each = 2),
           extra_ss_id = pair_index)
  
  # Assign extra_ss_id to OS (e.g. 800001,1600001 → 800002,1600002 → ...)
  meta_os <- meta_os %>%
    mutate(pair_index = rep(1:(n() / 2), each = 2),
           extra_ss_id = rep(c(800000, 1600000), times = n() / 2) + pair_index)
  
  # Combine and arrange
  meta_out <- bind_rows(meta_ss, meta_os) %>%
    arrange(id) %>%
    select(-pair_index)
  
  return(meta_out)
}


simulate_dx <- function(meta_df, icd10_codes, icd8_codes, K = 0.01, min_dx = 1, max_dx = 3, seed = 123) {
  set.seed(seed)
  
  n_cases <- round(nrow(meta_df) * K)
  case_ids <- sample(meta_df$id, size = n_cases, replace = FALSE)
  
  # Assign 1–3 diagnosis entries per case
  dx_counts <- sample(min_dx:max_dx, size = n_cases, replace = TRUE)
  dx_ids <- rep(case_ids, times = dx_counts)
  
  # Random ICD version per entry
  icd_sources <- sample(c("ICD10", "ICD8"), size = length(dx_ids), replace = TRUE)
  icd_codes <- purrr::map_chr(icd_sources, function(source) {
    if (source == "ICD10") sample(icd10_codes, 1) else sample(icd8_codes, 1)
  })
  
  # Diagnosis dates (icd8 codes)
  followup_start <- as.Date("1971-01-01")
  followup_end <- as.Date("1995-12-31")
  dx_dates <- sample(seq(followup_start, followup_end, by = "day"), size = length(dx_ids), replace = TRUE)
  
  dx_df <- tibble(
    id = dx_ids,
    icd_source = icd_sources,
    icd = icd_codes,
    dx_date = dx_dates
  ) %>%
    arrange(id, dx_date)
  
  return(dx_df)
}

# -------------------------------------------------------------------
# Run simulation and save
# -------------------------------------------------------------------

twin_df <- simulate_twins(n_pairs = n_pairs, seed = seed)
dx_df <- simulate_dx(twin_df, icd10, icd8, K = K)

# Save to disk
saveRDS(twin_df, file = "data/simulated_twin_df.rds")
saveRDS(dx_df, file = "data/simulated_dx.rds")

