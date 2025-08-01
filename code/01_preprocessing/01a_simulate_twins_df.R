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
K <- 0.02                # Population prevalence
seed <- 42
h2 = 0.8
p_mz_given_ss = 0.5

# ICD-10 and ICD-8 codes to be sampled (epilepsy in this case)
icd10 <- c("G40", "G400", "G401", "G403", "G409")
icd8  <- c("34511", "34530", "34519")

# -------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------

#' Simulate twin cohort metadata
#'
#' Generates metadata for a cohort of twin pairs, including sex configuration, age, and follow-up information. 
#' Pairs can be either same-sex (SS) or opposite-sex (OS), with configurable proportions. 
#' Each individual is assigned a unique identifier, and same-sex pairs are assigned sequential identifiers 
#' separate from those of opposite-sex pairs.
#'
#' @param n_pairs Total number of twin pairs to simulate (default: 3000)
#' @param prop_os Proportion of opposite-sex pairs among all twin pairs (default: 1/3)
#' @param seed Integer seed to ensure reproducibility (default: 42)
#'
#' @return A data frame with two rows per twin pair, including columns:
#'   \item{id}{Unique individual ID}
#'   \item{sib_id}{Unique family identifier for each twin pair}
#'   \item{sex_config}{Sex configuration of the pair ("SS" or "OS")}
#'   \item{sex}{Binary sex indicator (0 = female, 1 = male)}
#'   \item{age}{Shared age for each pair}
#'   \item{follow_up_end}{Static end of follow-up date ("2023-12-31")}
#'   \item{extra_ss_id}{Twin pair identifier, structured differently for SS and OS pairs}
#'

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
  meta_ss <- meta_data %>% dplyr::filter(sex_config == "SS")
  meta_os <- meta_data %>% dplyr::filter(sex_config == "OS")
  
  # Assign extra_ss_id to SS (sequential per pair: 1,1 → 2,2 → ...)
  meta_ss <- meta_ss %>%
    dplyr::mutate(extra_ss_id = rep(1:(n() / 2), each = 2))
  
  # Assign extra_ss_id to OS (e.g. 800001,1600001 → 800002,1600002 → ...)
  meta_os <- meta_os %>%
    dplyr::mutate(pair_index = rep(1:(n() / 2), each = 2),
           extra_ss_id = rep(c(800000, 1600000), times = n() / 2) + pair_index) %>%
    dplyr::select(-pair_index)
  
  # Combine and arrange
  meta_out <- bind_rows(meta_ss, meta_os) %>%
    dplyr::arrange(id) 
  
  return(meta_out)
}

#' Simulate diagnoses for twin data to match a target heritability (h2)
#'
#' @param meta_data A data frame with twin pairs including sex_config ("SS" or "OS") and unique id
#' @param h2_target Desired heritability on the observed scale
#' @param icd10_codes Vector of ICD-10 codes
#' @param icd8_codes Vector of ICD-8 codes
#' @param p_mz_given_ss Probability that a same-sex pair is monozygotic (default: 0.5)
#' @param prevalence Disorder prevalence (default: 0.01)
#' @param seed Random seed
#'
#' @return A data frame of diagnosis records (id, icd_source, icd, dx_date)

simulate_dx_from_h2 <- function(meta_data, h2_target, icd10_codes, icd8_codes,
                                p_mz_given_ss = 0.5, K = 0.01, seed = 42) {
  set.seed(seed)
  
  delta_r <- h2_target * p_mz_given_ss / 2
  r_os <- 0.55
  r_ss <- r_os + delta_r
  
  ss_data <- meta_data[meta_data$sex_config == "SS", ]
  os_data <- meta_data[meta_data$sex_config == "OS", ]
  
  # Helper: Assign pair types (concordant/discordant)
  assign_pair_type <- function(data, target_r) {
    n_pairs <- nrow(data) / 2
    n_conc <- round(n_pairs * target_r)
    pair_type <- sample(c(rep("concordant", n_conc), rep("discordant", n_pairs - n_conc)))
    rep(pair_type, each = 2)
  }
  
  ss_data$pair_type <- assign_pair_type(ss_data, r_ss)
  os_data$pair_type <- assign_pair_type(os_data, r_os)
  all_data <- dplyr::bind_rows(ss_data, os_data)
  
  # Target number of cases
  n_cases_target <- round(nrow(all_data) * K)
  
  # Shuffle pair order
  all_data <- all_data[order(sample(nrow(all_data))), ]
  
  # Build case list from pairs until reaching the target
  case_ids <- c()
  seen_pairs <- c()
  for (pair_id in unique(all_data$sib_id)) {
    pair <- all_data[all_data$sib_id == pair_id, ]
    type <- unique(pair$pair_type)
    
    if (type == "concordant") {
      case_val <- rbinom(1, 1, K)
      if (case_val == 1 && (length(case_ids) + 2) <= n_cases_target) {
        case_ids <- c(case_ids, pair$id)
        seen_pairs <- c(seen_pairs, pair_id)
      }
    } else {  # discordant
      sampled <- sample(pair$id, 1)
      if ((length(case_ids) + 1) <= n_cases_target) {
        case_ids <- c(case_ids, sampled)
        seen_pairs <- c(seen_pairs, pair_id)
      }
    }
    
    if (length(case_ids) >= n_cases_target) break
  }
  
  # Assign diagnoses to cases only
  n_cases <- length(case_ids)
  dx_counts <- sample(1:3, size = n_cases, replace = TRUE)
  dx_ids <- rep(case_ids, times = dx_counts)
  
  icd_sources <- sample(c("ICD10", "ICD8"), size = length(dx_ids), replace = TRUE)
  icd_codes <- purrr::map_chr(icd_sources, ~ if (.x == "ICD10") sample(icd10_codes, 1) else sample(icd8_codes, 1))
  
  dx_dates <- sample(seq(as.Date("1971-01-01"), as.Date("1995-12-31"), by = "day"),
                     size = length(dx_ids), replace = TRUE)
  
  tibble::tibble(
    id = dx_ids,
    icd_source = icd_sources,
    icd = icd_codes,
    dx_date = dx_dates
  ) |> dplyr::arrange(id, dx_date)
}



# -------------------------------------------------------------------
# Run simulation and save
# -------------------------------------------------------------------

twin_df <- simulate_twins(n_pairs = n_pairs, seed = seed)
dx_df <- simulate_dx_from_h2(
  meta_data = twin_df,
  h2_target = h2,
  icd10_codes = icd10,
  icd8_codes = icd8,
  K = K,
  seed = seed
)

# Save to disk
saveRDS(twin_df, file = "data/simulated_twin_df.rds")
saveRDS(dx_df, file = "data/simulated_dx.rds")

