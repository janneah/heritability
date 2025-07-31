simulate_twins_with_metadata <- function(n_pairs = 3000, prop_os = 1/3, h2_vec, K_vec, seed = 42) {
  stopifnot(length(h2_vec) == length(K_vec))
  set.seed(seed)

  library(MASS)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)

  n_os <- round(n_pairs * prop_os)
  n_ss <- n_pairs - n_os
  total_pairs <- n_os + n_ss

  # Set sex configurations (SS or OS)
  sex_config_vec <- c(rep("SS", n_ss), rep("OS", n_os))

  # Assign sex per pair
  sex_list <- map(sex_config_vec, function(cfg) {
    if (cfg == "OS") sample(c(0, 1))          # one male, one female
    else {
      sex <- sample(c(0, 1), 1, prob = c(0.51, 0.49))               # choose either male or female, 51% chance of male (as in approximately in paper)
      c(sex, sex)
    }
  })
  sex_vec <- do.call(rbind, sex_list)

  # Age range
  age_vec <- sample(18:40, total_pairs, replace = TRUE)

  # Dates for follow_up_end - arbitrarily set to 2023-12-31
  follow_up_end <- as.Date("2023-12-31")

  # Build pairwise metadata
  meta_df <- tibble(
    sib_id = rep(1:total_pairs, each = 2),
    id = 1:(2 * total_pairs),
    sex_config = rep(sex_config_vec, each = 2),
    sex = as.vector(sex_vec),
    age = rep(age_vec, each = 2),
    follow_up_end = rep(follow_up_dates, each = 2)
  )

  # Simulate each phenotype using liability model
  phenotypes <- map2_dfc(h2_vec, K_vec, function(h2, K) {
    Vg <- h2
    Ve <- 1 - h2

    sigma_ss <- matrix(c(1, Vg, Vg, 1), nrow = 2)
    sigma_os <- matrix(c(1, 0.5 * Vg, 0.5 * Vg, 1), nrow = 2)

    liab_ss <- mvrnorm(n_ss, mu = c(0, 0), Sigma = sigma_ss)
    liab_os <- mvrnorm(n_os, mu = c(0, 0), Sigma = sigma_os)
    liab <- rbind(liab_ss, liab_os)

    threshold <- qnorm(1 - K)
    y1 <- as.integer(liab[, 1] > threshold)
    y2 <- as.integer(liab[, 2] > threshold)

    tibble(pheno = c(y1, y2))
  })

  colnames(phenotypes) <- paste0("pheno_", seq_along(h2_vec), "_I")

  # Combine metadata and phenotypes
  df <- bind_cols(meta_df, phenotypes) %>%
    relocate(id, age, sex, sex_config, follow_up_end, sib_id)

  return(df)
}
