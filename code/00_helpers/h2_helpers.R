# Requires: dplyr, tidyr, lme4, broom.mixed

#' Function to widen data - i.e., have each pair be one row, not two
#'
#' @name widen_df
#' @param data twin or sibling df where each row is an individual with a link to its twin or sibling
#' @return twin or sibling df where each row is a twin or sibling pair
#' 

widen_df <- function(data, phenotype) {
  cols_to_keep <- c("id", "age", "sex", "sex_config", "follow_up_end", "sib_id", "extra_ss_id", phenotype) # Limit number of columns, when there are many phenotypes
  data <- data[, cols_to_keep]
  
  odd_idx <- seq(1, nrow(data), by = 2)
  even_idx <- seq(2, nrow(data), by = 2)
  
  twin1 <- data[odd_idx, ]
  twin2 <- data[even_idx, ]

  colnames(twin2) <- paste(colnames(twin2), "2", sep = "_")

  paired_data <- cbind(twin1, twin2)

  return(paired_data)
}


#' Function to create the random effects for the linear mixed model
#'
#' @name random_effects
#' @param wide_df takes a wide twin or sibling df and creates the random effects for each pair
#' @return a long df where each row is an individual, now with the random effect
#'

random_effects <- function(wide_df) {
  df_id <- dplyr::mutate(wide_df, sib_id = 1:nrow(wide_df), sib_id_2 = sib_id)
  
  df_ss <- df_id %>%
    dplyr::filter(sex_config == "SS") %>%

    # Random effect `extra_ss_id` is the same for SS pairs
    dplyr::mutate(extra_ss_id = 1:nrow(.), extra_ss_id_2 = extra_ss_id)
  
  df_os <- df_id %>%
    dplyr::filter(sex_config == "OS") %>%

    # Random effect `extra_ss_id` is different between OS pairs
    dplyr::mutate(extra_ss_id = 1:nrow(.) + 8000000,     
                  extra_ss_id_2 = 1:nrow(.) + 16000000)  
  
  new_wide_df <- dplyr::bind_rows(df_ss, df_os) %>%
    dplyr::arrange(sib_id)
  
  df_t1 <- dplyr::select(new_wide_df, !dplyr::ends_with("_2"))
  df_t2 <- dplyr::select(new_wide_df, dplyr::ends_with("_2"))
  colnames(df_t2) <- colnames(df_t1)
  
  longer_df <- dplyr::bind_rows(df_t1, df_t2) %>%
    dplyr::arrange(sib_id)
  
  return(longer_df)
}

#' Create parameters for every phenotype
#' 
#' @name parameters
#' @param data twin or sibling df in long format
#' @param phenotype_col the specific phenotype column in the data frame
#' @return df p containing some necessary parameters for h2 estimation
#'

parameters <- function(data, phenotype_col) {
  p <- data.frame(
    cases = sum(data[[phenotype_col]] == 1),
    controls = sum(data[[phenotype_col]] == 0)
  ) %>%
    dplyr::mutate(
      phenotype = stringr::str_remove(phenotype_col, "_I$"),
      total = cases + controls,
      K = cases / (cases + controls),
      Evary = K * (1 - K),
      T = stats::qnorm(1 - K),
      z = stats::dnorm(T),
      i = z / K
    )
  
  return(p)
}

#' Computing p, i.e. the probability of a pair being MZ given they are the same sex
#' 
#' @name prob_mz
#' @param data twin or sibling df (long format)
#' @return the probability p of a same-sex pair being monozygotic
#' 

prob_mz <- function(data) {
  ff <- data %>% dplyr::filter(sex_config == "SS", sex == 1) %>% dplyr::tally() %>% dplyr::pull()
  mm <- data %>% dplyr::filter(sex_config == "SS", sex == 0) %>% dplyr::tally() %>% dplyr::pull()
  mf <- data %>% dplyr::filter(sex_config == "OS") %>% dplyr::tally() %>% dplyr::pull()
  
  p_mz <- 1 - (2 * mf / (mm + ff + mf))
  p_ss <- (ff + mm) / (mm + ff + mf)
  p <- p_mz / p_ss
  
  return(p)
}

#' Internal function to fit LMM and compute components for h2 estimation
#'
#' @param data long-format twin/sibling data
#' @param phenotype the name of the phenotype column
#' @param p probability of a pair being MZ given same-sex
#' @return A data frame with r_ss_liab, r_os_liab, and h2_liab

compute_h2_components <- function(data, phenotype, p) {
  param_df <- parameters(data, phenotype_col = phenotype)
  
  data$sib_id <- as.factor(data$sib_id)
  data$sex <- as.factor(data$sex)
  y <- data[[phenotype]]
  
  fit <- lme4::lmer(y ~ age + sex + (1 | sib_id) + (1 | extra_ss_id), data = data)
  tidy_fit <- broom.mixed::tidy(fit)
  
  K <- param_df$K
  T <- param_df$T
  i <- param_df$i
  
  df <- tidy_fit[, c("group", "term", "estimate")] %>%
    dplyr::mutate(
      term = paste(group, term, sep = "."),
      phenotype = stringr::str_remove(phenotype, "_I$")
    ) %>%
    dplyr::select(phenotype, term, estimate) %>%
    tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
    dplyr::mutate(
      var_pair = `sib_id.sd__(Intercept)`^2,
      varextra_ss = `extra_ss_id.sd__(Intercept)`^2,
      var_res = `Residual.sd__Observation`^2,
      var_tot = var_pair + varextra_ss + var_res,
      var_os = var_pair,
      var_ss = var_pair + varextra_ss,
      rss_01 = var_ss / var_tot,
      ros_01 = var_os / var_tot,
      eb1 = K + (var_ss / K),
      eb2 = K + (var_os / K),
      T_ss = stats::qnorm(1 - eb1),
      T_os = stats::qnorm(1 - eb2),
      r_ss_liab = ((T - T_ss) * sqrt(1 - (T^2 - T_ss^2) * (1 - T/i))) / (i + T_ss^2 * (i - T)),
      r_os_liab = ((T - T_os) * sqrt(1 - (T^2 - T_os^2) * (1 - T/i))) / (i + T_os^2 * (i - T)),
      h2_liab = 2 / p * (r_ss_liab - r_os_liab)
    ) %>%
    dplyr::select(r_ss_liab, r_os_liab, h2_liab)
  
  return(df)
}

#' Returning only h2 (for the bootstrapping function)
#' 
#' @name estim_h2
#' @param data twin or sibling df (long format)
#' @param idx parameter for the boot function so it can resample the data for the bootstrapping
#' @param p prob of SS twins being MZ, output from `prob_mz`
#' @param phenotype column name of the phenotype to estimate h2 for
#' @return h2 estimate
#' 

estim_h2 <- function(data, idx, phenotype, p) {
  wide_data <- widen_df(data, phenotype)
  sampled_data <- wide_data[idx, ]
  sampled_data_long <- random_effects(sampled_data)
  
  df <- compute_h2_components(data = sampled_data_long, phenotype = phenotype, p = p)
  return(df$h2_liab)
}

#' Returning h2, ros, rss
#' 
#' @name estim_h2_ros_rss
#' @param data twin or sibling df (long format)
#' @param p prob of SS twins being MZ, output from `prob_mz`
#' @param phenotype column name of the phenotype to estimate h2 for
#' @return h2, rss, ros estimates
#' 

estim_h2_ros_rss <- function(data, p, phenotype) {
  wide_data <- widen_df(data, phenotype)
  sampled_data <- wide_data[1:nrow(data), ]  # No resampling
  sampled_data_long <- random_effects(sampled_data)
  
  df <- compute_h2_components(data = sampled_data_long, phenotype = phenotype, p = p)
  return(df)
}