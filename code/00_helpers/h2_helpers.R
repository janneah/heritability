library(tidyr)
library(lme4)
library(broom.mixed)
library(boot)
library(fMultivar)
library(stringr)
library(dplyr)

# Some functions for working with the data

#' Function to widen data - i.e., have each pair be one row, not two
#'
#' @name widen_df
#' @param data twin or sibling df where each row is an individual with a link to its twin or sibling
#' @returns twin or sibling df where each row is a twin or sibling pair
#' 

widen_df <- function(data, phenotype) {
  data <- data[, c(1:8, which(names(data) == phenotype))]
  
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
#' @returns a long df where each row is an individual, now with the random effect
#'

random_effects <- function(wide_df) {
  df_id <- wide_df %>%
    mutate(twin_id = 1:nrow(.),
           twin_id_2 = twin_id)

  df_ss <- df_id %>%
    filter(sex_config == "SS") %>%
    mutate(extra_ss_id = 1:nrow(.),
           extra_ss_id_2 = extra_ss_id)

  df_os <- df_id %>%
    filter(sex_config == "OS") %>%
    mutate(extra_ss_id = 1:nrow(.) + 8000000,
           extra_ss_id_2 = 1:nrow(.) + 16000000)

  new_wide_df <- bind_rows(df_ss, df_os) %>%
    arrange(twin_id)

  df_t1 <- new_wide_df %>% select(!ends_with("_2"))
  df_t2 <- new_wide_df %>% select(ends_with("_2"))

  colnames(df_t2) <- colnames(df_t1)

  longer_df <- bind_rows(df_t1, df_t2) %>%
    arrange(twin_id)

  return(longer_df)
}


#' Create parameters for every phenotype
#' 
#' @name parameters
#' @param data twin or sibling df in long format
#' @param phenotype_col the specific phenotype column in the data frame
#' @returns df p containing some necessary parameters for h2 estimation
#'

parameters <- function(data, phenotype_col) {
  p <- data.frame(
    cases = length(which(data[[phenotype_col]] == 1)),
    controls = length(which(data[[phenotype_col]] == 0))) %>% 
    mutate(
      phenotype = substr(phenotype_col, 1, nchar(phenotype_col)-2),
      total = cases + controls,
      K = cases / (cases + controls),
      Evary = K * (1 - K),
      T = qnorm(1 - K),
      z = dnorm(T),
      i = z / K)
  
  return(p)
}

#' Computing p, i.e. the probability of a pair being MZ given they are the same sex
#' 
#' @name prob_mz
#' @param data twin or sibling df (long format)
#' @returns the probability p of a same-sex pair being monozygotic
#' 

prob_mz <- function(data) {
  ff <- data %>% filter(sex_config == "SS" & sex == 1) %>% tally() %>% pull()
  mm <- data %>% filter(sex_config == "SS" & sex == 0) %>% tally() %>% pull()
  mf <- data %>% filter(sex_config == "OS") %>% tally() %>% pull()
  
  p_mz <- 1 - (2 * mf/(mm + ff + mf))
  p_ss <- (ff + mm)/(mm + ff + mf)
  p <- p_mz/p_ss
  
  return(p)
}

#' Fitting linear mixed model and returning h2 (for the bootstrapping)
#' 
#' @name estim_h2
#' @param data twin or sibling df (long format)
#' @param idx parameter for the boot function so it can resample the data for the bootstrapping
#' @param phenotype column name of the phenotype to estimate h2 for
#' @returns h2 estimate
#' 

estim_h2 <- function(data, idx, phenotype, p) {
  # Widening data for sampling of twin pairs
  wide_data <- widen_df(data, phenotype)
  
  # Sampling
  sampled_data <- wide_data[idx,]
  
  # Back to long format
  sampled_data_long <- random_effects(sampled_data)
  
  parameter_df <- parameters(data = sampled_data_long,
                             phenotype_col = phenotype)
  
  sampled_data_long$twin_id <- as.factor(sampled_data_long$twin_id)
  sampled_data_long$sex <- as.factor(sampled_data_long$sex)
  y <- sampled_data_long[[phenotype]]
  
  fit <- lmer(y ~ age + sex + (1 | twin_id) + (1 | extra_ss_id),
              data = sampled_data_long)
  
  tidy_fit <- tidy(fit)
  K <- parameter_df$K
  T <- parameter_df$T
  i <- parameter_df$i
  
  df <- tidy_fit[, c("group", "term", "estimate")] %>% 
    mutate(
      term = paste(group, term, sep="."),
      phenotype = substr(phenotype, 1, nchar(phenotype)-2)
    ) %>%
    select(phenotype, term, estimate) %>% 
    spread(term, estimate) %>% 
    mutate(
      var_pair = `twin_id.sd__(Intercept)` * `twin_id.sd__(Intercept)`,
      varextra_ss = `extra_ss_id.sd__(Intercept)` * `extra_ss_id.sd__(Intercept)`,
      var_res = `Residual.sd__Observation` * `Residual.sd__Observation`,
      var_tot = var_pair + varextra_ss + var_res,
      var_os = var_pair,
      var_ss = var_pair + varextra_ss,
      rss_01 = var_ss/var_tot,
      ros_01 = var_os/var_tot,
      eb1 = K + (var_ss/K),
      eb2 = K + (var_os/K),
      K = K,
      T_ss = qnorm(1 - eb1),
      T_os = qnorm(1 - eb2),
      r_ss_liab = ((T - T_ss) * sqrt(1 - (T^2 - T_ss^2) * (1 - T/i))) / (i + T_ss^2 * (i - T)),
      r_os_liab = ((T - T_os) * sqrt(1 - (T^2 - T_os^2) * (1 - T/i))) / (i + T_os^2 * (i - T)),
      h2_liab = 2/p * (r_ss_liab - r_os_liab)
    ) %>% 
    select(r_ss_liab, r_os_liab, h2_liab)
  
  return(df$h2_liab)
}

#' Function for getting h2 as well as ros and rss 

estim_h2_ros_rss <- function(data, p, phenotype) {
  # Widening data for sampling of pairs
  wide_data <- widen_df(data, phenotype)
  
  # Sampling
  sampled_data <- wide_data[1:nrow(data),]
  
  # Back to long format
  sampled_data_long <- random_effects(sampled_data)
  
  
  parameter_df <- parameters(data = sampled_data_long,
                             phenotype_col = phenotype)
  
  sampled_data_long$twin_id <- as.factor(sampled_data_long$twin_id)
  sampled_data_long$sex <- as.factor(sampled_data_long$sex)
  y <- sampled_data_long[[phenotype]]
  
  fit <- lmer(y ~ age + sex + (1 | twin_id) + (1 | extra_ss_id),
              data = sampled_data_long)
  
  tidy_fit <- tidy(fit)
  K <- parameter_df$K
  T <- parameter_df$T
  i <- parameter_df$i
  
  df <- tidy_fit[, c("group", "term", "estimate")] %>% 
    mutate(
      term = paste(group, term, sep="."),
      phenotype = substr(phenotype, 1, nchar(phenotype)-2)
    ) %>%
    select(phenotype, term, estimate) %>% 
    spread(term, estimate) %>% 
    mutate(
      var_pair = `twin_id.sd__(Intercept)` * `twin_id.sd__(Intercept)`,
      varextra_ss = `extra_ss_id.sd__(Intercept)` * `extra_ss_id.sd__(Intercept)`,
      var_res = `Residual.sd__Observation` * `Residual.sd__Observation`,
      var_tot = var_pair + varextra_ss + var_res,
      var_os = var_pair,
      var_ss = var_pair + varextra_ss,
      rss_01 = var_ss/var_tot,
      ros_01 = var_os/var_tot,
      eb1 = K + (var_ss/K),
      eb2 = K + (var_os/K),
      K = K,
      T_ss = qnorm(1 - eb1),
      T_os = qnorm(1 - eb2),
      r_ss_liab = ((T - T_ss) * sqrt(1 - (T^2 - T_ss^2) * (1 - T/i))) / (i + T_ss^2 * (i - T)),
      r_os_liab = ((T - T_os) * sqrt(1 - (T^2 - T_os^2) * (1 - T/i))) / (i + T_os^2 * (i - T)),
      h2_liab = 2/p * (r_ss_liab - r_os_liab)
    ) %>% 
    select(r_ss_liab, r_os_liab, h2_liab)
  
  return(df)
}
