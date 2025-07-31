#' @title Heritability Estimation Pipeline for Multiple Phenotypes
#' @description Main script to estimate heritability (h2) with standard errors and confidence intervals across multiple phenotypes using bootstrapping.
#' @details Designed for simulated twin data with parallel processing and progress tracking. Requires prior formatting of data and helper functions.

# ------------------------------------------------------------------
# Load data --------------------------------------------------------
# ------------------------------------------------------------------

simdata <- readRDS(here::here("steps", "simulated_twins.rds"))

# ------------------------------------------------------------------
# Heritability Estimation Function ---------------------------------
# ------------------------------------------------------------------

multi_h2_se <- function(data, phenotypes, iterations, cores) {
  
  # Register future plan for parallel execution
  doFuture::registerDoFuture()
  future::plan(future::multisession, workers = cores)
  
  # Enable progress bar
  progressr::handlers(global = TRUE)
  progressr::handlers(progressr::handler_progress(
    format = ":percent [:bar] :elapsed | eta: :eta",
    width = 60,
    complete = "="
  ))

  # Estimate p (probability same-sex pair is MZ)
  prob <- prob_mz(data)

  results <- data.frame()

  progressr::with_progress({
    p <- progressr::progressor(along = 1:length(phenotypes))

    results <- foreach::foreach(i = seq_along(phenotypes), .packages = c("boot", "stringr")) %dopar% {
      source(here::here("code", "helper_functions", "h2_helper_functions.R"))

      phenotype <- phenotypes[i]
      name <- stringr::str_remove(phenotype, "_I$")

      parameter_df <- parameters(data = data, phenotype_col = phenotype)
      corrs <- estim_h2_ros_rss(data, prob, phenotype)

      bootstrap <- boot::boot(
        data      = data,
        statistic = estim_h2,
        R         = iterations,
        parallel  = "snow",
        phenotype = phenotype,
        p         = prob
      )

      boot_ci <- boot::boot.ci(
        boot.out = bootstrap,
        conf     = 0.95,
        type     = "norm"
      )

      p(sprintf("Phenotype %s", name))

      row <- cbind(
        name,
        parameter_df$cases,
        parameter_df$controls,
        bootstrap$t0,
        sd(bootstrap$t),
        boot_ci$normal[2],
        boot_ci$normal[3],
        corrs$r_ss_liab,
        corrs$r_os_liab
      )
      row
    }
  })

  df <- as.data.frame(do.call(rbind, results))
  colnames(df) <- c("phenotype", "cases", "controls", "h2_liab", "h2_se", "ci_95_lower", "ci_95_upper", "r_ss_liab", "r_os_liab")

  return(df)
}

# ------------------------------------------------------------------
# Run heritability estimation --------------------------------------
# ------------------------------------------------------------------

phenotype_cols <- setdiff(colnames(simdata), c("id", "birth_date", "age", "sex", "sex_config", "follow_up_end", "sib_id", "extra_ss_id"))

iterations <- 500
result <- multi_h2_se(
  data       = simdata,
  phenotypes = phenotype_cols,
  iterations = iterations,
  cores      = 25
)

# ------------------------------------------------------------------
# Save results -----------------------------------------------------
# ------------------------------------------------------------------

saveRDS(result, here::here("results", sprintf("simulated_twin_10_phenotypes_%d_iter.rds", iterations)))
