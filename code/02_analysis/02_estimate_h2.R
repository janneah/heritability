#' @title Heritability Estimation Pipeline (Simulated Data)
#' @description Runs bootstrapped heritability estimation on multiple binary phenotypes.

# ------------------------------------------------------------------
# Setup ------------------------------------------------------------
# ------------------------------------------------------------------

# Required packages (scoped usage preferred in parallel)
library(here)

# Load data --------------------------------------------------------
simdata <- readRDS(here::here("steps", "twins_dx_status_400min.rds"))

# Progress bar handler
progressr::handlers(global = TRUE)

# ------------------------------------------------------------------
# Main heritability function ---------------------------------------
# ------------------------------------------------------------------

multi_h2_se <- function(data, phenotypes, iterations, cores) {
  library(foreach)
  
  # Set up parallel backend
  doFuture::registerDoFuture()
  future::plan("multisession", workers = cores)
  
  # Set up progress bar
  progressr::handlers(global = TRUE)
  progressr::handlers(progressr::handler_progress(
    format = ":percent [:bar] :elapsed | eta: :eta",
    width = 60,
    complete = "="
  ))
  
  source(here::here("code", "00_helpers", "h2_helpers.R"))
  # Estimate p (probability of MZ given SS)
  prob <- prob_mz(data)
  
  # Loop over phenotypes
  results <- data.frame()
  
  progressr::with_progress({
    p <- progressr::progressor(along = phenotypes)
    
    results <- foreach(i = seq_along(phenotypes), .packages = c("boot", "stringr")) %dopar% {
      
      source(here::here("code", "00_helpers", "h2_helpers.R"))
      phenotype <- phenotypes[i]
      
      name <- stringr::str_remove(phenotype, "_I$")
      parameter_df <- parameters(data = data, phenotype_col = phenotype)
      
      corrs <- estim_h2_ros_rss(data, prob, phenotype)
      
      bootstrap <- boot::boot(
        data = data,
        statistic = estim_h2,
        R = iterations,
        parallel = "snow",
        phenotype = phenotype,
        p = prob
      )
      
      boot_ci <- boot::boot.ci(
        boot.out = bootstrap,
        conf = 0.95,
        type = "norm"
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
    }
  })
  
  df <- as.data.frame(do.call(rbind, results))
  colnames(df) <- c("phenotype", "cases", "controls", "h2_liab", "h2_se", "ci_95_lower", "ci_95_upper", "r_ss_liab", "r_os_liab")
  
  return(df)
}

# ------------------------------------------------------------------
# Run estimation ---------------------------------------------------
# ------------------------------------------------------------------

phenotype_cols <- setdiff(colnames(simdata), c("id", "age", "sex", "sex_config", "follow_up_end", "sib_id", "extra_ss_id"))
iterations <- 5
res <- multi_h2_se(simdata, phenotype_cols, iterations, 4)

saveRDS(res, here::here("results", paste0("h2_estimates_", iterations, "_iter.rds")))
