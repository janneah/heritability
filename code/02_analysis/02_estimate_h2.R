# Main script wherefrom the general script can be run on many phenotypes
library(foreach)
library(doFuture)
library(progressr)
library(profvis)


simdata <- readRDS("steps/simulated_twins.rds"))

handlers(global = TRUE)

multi_h2_se <- function(data, phenotypes, iterations, cores) {
  options(future.globals.maxSize = 7.6 * 1024^3)
  
  # Setting up cluster for parallel runs
  registerDoFuture()
  plan(multisession, workers = cores)
  
  # For progress bar
  handlers(global = TRUE)
  handlers(handler_progress(
    format = ":percent [:bar] :elapsed | eta: :eta",
    width = 60,
    complete = "="
  ))
  
  # Estimating p
  prob <- prob_mz(data)
  
  # Estimating h2 with standard errors and confidence intervals for each of the
  # phenotypes provided
  results <- data.frame()
  
  progressr::with_progress({
    p <- progressor(along = 1:length(phenotypes))
    
    results <- foreach(i = 1:length(phenotypes), .packages = "boot") %dopar% {
  
      source("code/helper_functions/h2_helper_functions.R")
      phenotype <- phenotypes[i]
      
      
      name <- stringr::str_remove(phenotype, "_I$")             # Removing "_I" from phenotype name
      parameter_df <- parameters(data = data,                   # Creating dataframe with parameters used in h2 estimation (K, Evary, T, z, I)        
                                 phenotype_col = phenotype)     # This function will also be called inside estim_h2, so this is just for saving
      
      # Running estim_h2 on its own here to get r_os and r_ss
      corrs <- estim_h2_ros_rss(data, prob, phenotype)
      
      bootstrap <- boot(data = data, 
                        statistic = estim_h2,
                        R = iterations,
                        parallel = "snow",
                        phenotype = phenotype,
                        p = prob)
      
      boot_ci <- boot.ci(boot.out = bootstrap,
                         conf = 0.95,
                         type = "norm")
      
      p(sprintf("Phenotype %s", name))
      row <- cbind(name, parameter_df$cases, parameter_df$controls, bootstrap$t0, sd(bootstrap$t), boot_ci$normal[2], boot_ci$normal[3], corrs$r_ss_liab, corrs$r_os_liab)
    }
  })
  
  df <- as.data.frame(do.call(rbind, results))
  colnames(df) <- c("phenotype", "cases", "controls", "h2_liab", "h2_se", "ci_95_lower", "ci_95_upper", "r_ss_liab", "r_os_liab")
  
  return(df)
}


iterations <- 500
res <- multi_h2_se(data, colnames(data)[9:length(colnames(data))], iterations, 25)
saveRDS(res, paste0("results/simulated_twin_10_phenotypes_", iterations, "_iter.rds"))




