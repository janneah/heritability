#' Main script
#' 
#' @author JAH
#' @date 2024-02-29
#'
#' @description
#'

source("code/aux/functions.R")

multi_h2 <- function(twin_data, phenotype_list, bootstrap_iter){
    
    # Estimating p - the probability of being MZ given SS
    p <- prob_mz(twin_data)

    # Initiating output df
    res <- data.frame(
        phenotype = NULL,
        cases = NULL,
        controls = NULL,
        K = NULL,
        h2_liab = NULL,
        h2_se = NULL,
        ci_95_lower = NULL,
        ci_95_upper = NULL
    )

    # Estimating h2 and bootstrapping to get 95% CIs for each phenotype
    for (i in 1:length(phenotype_list)){

        phenotype <- phenotype_list[i]
        name <- substr(phenotype, 1, nchar(phenotype)-2)    # Removing "_I" from phenotype name

        # Computing parameters for estim_h2
        h2_parameters <- parameters(
            data = twin_data,
            phenotype_col = phenotype
        )

        # Bootstrapping h2 using estim_h2 and boot package
        bootstrap <- boot(
            data = twin_data,                               # Input data
            statistic = estim_h2,                           # Own function to estimate h2
            R = bootstrap_iter,                             # Number of iterations
            parallel = "snow",                              # parallel
            phenotype = phenotype                           # Argument for estim_h2 function
        )

        # 95% confidence intervals
        boot_ci <- boot.ci(
            boot.out = bootstrap,
            conf = 0.95,
            norm = "norm"
        )

        # Saving relevant info
        pheno_res <- data.frame(
            phenotype = name,
            cases = h2_parameters$cases,
            controls = h2_parameters$controls,
            K = h2_parameters$K,
            h2_liab = bootstrap$t0,
            h2_se = sd(bootstrap$t),
            ci_95_lower = boot_ci$normal[2],
            ci_95_upper = boot_ci$normal[3]
        )
    }
    # Concatenating
    res <- rbind(res, pheno_res)
}

data <- readRDS("")
h2_res <- multi_h2(
    twin_data = data,
    phenotype_list = colnames(data)[7, length(colnames(data))],  # The first 6 columns contain other info (date of birth, sex, etc.)
    bootstrap_iter = 500
)