
#'
#' @name widen_df
#'
#' @param data a twin df where each row is an individual with a link to the twin pnr
#'
#' @returns a twin df where each row is a twin pair for sampling in the bootstrapping
#'
#' @author JAH
#' @date 2023
#'

widen_df <- function(data){
    twin_1 <- data %>%
        mutate(rows_odd = seq_len(nrow(.)) %% 2) %>%
        filter(rows_odd == 1) %>%
        select(!rows_odd)
    
    twin_2 <- data %>%
        mutate(rows_odd = seq_len(nrow(.)) %% 2) %>%
        filter(rows_odd == 0) %>%
        select(!rows_odd)
    
    colnames(twin_2) <- paste(colnames(twin_2), "2", sep = "_")
    paired_twins <- bind_cols(twin_1, twin_2)

    return(paired_twins)
}

#'
#' @name random_effects
#'
#' @param wide_df takes a wide df created using the above function, each row is a twin pair
#'
#' @returns a twin df where eachs row is an individual - now with the random effects twin_id and extra_ss_id
#'
#' @author JAH
#' @date 2023
#'

random_effects <- function(wide_df){
    require(tidyr)

    # Creating he first random effect, twin_id, a number shared in a twin pair
    twin_id <- wide_df %>%
        mutate(twin_id = 1:nrow(.),
               twin_id2 = 1:nrow(.))
    
    # Creating the second random effect, extra_ss_id, which depends on the sex configuration of the pair
    # If SS, variable is similar to twin_id
    ss_twin <- twin_id %>%
        filter(sex_config == "SS") %>%
        mutate(extra_ss_id = 1:nrow(.),
               extra_ss_id2 = 1:nrow(.)) %>%
        pivot_longer()

    # If OS, variable is different for each ind in a twin pair 80000X for one and 160000X for the other
    os_twins <- twin_id %>%
        filter(sex_config == "OS") %>%
        mutate(extra_ss_id = 1:nrow(.) + 800000,
               extra_ss_id2 = 1:nrow(.) + 1600000) %>%
        pivot_longer()
    
    # Binding all inds back together
    all_twins <- bind_rows(ss_twin, os_twins) %>%
        arrange(twin_id)

}

#'
#' @name parameters
#'
#' @param data
#' @param phenotype_col
#'
#' @returns
#'
#' @author JAH
#' @date 2023
#'

parameters <- function(data, phenotype_col){

}

#'
#' @name prob_mz
#'
#' @param data
#'
#' @returns
#'
#' @author JAH
#' @date 2023
#'

prob_mz <- function(data){

}

#'
#' @name estim_h2
#'
#' @param data
#' @param idx
#' @param phenotype
#'
#' @returns
#'
#' @author JAH
#' @date 2023
#'

 <- function(data, idx, phenotype)