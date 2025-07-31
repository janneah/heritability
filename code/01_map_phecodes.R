#' Script for translating ICD codes into phecodes.
#'
#' @author Janne AH
#' @date 01-03-2024
#' 
#' @description In the code, any ICD8 codes in the cohort of interest are first translated to ICD10, and then all ICD10 codes are translated into phecodes
#' 


# Libs -----

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

# Loading in data ---------------------------------------------------------
lpr <- readRDS("simulated_twin_lpr.rds")) 


# Loading in aux code and maps ----------------------------------------------
source("code/helper_functions/from_icd_to_phecodes.R")

# Map to translate ICD8 to ICD10
icd8_map <- fread("data/maps/mapICD8ToICD10.txt", sep = "\t") %>%
  # Doing some basic reformatting of some of the cols
  mutate(START_DATE = as.Date(as.character(START_DATE), format = "%Y%m%d"),
         END_DATE = as.Date(as.character(END_DATE), format = "%Y%m%d"),
         ICD8_CODE = as.character(ICD8_CODE),
         ICD10_CODE = ifelse(ICD10_CODE == "", NA, ICD10_CODE)) %>% 
  filter(!is.na(ICD10_CODE))                                # Some ICD8s are in the map but without any ICD10, these are removed

# Phecode map and metainfo
load("data/maps/pheinfo.rda")
load("data/maps/sex_restriction.RData")
phecode_map <- readRDS("data/maps/phecode_icd10.rds")      # This .rds version is slightly edited from the original one

# Pivot longer ----------------------------------------------------------------
# Going from having three columns with diagnoses to one. Removing NAs and duplicate diagnoses
icd_col <- "icd"

long_lpr <- pivot_longer_lpr(lpr,
                             new_col = icd_col,
                             cols = c("c_diag", "c_adiag", "c_tildiag"))

head(long_lpr)

# ICD8 to ICD10 ---------------------------------------------------------------
icd10_col <- "icd10"
# debugonce(translate_icd8) # For debugging
lpr_icd10 <- translate_icd8(lpr = long_lpr, 
                            col = icd_col,
                            new_col = icd10_col,
                            map = icd8_map)

head(lpr_icd10)

# ICD10 to phecodes -----------------------------------------------------------

# Removing trailing letters in the ICD10s (mainly a Danish convention)
# E.g. R10.2C becomes R10.2

# Reformatting the ICD10s - inserting . and removing trailing letters

##### Consider moving this part to before ICD8-ICD10 as there is a mixture of ICD10 codes with different format (see Betinas comment)
##### Consider removing "." from map rather than inserting it in dataset

lpr_reformatted_icd10 <- format_icd10(lpr = lpr_icd10,
                                      col = icd10_col)

head(lpr_reformatted_icd10)

# Translating the formatted ICD codes to phecodes
# This is done in the following steps:

# 1. First naive mapping of ICD10 codes to phecodes:
phecode_col <- "phecode"
# debugonce(naive_mapping)
lpr_phecodes <- naive_mapping(lpr = lpr_reformatted_icd10, 
                              icd10_col = icd10_col,
                              phecode_col = phecode_col,
                              map = phecode_map)

# 2. Some of the ICD10 do not map directly, however. Possibly due to some specific ICD10 codes
# being used mainly in the Danish system and not in the American. Such codes (as e.g. many of
# the .9 codes, are simply mapped to a parental code, e.g., F30.9 being mapped to F30)

lpr_phecodes <- broader_mapping(lpr = lpr_phecodes,
                                icd10_col = icd10_col,
                                phecode_col = phecode_col,
                                map = phecode_map)

head(lpr_phecodes)


# 3. Now, all diagnoses should be counted in all levels of the phecode, so if an individual 
# has the phecode 803.12, it should also be mapped to 803.1 and 803 
lpr_parental_phecodes <- parental_phecodes(lpr = lpr_phecodes,
                                           phecode_col = phecode_col,
                                           new_col = "all_phecodes",
                                           map = phecode_map)
  
head(lpr_parental_phecodes)

# Now creating the df with twin/sibling info as well as case/control status for all phecodes with at least 400 cases (0.5%) -----

# Filtering
sex_restriction <- sex_restriction %>% 
  left_join(pheinfo[, c("phecode", "description")], by = "phecode") %>% 
  filter(male_only == TRUE | female_only == TRUE)



# Making case control matrix ------------------------------------------------------------------------------
if(pair == "twins") {
  stam <- readRDS(paste0("steps/stam_lpr/linked_twins_", birth_cohort, ".rds"))
} else {
  stam <- readRDS(paste0("steps/stam_lpr/all_linked_sibs_", birth_cohort, ".rds"))
}

# Finding phenotypes with minimally 400 cases
pheno_400 <- lpr_parental_phecodes %>% 
  group_by(Phenotype, phecode) %>% 
  summarise(n = n()) %>%
  arrange(desc(n)) %>% 
  filter(!Phenotype %in% sex_restriction$description) %>% # Sex restriction def from PheWas
  filter(n >= 400)  

# Making case control status df
# debugonce(case_control_status)
dx_status <- case_control_status(stam, lpr_parental_phecodes, pheno_400) # Some issue with 1955 sibs, CHECK that inds are counted correctly

# Filtering away phenotypes with imbalanced sex ratios
dx_status_flt <- remove_imbalanced_pheno(dx_status)


saveRDS(dx_status_flt, paste0("steps/pheno/", pair, "_", birth_cohort, "_phecode_status_min400.rds"))

