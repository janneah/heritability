#' Converting ICD-8 to ICD-10 
#' 
#' @author Janne Auning Hansen
#' @date 18-01-2024
#' 
#' @description
#' Script for converting ICD-8 codes to ICD-10 using the mapping from https://pubmed.ncbi.nlm.nih.gov/37555907/
#' 
#' @Todo
#'    - In left_join: account for date of diagnosis
#' 

devtools::install_github("PheWAS/PheWAS")
library(PheWAS)
library(dplyr)
library(data.table)
library(stringr)

mapping <- fread("data/10654_2023_1027_MOESM1_ESM.txt", colClasses = "character") %>% 
  mutate(START_DATE = as.Date(as.character(START_DATE), format = "%Y%m%d"),               # Fixing date format
         END_DATE = as.Date(as.character(END_DATE), format = "%Y%m%d"))                  

icd8 <- fread("data/SGDklass_ICD8.txt", header = FALSE, col.names = c("ICD8_CODE", "W/E")) %>% 
  mutate(ICD8_CODE = str_replace(ICD8_CODE, "^dia", "")) %>%                             # Removing prefix
  filter(!(str_detect(ICD8_CODE, "^E") | str_detect(ICD8_CODE, "^Y")))                   # Removing E and Y diagnoses (not in map)

# Creating dummy data --------------------------------------------------------

dummy <- generateExample(n = 10000, phenotypes.per = 10)

inds <- dummy$id.sex %>% 
  mutate(icd8_diag = sample(x = icd8$ICD8_CODE, size = 10000, replace = TRUE),
         dinddto = sample(x = seq(as.Date('1965-01-01'), as.Date('1996-01-01'), by="day"), size = 10000, replace = TRUE))

# Converting ICD code -------------------------------------------------------

converted_diag <- inds %>%                                   
  left_join(mapping[, c("ICD8_CODE", "ICD10_CODE")], by = c("icd8_diag" = "ICD8_CODE"))




