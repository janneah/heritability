#' Converting ICD-8 to ICD-10 
#' 
#' @author Janne Auning Hansen
#' @date 18-01-2024
#' 
#' @description
#' Script for converting ICD-8 codes to ICD-10 using the mapping from https://pubmed.ncbi.nlm.nih.gov/37555907/
#' 
#' @Todo
#' 

devtools::install_github("PheWAS/PheWAS")
library(PheWAS)
library(dplyr)
library(data.table)

mapping <- fread("data/10654_2023_1027_MOESM1_ESM.txt")
icd8 <- fread("data/icd8.csv")

# Creating dummy data --------------------------------------------------------

dummy <- generateExample(n = 10000, phenotypes.per = 10)

inds <- dummy$id.sex %>% 
  mutate(diag = sample(x = icd8$Code, size = 10000, replace = TRUE))

# Converting ICD code -------------------------------------------------------








