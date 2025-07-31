#' Function for LPR2 oprensning
#' 
#' @author Janne Auning
#' @date 2/20/2025
#' 
#' @name lpr2 adm, diag, linked_sib_df
#' @param adm 
#' @param diag
#' @param linked_sib_df
#' 

lpr2 <- function(adm, diag, linked_sib_df) {
  diag <- diag %>% rename_with(tolower) # renaming col names in case some diag dfs contain upper case names
  adm <- adm %>% rename_with(tolower)
  
  df <- adm %>%
    filter(
      pnr %in% linked_sib_df$pnr
      # hvis der skulle v?re brugt sm? bogstaver g?r den? (T?nker der n?sten altid i vores datas?t er 
      # cleaned up s? alle er konverteret til store bogstaver, men jeg plejer ofte at have det med for 
      # en sikkerhedsskyld, s? det kan du overveje til fremtidige projekter)
      & !stringr::str_starts(c_adiag, "^DZ")       # Fjerner hvis aktionsdiagnose er DZ
      & !stringr::str_starts(c_adiag, "^Y719")     # Tilsvarende ICD8
    ) %>%
    inner_join(diag, by = "recnum") %>%
    arrange(pnr, d_inddto) %>%
    filter(
      !c_diagtype == "H"
      & !c_diagmod %in% c(1, 2, 3, 4)
      & !str_starts(c_diag, "^DZ")
    ) %>%
    left_join(linked_sib_df[c("pnr", "fdato", "follow_up_end")], by = "pnr") %>%   # Joining with the twins dataset to access the end of follow up variable
    filter(d_inddto < follow_up_end) %>%              # Beholder kun de diagnoser, der er blevet stillet inden end of follow up
    select(pnr, fdato, recnum, d_inddto, c_adiag, c_diag, c_tildiag)

  return(df)
}

#' Function for LPR3 oprensning
#' 
#' @author Janne Auning
#' @date 2/20/2025
#' 
#' @name lpr3 adm, diag, linked_sib_df
#' @param adm 
#' @param diag
#' @param linked_sib_df
#' 

lpr3 <- function(adm, diag, linked_sib_df) {

  df <- adm %>% 
    filter(
      pnr %in% linked_sib_df$pnr
      & !stringr::str_starts(AKTIONSDIAGNOSE, "^DZ")       # Fjerner hvis aktionsdiagnose er DZ
    ) %>%
    inner_join(diag, by = "DW_EK_KONTAKT") %>%
    arrange(pnr, DATO_START) %>%
    filter(
      !DIAGNOSETYPE == "H"
      # & !c_diagmod %in% c(1, 2, 3, 4) # Not in LPR3 as it stems from the ICD8 system
      & !SENERE_AFKRAEFTET == "JA"    
      & !str_starts(DIAGNOSEKODE, "^DZ")
    ) %>%
    left_join(linked_sib_df[c("pnr", "fdato", "follow_up_end")], by = "pnr") %>%
    filter(DATO_START < follow_up_end) %>% 
    select(pnr, fdato, DW_EK_KONTAKT, DATO_START, AKTIONSDIAGNOSE, DIAGNOSEKODE) %>% # No till?gsdiagnoser?
    mutate(c_tildiag = NA) %>% 
    rename(recnum = DW_EK_KONTAKT, d_inddto = DATO_START, c_adiag = AKTIONSDIAGNOSE, c_diag = DIAGNOSEKODE)
  
  return(df)
}


#' Function for pivoting lpr df into long format
#' 
#' @author Janne Auning
#' @date 02/20/2024
#' 
#' @name translate_icd8
#' @param lpr lpr df
#' @param cols_in_lpr list of columns to convert (typically c_adiag, c_diag, c_tildiag)


pivot_longer_lpr <- function(lpr, new_col, cols) {
  long_lpr <- lpr %>% 
    pivot_longer(cols = all_of(cols), 
                 names_to = "diagnosis_type", 
                 values_to = new_col) %>% 
    filter(!is.na(!!as.name(new_col))) %>% 
    distinct(pnr, !!as.name(new_col), .keep_all = TRUE) # This removes any duplicate diagnoses
  
  return(long_lpr)
}


#' Function to translate ICD8 codes in a specified column
#' 
#' @name translate_icd8
#' @param lpr lpr df
#' @param cols_in_lpr list of columns to convert (typically c_adiag, c_diag, c_tildiag)
#' @param map map to convert ICD8 to ICD10, using the one from Brunak et al
#' 

translate_icd8 <- function(lpr, col, new_col, map){
  
  lpr <- lpr %>% 
    left_join(map, by = setNames(nm = col, object = "ICD8_CODE"), relationship = "many-to-many") %>%               # Left joining with the map
    mutate(
      !!col :=                                                                      # Making some small edits to the original col
        ifelse(stringr::str_starts(!!as.name(col), "^Y|^E"), NA, !!as.name(col)),   # If Y or E diagnosis remove (as these are not in the map), otherwise keep
      !!new_col :=                                                             # Creating the new column with only ICD10
        ifelse(is.na(ICD10_CODE), !!as.name(col),                                   # If ICD10_CODE col is NA, then assume the diagnosis code is an ICD10 code
               ifelse(
                 !is.na(ICD10_CODE) & d_inddto > START_DATE & d_inddto < END_DATE,  # If ICD10_CODE col is not NA, assume the diagnosis code is an ICD8 code, and make sure d_inddto is within the time range Brunak et al have specified
                 ICD10_CODE, Inf                                                    # Otherwise, Inf - this would be a mapping of an ICD8 code to a code in the map, which isn't within the time range Brunak et al have specified
               )
        )
    ) %>%
    # When a mapping of an ICD8 code to a code in the map, which isn't within the time range Brunak et al have specified
    # the row is removed
    # 
    # This does however seem to sometimes remove some diagnoses completely, as an ICD8 has been mapped to one or several ICD10 codes 
    # in the map, but d_inddto is not within any of the start and end dates in the map. As it concerns very few diagnoses, these are simply removed
    filter(!grepl("Inf", !!sym(new_col)), !is.na(!!as.name(col))) %>%   
    
    # Are there still ICD8 diagnoses which are mapped to several ICD10 codes? I checked this by seeing if there were any mappings that were identical
    # rows across the three cols recnum, col, ICD10_CODE, but this is not the case in any of my runs
    # filter(!is.na(ICD10_CODE)) %>% group_by(recnum, !!as.name(col), ICD10_CODE) %>% filter(n() > 1) %>%  ungroup()
    
    # Removing redundant cols
    select(-c(ICD10_CODE, START_DATE, END_DATE, MAPPING_SCORE, MAPPING_FLAG, ICD8_TEXT))
  
  return(lpr)
}

#' Function to insert . in ICD10, removing D, terminal letters
#'
#' @name format_icd10
#' @param lpr lpr df
#' @param cols_in_lpr list of columns to convert (typically c_adiag, c_diag, c_tildiag)
#' 

format_icd10 <- function(lpr, col){

  lpr <- lpr %>% 
    mutate(
      !!as.name(col) := str_replace(!!as.name(col), "(\\d.*?\\d)", "\\1."),         # Inserts "."
      !!as.name(col) := str_remove(!!as.name(col), "\\.$"),                         # As the above line of code also introduced "." at the very end of the codes
      !!as.name(col) := str_remove(!!as.name(col), "^D(?=[a-zA-Z])"),               # Removing the leading D's - but only if it is followed by another letter
      !!as.name(col) := str_replace_all(!!as.name(col), "(\\d+)[a-zA-Z]+", "\\1")   # Removes trailing letters
    ) %>% 
    distinct(pnr, !!as.name(col), .keep_all = TRUE)
  
  return(lpr)
}

#' Naive matching ICD10s to phecodes on several cols
#' 
#' @name naive_mapping
#' @param lpr lpr df
#' @param cols_in_lpr list of columns to convert (typically c_adiag, c_diag, c_tildiag)
#' 

naive_mapping <- function(lpr, icd10_col, phecode_col, map){

  lpr <- lpr %>% 
    left_join(map, by = setNames(nm = icd10_col, object = "ICD10")) %>% 
    select(-c(`ICD10 String`, Phenotype, `Excl. Phecodes`, `Excl. Phenotypes`)) %>% 
    rename(!!as.name(phecode_col) := PheCode) %>% 
    distinct(pnr, !!as.name(icd10_col), !!as.name(phecode_col), .keep_all = TRUE)
  
  return(lpr)

}

#' Mapping unmapped ICD10s to a broader category, when the naive matching fails, 
#' because the specific ICD10 in question is not in the map
#' 
#' This will be in situations with codes that are mainly used in DK and not in America 
#' where the map is created
#' 
#' @name broader_mapping
#' @param lpr lpr df
#' @param cols_in_lpr list of columns to convert (typically c_adiag, c_diag, c_tildiag)
#' 

broader_mapping <- function(lpr, icd10_col, phecode_col, map) {
  
  # First removing one decimal from all the unmapped ICD10 codes
  broader_mapping_lpr <- lpr %>% 
    mutate(!!as.name(icd10_col) := ifelse(
      !is.na(!!as.name(icd10_col)) & is.na(!!as.name(phecode_col)),          # The instance when there is an unmapped ICD10 code
      str_replace(!!as.name(icd10_col), "(\\d+\\.\\d)\\d+", "\\1"),          # Then make sure it only has 1 decimal
      !!as.name(icd10_col)
    )
    ) %>% 
    left_join(map[, c("ICD10", "PheCode")], by = setNames(nm = icd10_col, object = "ICD10")) %>%     # Try to map again
    mutate(!!as.name(phecode_col) := PheCode) %>%                           # Saving the phecode col as the new one
    select(!PheCode) %>% 
    
    # This maps some of the unmapped codes. Removing all decimals for the codes that are still unmapped
    mutate(!!as.name(icd10_col) := ifelse(
      !is.na(!!as.name(icd10_col)) & is.na(!!as.name(phecode_col)),         # The instance when there is an unmapped ICD10 code
      str_replace(!!as.name(icd10_col), "(\\.\\d*)", ""),                   # Removing all decimals and "."
      !!as.name(icd10_col)
    )
    ) %>% 
    left_join(map[, c("ICD10", "PheCode")], by = setNames(nm = icd10_col, object = "ICD10")) %>%     # Try to map again
    mutate(!!as.name(phecode_col) := PheCode) %>% 
    select(!PheCode) %>% 
    distinct(pnr, !!as.name(icd10_col), !!as.name(phecode_col), .keep_all = TRUE)
  
  return(broader_mapping_lpr)
}

#' Getting the one and zero decimal phecodes as well
#' 
#' @name parental_phecodes
#' @param lpr lpr df
#' @param cols list of columns to convert (typically c_adiag, c_diag, c_tildiag)
#' 

parental_phecodes <- function(lpr, phecode_col, new_col, map){

  zero_dec_code <- paste0("zero_dec_", phecode_col)
  one_dec_code <- paste0("one_dec_", phecode_col)
  
  map$PheCode <- as.character(map$PheCode)
  
  lpr <- lpr %>% 
    mutate( 
      !!as.name(zero_dec_code) := as.character(gsub("\\..*", "", !!as.name(phecode_col))),         # phecode with 0 decimals
      !!as.name(one_dec_code) := as.character(ifelse(grepl("\\.\\d{2}", !!as.name(phecode_col)),   # phecode with 1 decimal
                                            substr(!!as.name(phecode_col), 1, nchar(!!as.name(phecode_col))-1),
                                            !!as.name(phecode_col))),
      !!as.name(phecode_col) := as.character(!!as.name(phecode_col))                               # changing class
    ) %>% 
    pivot_longer(cols = ends_with("phecode"),                                                      # Pivoting df so new cols are combined with phecode col
                 names_to = "decimal_code",                                                        # and each ind has only one instance of each unqie phecode
                 values_to = new_col) %>% 
    filter(!is.na(!!as.name(new_col))) %>%
    select(!decimal_code) %>% 
    rename(phecode = !!as.name(new_col)) %>%
    left_join(map[, c("PheCode", "Phenotype")], by = c("phecode" = "PheCode")) %>%        # Combining with phenotype info
    distinct(pnr, phecode, .keep_all = TRUE) %>%  # This removes any duplicate diagnoses
    filter(!is.na(Phenotype)) # If I've created a parental code which isn't in the map
  
  
  return(lpr)
}

#' 
#' 
#' @name case_control_status
#' @param stam
#' @param lpr lpr df
#' @param phenotypes phenotype df containing name and code
#' 

# case_control_status <- function(stam, lpr, phenotypes){
#   status_df <- stam %>% 
#     left_join(lpr[, c("pnr", "phecode")], by = "pnr") %>% 
#     distinct(pnr, phecode, .keep_all = TRUE) 
#   
#   pheno_list <- phenotypes$Phenotype
#   
#   for(pheno in pheno_list){
#     pheno_col <- paste0(pheno, "_I")
#     
#     current_code <- phenotypes[phenotypes$Phenotype == pheno, "phecode"]$phecode # Extracting the current phecode
#     
#     status_df <- status_df %>% 
#       mutate(!!as.name(pheno_col) := ifelse(phecode == current_code, 1, 0),
#              !!as.name(pheno_col) := ifelse(is.na(!!as.name(pheno_col)), 0, !!as.name(pheno_col)))
#     
#   }
#   status_df <- status_df %>%
#     group_by(pnr) %>%
#     summarise(across(ends_with("_I"), max))
#   
#   output <- stam %>% 
#     left_join(status_df, by = "pnr")
#   
#   return(output)
# }

case_control_status <- function(stam, lpr, phenotypes){
  library(data.table)
  
  stam <- as.data.table(stam)
  lpr <- as.data.table(lpr)
  pheno_list <- phenotypes$Phenotype
  
  
  # Subset to relevant codes
  lpr_filtered <- unique(lpr[Phenotype %in% pheno_list, .(pnr, Phenotype)])
  lpr_filtered[, status := 1]
  
  # Wide format: id ~ phecode, binary 1/0
  wide <- dcast(lpr_filtered, pnr ~ Phenotype, value.var = "status", fill = 0)
  setnames(wide, old = pheno_list, new = paste0(pheno_list, "_I"))
  
  # Merge into stam
  out <- left_join(stam, wide, by = "pnr")
  
  # replaced any NA with 0
  for (col in pheno_list) set(out, which(is.na(out[[col]])), col, 0)
  
  return(out)
}

#' 
#' 
#' @name compute_sex_ratio
#' @param data df with sex and phenotype cols
#' @param phenotype column containing the current 0/1 phenotype status
#' 

compute_sex_ratio <- function(data, phenotype){
  female_count <- sum(data[[phenotype]][data$sex == 1], na.rm = TRUE)
  male_count <- sum(data[[phenotype]][data$sex == 0], na.rm = TRUE)
  
  # Avoiding dividing by 0
  if (female_count == 0 | male_count == 0) return(Inf)
  
  ratio <- max(female_count/male_count, male_count/female_count)  # Returning the large ratio of the two, so it always represents the larger-to-samaller count, avoiding the irrelevant values < 1
  
  return(ratio)
}

#' 
#' 
#' @name compute_sex_ratio
#' @param data df with sex and phenotype cols
#' @param phenotype column containing the current 0/1 phenotype status

remove_imbalanced_pheno <- function(data){
  phenotype_cols <- grep("_I$", colnames(data), value = TRUE)
  non_phenotype_cols <- setdiff(colnames(data), phenotype_cols)
  
  ratios <- sapply(phenotype_cols, compute_sex_ratio, data = data)
  phenotypes_keep <- phenotype_cols[ratios <= 5]
  
  data_filtered <- data %>% select(all_of(non_phenotype_cols), all_of(phenotypes_keep))
  cat("Removed", length(phenotype_cols) - length(phenotypes_keep), "phenotypes due to high sex ratio. \n")
  
  return(data_filtered)
}


#' @name dirty_polish
#' @param data df with results
#'

dirty_polish <- function(data, birth_cohort){
  
  pretty_res <- data %>% 
    mutate(cases = round(as.numeric(cases), digits = -1),
           controls = round(as.numeric(controls), digits = -1),
           h2_liab = round(as.numeric(h2_liab), 3),
           ci_95_lower = round(as.numeric(ci_95_lower), 3),
           ci_95_upper = round(as.numeric(ci_95_upper), 3),
           h2_se = round(as.numeric(h2_se), 3),
           r_ss_liab = round(as.numeric(r_ss_liab), 3),
           r_os_liab = round(as.numeric(r_os_liab), 3),
           pval = 2 * (1 - pnorm(abs(h2_liab / h2_se))) # two-sided, pnorm of z score
    ) %>% 
    arrange(desc(cases)) %>% 
    select(!c(cases, controls))
  
  if( "K" %in% names(pretty_res)) {pretty_res <- pretty_res %>% select(-K)}
  
  
  return(pretty_res)
}
