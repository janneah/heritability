# Requires: dplyr, tidyr, stringr, data.table
# This script assumes required packages are loaded externally or functions are referenced using :: notation.


#' Translate ICD-8 codes to ICD-10 using a mapping table
#'
#' @param lpr A data frame containing registry diagnoses
#' @param col The name of the column in `lpr` with ICD-8 codes
#' @param new_col Name of the new column to hold ICD-10 codes
#' @param map Mapping table with columns ICD8_CODE, ICD10_CODE, START_DATE, END_DATE
#'
#' @return Modified `lpr` data frame with mapped codes and unmapped/invalid codes removed

translate_icd8 <- function(lpr, col, new_col, map) {
  lpr <- lpr %>%
    dplyr::left_join(map, by = setNames(nm = col, object = "ICD8_CODE"), relationship = "many-to-many") %>%
    dplyr::mutate(
      !!col := ifelse(stringr::str_starts(!!as.name(col), "^Y|^E"), NA, !!as.name(col)),
      !!new_col := ifelse(is.na(ICD10_CODE), !!as.name(col),
                          ifelse(!is.na(ICD10_CODE) & dx_date > START_DATE & dx_date < END_DATE, ICD10_CODE, Inf))
    ) %>%
    dplyr::filter(!grepl("Inf", !!rlang::sym(new_col)), !is.na(!!as.name(col))) %>%
    dplyr::select(-c(ICD10_CODE, START_DATE, END_DATE, MAPPING_SCORE, MAPPING_FLAG, ICD8_TEXT))

  return(lpr)
}


#' Format ICD-10 codes: insert dot, remove 'D', remove trailing letters
#'
#' @param lpr Data frame with ICD-10 codes
#' @param col Column name to format
#'
#' @return Modified data frame with formatted ICD-10 codes

format_icd10 <- function(lpr, col) {
  lpr <- lpr %>%
    dplyr::mutate(
      !!as.name(col) := stringr::str_replace(!!as.name(col), "(\\d.*?\\d)", "\\1."),
      !!as.name(col) := stringr::str_remove(!!as.name(col), "\\.$"),
      !!as.name(col) := stringr::str_remove(!!as.name(col), "^D(?=[a-zA-Z])"),
      !!as.name(col) := stringr::str_replace_all(!!as.name(col), "(\\d+)[a-zA-Z]+", "\\1")
    ) %>%
    dplyr::distinct(id, !!as.name(col), .keep_all = TRUE)

  return(lpr)
}


#' Naively map ICD-10 codes to phecodes
#'
#' @param lpr Data frame with ICD-10 codes
#' @param icd10_col Column name with ICD-10 codes
#' @param phecode_col Desired output column name for phecodes
#' @param map Mapping table with ICD10 and PheCode columns
#'
#' @return Data frame with matched phecodes

naive_mapping <- function(lpr, icd10_col, phecode_col, map) {
  lpr <- lpr %>%
    dplyr::left_join(map, by = setNames(nm = icd10_col, object = "ICD10")) %>%
    dplyr::select(-c(`ICD10 String`, Phenotype, `Excl. Phecodes`, `Excl. Phenotypes`)) %>%
    dplyr::rename(!!as.name(phecode_col) := PheCode) %>%
    dplyr::distinct(id, !!as.name(icd10_col), !!as.name(phecode_col), .keep_all = TRUE)

  return(lpr)
}


#' Broader mapping: fallback when ICD-10 code doesn't match directly
#'
#' @param lpr Data frame with ICD-10 codes
#' @param icd10_col ICD-10 column name
#' @param phecode_col Output phecode column name
#' @param map Mapping table with ICD10 and PheCode columns
#'
#' @return Data frame with additional mapped phecodes using generalised codes

broader_mapping <- function(lpr, icd10_col, phecode_col, map) {
  broader_mapping_lpr <- lpr %>%
    dplyr::mutate(!!as.name(icd10_col) := ifelse(
      !is.na(!!as.name(icd10_col)) & is.na(!!as.name(phecode_col)),
      stringr::str_replace(!!as.name(icd10_col), "(\\d+\\.\\d)\\d+", "\\1"),
      !!as.name(icd10_col)
    )) %>%
    dplyr::left_join(map[, c("ICD10", "PheCode")], by = setNames(nm = icd10_col, object = "ICD10")) %>%
    dplyr::mutate(!!as.name(phecode_col) := PheCode) %>%
    dplyr::select(!PheCode) %>%
    dplyr::mutate(!!as.name(icd10_col) := ifelse(
      !is.na(!!as.name(icd10_col)) & is.na(!!as.name(phecode_col)),
      stringr::str_replace(!!as.name(icd10_col), "(\\.\\d*)", ""),
      !!as.name(icd10_col)
    )) %>%
    dplyr::left_join(map[, c("ICD10", "PheCode")], by = setNames(nm = icd10_col, object = "ICD10")) %>%
    dplyr::mutate(!!as.name(phecode_col) := PheCode) %>%
    dplyr::select(!PheCode) %>%
    dplyr::distinct(id, !!as.name(icd10_col), !!as.name(phecode_col), .keep_all = TRUE)

  return(broader_mapping_lpr)
}


#' Generate parental phecodes (zero or one decimal)
#'
#' @param lpr Data frame with mapped phecodes
#' @param phecode_col Name of the original phecode column
#' @param new_col Desired output column for merged parental phecodes
#' @param map Phecode mapping table (with Phenotype names)
#'
#' @return Long-form data frame with 0- and 1-decimal phecodes merged

parental_phecodes <- function(lpr, phecode_col, new_col, map) {
  zero_dec_code <- paste0("zero_dec_", phecode_col)
  one_dec_code <- paste0("one_dec_", phecode_col)

  map$PheCode <- as.character(map$PheCode)

  lpr <- lpr %>%
    dplyr::mutate(
      !!zero_dec_code := gsub("\\..*", "", !!as.name(phecode_col)),
      !!one_dec_code := ifelse(
        grepl("\\.\\d{2}", !!as.name(phecode_col)),
        substr(!!as.name(phecode_col), 1, nchar(!!as.name(phecode_col)) - 1),
        !!as.name(phecode_col)
      ),
      !!as.name(phecode_col) := as.character(!!as.name(phecode_col))
    ) %>%
    tidyr::pivot_longer(cols = dplyr::ends_with("phecode"),
                        names_to = "decimal_code",
                        values_to = new_col) %>%
    dplyr::filter(!is.na(!!as.name(new_col))) %>%
    dplyr::select(!decimal_code) %>%
    dplyr::rename(phecode = !!as.name(new_col)) %>%
    dplyr::left_join(map[, c("PheCode", "Phenotype")], by = c("phecode" = "PheCode"), relationship = "many-to-many") %>%
    dplyr::distinct(id, phecode, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(Phenotype))

  return(lpr)
}


#' Generate binary case/control matrix for phenotypes
#'
#' @param stam Base cohort data with id identifiers
#' @param lpr Long-format df with mapped phenotypes per id
#' @param phenotypes Data frame with phenotype names and phecodes
#'
#' @return Wide-format case/control matrix merged with `stam`

case_control_status <- function(stam, lpr, phenotypes){
    status_df <- stam %>%
      left_join(lpr[, c("id", "phecode")], by = "id") %>%
      distinct(id, phecode, .keep_all = TRUE)

    pheno_list <- phenotypes$Phenotype

    for(pheno in pheno_list){
      pheno_col <- paste0(pheno, "_I")

      current_code <- phenotypes[phenotypes$Phenotype == pheno, "phecode"]$phecode # Extracting the current phecode

      status_df <- status_df %>%
        mutate(!!as.name(pheno_col) := ifelse(phecode == current_code, 1, 0),
               !!as.name(pheno_col) := ifelse(is.na(!!as.name(pheno_col)), 0, !!as.name(pheno_col)))

    }
    status_df <- status_df %>%
      group_by(id) %>%
      summarise(across(ends_with("_I"), max))

    output <- stam %>%
      left_join(status_df, by = "id")

    return(output)
  }



#' Compute the sex ratio (larger/smaller) for a binary phenotype
#'
#' @param data Data frame with `sex` and phenotype columns
#' @param phenotype Character name of phenotype column (0/1)
#'
#' @return Numeric sex ratio or Inf if either sex count is 0

compute_sex_ratio <- function(data, phenotype) {
  female_count <- sum(data[[phenotype]][data$sex == 1], na.rm = TRUE)
  male_count <- sum(data[[phenotype]][data$sex == 0], na.rm = TRUE)

  if (female_count == 0 | male_count == 0) return(Inf)

  ratio <- max(female_count / male_count, male_count / female_count)
  return(ratio)
}


#' Remove phenotypes with sex ratio > 5
#'
#' @param data Case/control matrix with phenotype columns (ending in _I)
#'
#' @return Data frame excluding high-imbalance phenotypes

remove_imbalanced_pheno <- function(data) {
  phenotype_cols <- grep("_I$", colnames(data), value = TRUE)
  non_phenotype_cols <- setdiff(colnames(data), phenotype_cols)

  ratios <- sapply(phenotype_cols, compute_sex_ratio, data = data)
  phenotypes_keep <- phenotype_cols[ratios <= 5]

  data_filtered <- dplyr::select(data, dplyr::all_of(non_phenotype_cols), dplyr::all_of(phenotypes_keep))
  cat("Removed", length(phenotype_cols) - length(phenotypes_keep), "phenotypes due to high sex ratio. \n")

  return(data_filtered)
}

