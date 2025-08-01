#' @title ICD8 → ICD10 → Phecode Mapping Pipeline
#' @description Loads diagnosis data, translates ICD-8 codes to ICD-10, then maps to phecodes.
#' @details This script is designed for Danish registry data with simulated twin/sibling structure.

# ------------------------------------------------------------------
# Load data --------------------------------------------------------
# ------------------------------------------------------------------

stam  <- readRDS(here::here("data", "simulated_twin_df.rds"))
lpr <- readRDS(here::here("data", "simulated_dx.rds"))

# ------------------------------------------------------------------
# Load helpers and mappings ----------------------------------------
# ------------------------------------------------------------------

source(here::here("code", "00_helpers", "mapping_helpers.R"))

# ICD8 → ICD10 mapping
icd8_map <- data.table::fread(here::here("data", "maps", "map_icd8_2_icd10.txt"), sep = "\t") |>
  dplyr::mutate(
    START_DATE  = as.Date(as.character(START_DATE), format = "%Y%m%d"),
    END_DATE    = as.Date(as.character(END_DATE), format = "%Y%m%d"),
    ICD8_CODE   = as.character(ICD8_CODE),
    ICD10_CODE  = dplyr::na_if(ICD10_CODE, "")
  ) |>
  dplyr::filter(!is.na(ICD10_CODE))

# Phecode map and phenotype info
load(here::here("data", "maps", "pheinfo.rda"))
load(here::here("data", "maps", "sex_restriction.RData"))
phecode_map <- fread(here::here("data", "maps", "phecode_icd10.csv"))

# ------------------------------------------------------------------
# ICD8 → ICD10 -----------------------------------------------------
# ------------------------------------------------------------------
icd_col <- "icd"
icd10_col <- "icd10"

lpr_icd10 <- translate_icd8(
  lpr      = lpr,
  col      = icd_col,
  new_col  = icd10_col,
  map      = icd8_map
)

# ------------------------------------------------------------------
# Reformat ICD10 ---------------------------------------------------
# ------------------------------------------------------------------

lpr_reformatted_icd10 <- format_icd10(
  lpr = lpr_icd10,
  col = icd10_col
)

# ------------------------------------------------------------------
# ICD10 → Phecodes -------------------------------------------------
# ------------------------------------------------------------------

phecode_col <- "phecode"

# Step 1: Naive mapping
lpr_phecodes <- naive_mapping(
  lpr          = lpr_reformatted_icd10,
  icd10_col    = icd10_col,
  phecode_col  = phecode_col,
  map          = phecode_map
)

# Step 2: Broader mapping (e.g. F30.9 → F30)
lpr_phecodes <- broader_mapping(
  lpr          = lpr_phecodes,
  icd10_col    = icd10_col,
  phecode_col  = phecode_col,
  map          = phecode_map
)

# Step 3: Add all parental phecodes (e.g. 803.12 → 803.1 and 803)
lpr_parental_phecodes <- parental_phecodes(
  lpr          = lpr_phecodes,
  phecode_col  = phecode_col,
  new_col      = "all_phecodes",
  map          = phecode_map
)

# ------------------------------------------------------------------
# Build phenotype matrix -------------------------------------------
# ------------------------------------------------------------------

# Filter sex-restricted phenotypes
sex_restriction <- dplyr::left_join(
  sex_restriction,
  pheinfo[, c("phecode", "description")],
  by = "phecode"
) |>
  dplyr::filter(male_only == TRUE | female_only == TRUE)

# Define phenotypes with ≥ 400 cases (only any category which the epilepsy icds might fall into in this case)
pheno_400 <- lpr_parental_phecodes |>
  dplyr::group_by(Phenotype, phecode) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n)) |>
  dplyr::filter(!Phenotype %in% sex_restriction$description) |>
  dplyr::filter(n >= 400)

head(pheno_400)

# ------------------------------------------------------------------
# Case-control and filtering ---------------------------------------
# ------------------------------------------------------------------

dx_status <- case_control_status(
  stam       = stam,
  lpr        = lpr_parental_phecodes,
  phenotypes = pheno_400
)

dx_status_flt <- remove_imbalanced_pheno(dx_status)

# ------------------------------------------------------------------
# Save results -----------------------------------------------------
# ------------------------------------------------------------------

saveRDS(dx_status_flt, here::here("steps", "twins_dx_status_400min.rds"))
