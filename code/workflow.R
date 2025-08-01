# Main script for running the full analysis pipeline:
# 1. Simulate data
# 2. Preprocess data (map ICD codes to phecodes)
# 3. Estimate heritability

# ---------------------------
# 📦 Load required packages
# ---------------------------
required_packages <- c("tidyverse", "data.table", "here", "stringr")  # Add others as needed
lapply(required_packages, require, character.only = TRUE)

# Use project root as working directory
setwd(here::here())

# ---------------------------
# 🔧 Load helper functions
# ---------------------------
source("code/00_helpers/h2_helpers.R")
source("code/00_helpers/mapping_helpers.R")

# ---------------------------
# 1️⃣ Simulate twin dataset ~ < 1 min
# ---------------------------
source("code/01_preprocessing/01a_simulate_twins_df.R")
# Writes output to: data/simulated_twins_df.rds and data/simulated_dx.rds

# ---------------------------
# 2️⃣ Map ICD codes to phecodes ~ 1 min
# ---------------------------
source("code/01_preprocessing/01b_map_phecodes.R")
# Writes output to steps/twins_dx_status_min400.rds

# ---------------------------
# 3️⃣ Estimate heritability ~ 5-10 min
# ---------------------------
source("code/02_analysis/02_estimate_h2.R")
# Output: results/h2_estimates.csv

# ---------------------------
# ✅ Done
# ---------------------------
message("✅ Pipeline complete.")
