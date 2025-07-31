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
# 1️⃣ Simulate twin dataset
# ---------------------------
source("code/01_preprocessing/01a_simulate_twins_df.R")
# Assumes it writes output to: data/simulated_twins.csv

# ---------------------------
# 2️⃣ Map ICD codes to phecodes
# ---------------------------
source("code/01_preprocessing/01b_map_phecodes.R")
# Assumes it creates or modifies: data/processed_phenotypes.csv

# ---------------------------
# 3️⃣ Estimate heritability
# ---------------------------
source("code/02_analysis/02_estimate_h2.R")
# Output: results/h2_estimates.csv

# ---------------------------
# ✅ Done
# ---------------------------
message("✅ Pipeline complete.")
