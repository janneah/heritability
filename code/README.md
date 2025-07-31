# 📁 code/

This folder contains all R scripts used for data processing, simulation, and analysis in the project. The structure follows a workflow, divided into helpers, preprocessing, and analysis stages, which are all run from the master script `workflow.R`.

## Folder Structure

### 📁 `00_helpers/`

General-purpose utility functions used across scripts.

- `h2_helpers.R` — Functions for estimating heritability.
- `mapping_helpers.R` — Functions for working with ICD and phecode mappings.

### 📁 `01_preprocessing/`

Scripts for data preparation and harmonisation.

- `01a_simulate_twins_df.R` — Simulates a twin dataset for heritability estimation.
- `01b_map_phecodes.R` — Maps ICD codes to phecodes.

### 📁 `02_analysis/`

Scripts for running the linear mixed model with 500 bootstrap iterations.

- `02_estimate_h2.R` — Estimates heritability using the harmonised data.

---

## 📌 Notes

- `01_preprocessing/` does not include the full preprocessing pipeline applied to the data extracted from the Danish national registries. Only components relevant for simulation or reproducible subsets (e.g. phecode mapping) are provided.
- Scripts are numbered to reflect the recommended execution order.
- Shared functions are modularised under `00_helpers`.
- Data files referenced by these scripts are located in the `data/` or `steps/` folder.

