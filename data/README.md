# 📁 data/

This directory contains project-specific datasets and mapping resources used for phenotype definition and simulation.

## Subfolders and Files

### 📁 `data/maps/`

Contains code and phenotype mapping tables used for harmonisation and analysis:

- `icd8_to_icd10.csv` — Mapping table for converting ICD-8 codes to ICD-10.
- `icd10_to_phecodes.csv` — Mapping table for converting ICD-10 codes to Phecodes.

These mappings are used to translate registry-based diagnostic codes into phenotypic groupings suitable for downstream genetic analyses.

### 📄 `data/simulated_twins.csv`

Simulated twin dataset generated using the script `01a_simulate_twins_data.R`.  
This dataset is used for heritability estimation and model validation.

