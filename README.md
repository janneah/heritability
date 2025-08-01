# Heritability of Disorders in the Danish Health Registers

Project for estimating heritability of all diagnoses with sufficient prevalence in the Danish health registers using twins and siblings.

The project is divided into three parts:

1. **Simulation of twin data**
2. **Definition of phenotypes** through conversion of ICD-8 and ICD-10 to phecodes
3. **Heritability estimation** using the method from [Lakhani et al. (2019)](https://www.nature.com/articles/s41588-018-0313-7)


## 1. System Requirements

* **Dependencies**:

  * R ≥ 4.2
  * Required R packages:
    `dplyr`, `tidyr`, `purrr`, `boot`, `stringr`, `data.table`,
    `here`, `progressr`, `doFuture`, `foreach`, `future`, `MASS`
* **Hardware**: No special hardware required, but multi-core CPU recommended


## 2. Demo

The full pipeline can be run with the script `code/workflow.R`.

**Runtime** (approx.):

* Simulated data generation: <10 sec
* Mapping and preprocessing: \~30–60 sec
* Heritability estimation: \~5–10 min (25 cores, 500 bootstrap samples)


## 3. Output Files

| File                                                | Description                                       |
| --------------------------------------------------- | ------------------------------------------------- |
| `data/simulated_twins.rds `                         | Twin cohort with age, sex, zygosity proxy         |
| `data/simulated_dx.rds`                             | Diagnosis list for simulated cases                |
| `steps/twins_dx_status_400min.rds`                  | Phenotype matrix of binary traits (min 400 cases) |
| `results/h2_estimates_500_iter.rds`                 | Bootstrapped heritability estimates               |

