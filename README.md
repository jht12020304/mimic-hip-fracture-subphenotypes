# Emergency Pre-Admission Host Physiological Subphenotypes in Older Adults with Acute Hip Fracture

Official analysis code, closed-form 9-biomarker subphenotype posterior probability calculator, and aggregate summary tables accompanying the study:

> **Emergency Pre-Admission Host Physiological Subphenotypes and Prognostic Stratification in Older Adults with Acute Hip Fracture**

---

## Overview

Using first-draw emergency department (ED) complete blood count (CBC) and basic metabolic panel (BMP) biomarkers obtained at or before hospital admission (`[admittime - 24 h, admittime]`) across **2,488 older adults (age ≥ 65 years)** with acute hip fracture (and **5,508 same-hospital older adults** with non-hip major fractures) in **MIMIC-IV v3.1 (2008–2022)**, diagonal unequal-variance Gaussian mixture modeling (`VVI`, `K = 3`) resolves three reproducible host physiological subphenotypes while modeling baseline anticoagulation separately:

- **P1: Homeostatic (38.9%)** — Preserved hemoglobin, renal filtration, electrolytes, and acid–base balance.
- **P2: Stress–Inflammatory (35.3%)** — Leukocytosis and stress hyperglycemia with preserved baseline filtration.
- **P3: Cardiorenal–Anemic (25.8%)** — Pre-admission anemia, elevated red cell distribution width (RDW), azotemia, higher potassium, and lower bicarbonate.

---

## Repository Structure

```text
mimic-hip-fracture-subphenotypes/
├── code/
│   ├── 01_extract_preadmit_cohorts.py      # DuckDB/SQL cohort, first-draw ED lab, & incident AKI extraction
│   ├── 02_run_subphenotype_pipeline.R      # GMM (20-model grid, BLRT, VVI K=3), poLCA, Cox, Fine-Gray AKI, & CCW
│   ├── 03_build_publication_figures.R      # Publication figures (Figures 1–5 & Supplementary Figures S1–S3)
│   └── 04_posterior_calculator_9lab.R      # Standalone closed-form 9-biomarker posterior probability calculator
├── parameters/
│   └── SuppTable_S5_Exact_9Lab_Posterior_Parameters.csv  # Exact Pre-2017 winsorization, scaling & Gaussian parameters
├── summary_tables/                         # Aggregate summary tables (Tables 1–3, Supplementary Tables S1–S7)
├── LICENSE                                 # MIT License
└── .gitignore
```

---

## Quick Start: Closed-Form 9-Biomarker Subphenotype Calculator (No MIMIC-IV Required)

Any researcher or clinician can compute posterior probabilities `(Pr[P1], Pr[P2], Pr[P3])` directly from routine emergency CBC + BMP measurements without needing access to MIMIC-IV:

```bash
Rscript code/04_posterior_calculator_9lab.R
```

Or inside R:

```r
source("code/04_posterior_calculator_9lab.R")

res <- predict_hip_subphenotype_9lab(
  hb          = 10.4,   # Hemoglobin (g/dL)
  rdw         = 15.7,   # Red cell distribution width (%)
  platelets   = 204,    # Platelets (x10^3/uL)
  wbc         = 9.8,    # White blood cell count (x10^3/uL)
  creatinine  = 1.40,   # Serum creatinine (mg/dL)
  glucose     = 125,    # Blood glucose (mg/dL)
  sodium      = 138,    # Serum sodium (mEq/L)
  potassium   = 4.7,    # Serum potassium (mEq/L)
  bicarbonate = 23,     # Serum bicarbonate (mEq/L)
  param_csv   = "parameters/SuppTable_S5_Exact_9Lab_Posterior_Parameters.csv"
)
print(res)
```

---

## Reproducing the Full MIMIC-IV v3.1 Analysis Pipeline

### 1. Data Access & Environment Requirements
- **MIMIC-IV (v3.1):** Available to credentialed researchers via PhysioNet (`https://doi.org/10.13026/kpb9-mt58`) under the PhysioNet Credentialed Health Data Use Agreement. Load the core `hosp` and `icu` modules into a local DuckDB database (`mimiciv.duckdb`). *In accordance with PhysioNet governance, zero patient-level records are stored in this repository.*
- **Python (≥ 3.10):** `duckdb`
- **R (≥ 4.3):** `mclust`, `poLCA`, `lcmm`, `survival`, `cmprsk`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `scales`

### 2. Step-by-Step Execution

```bash
# Step 1: Extract pre-admission cohorts & reference-interval audit from MIMIC-IV v3.1 DuckDB
export MIMIC_DUCKDB_PATH="/path/to/mimiciv.duckdb"
python3 code/01_extract_preadmit_cohorts.py

# Step 2: Run Gaussian mixture modeling, temporal validation, competing-risk AKI, & CCW target trial emulation
Rscript code/02_run_subphenotype_pipeline.R

# Step 3: Generate all main and supplementary figures (600 dpi PNG & vector PDF/SVG)
Rscript code/03_build_publication_figures.R
```

---

## License

Released under the [MIT License](LICENSE).
