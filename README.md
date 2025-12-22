# Alzheimer's Drug Discovery using Machine Learning

A computational drug discovery project using QSAR (Quantitative Structure-Activity Relationship) modeling to predict bioactivity of compounds against Beta-amyloid A4 protein — a key target in Alzheimer's disease treatment.

## Project Goal

Build machine learning models that predict whether chemical compounds will be effective inhibitors of Beta-amyloid A4 protein, enabling virtual screening of drug candidates before expensive lab testing.

## Key Results

| Metric | Value |
|--------|-------|
| Best Test R² | 0.78 (RandomForest & HistGradientBoosting) |
| Best Test RMSE | 0.58 |
| Dataset Size | 1,319 compounds |
| Features | 881 PubChem fingerprints → 178 after variance filtering |

### Top Performing Models

| Model | Train R² | Test R² | RMSE | Overfitting |
|-------|----------|---------|------|-------------|
| RandomForestRegressor | 0.95 | 0.78 | 0.58 | Moderate |
| HistGradientBoostingRegressor | 0.94 | 0.78 | 0.59 | Moderate |
| XGBRegressor | 0.97 | 0.77 | 0.61 | High |
| MLPRegressor | 0.94 | 0.75 | 0.63 | High |
| SVR | 0.86 | 0.72 | 0.67 | Low |

## Pipeline Overview

```
ChEMBL Database → Raw Data (7,918) → Curated (1,319) → Fingerprints (881) → ML Model → Predictions
```

## Quick Start

```bash
# Install dependencies
pip install pandas numpy rdkit scikit-learn lazypredict seaborn matplotlib chembl_webresource_client padelpy openpyxl

# Run notebook
jupyter notebook ml_for_alzheimer_s_drug_discovery.ipynb
```

## Documentation

- [glossary.md](glossary.md) — Domain terminology and project phases
- [index.md](index.md) — Workflow index mapping theory to implementation
- [libraries.md](libraries.md) — Complete library reference
- [linux.md](linux.md) — Shell commands used
- [notebook-analysis.md](notebook-analysis.md) — Detailed notebook evaluation
- [README_DETAILED.md](README_DETAILED.md) — Detailed pipeline documentation
- [README_COMPREHENSIVE.md](README_COMPREHENSIVE.md) — Complete guide with biology/chemistry background

## Domain Context

- **Target:** Beta-amyloid A4 protein (ChEMBL ID: CHEMBL2487)
- **Data source:** ChEMBL database
- **Measurement:** IC50 (Half Maximal Inhibitory Concentration)
- **Classification:**
  - Active: IC50 ≤ 1,000 nM (pIC50 ≥ 6)
  - Intermediate: 1,000 < IC50 < 10,000 nM
  - Inactive: IC50 ≥ 10,000 nM (pIC50 ≤ 5)
