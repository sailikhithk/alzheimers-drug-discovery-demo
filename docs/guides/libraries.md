# Libraries Used in Project

## Python Version
- Python 3.7+ (tested on Python 3.12 in Google Colab)

## Complete Library Reference

| Category | Library | Import Statement | Purpose |
|----------|---------|------------------|---------|
| **Core Python** | `os` | `import os` | File system operations |
| **Data Processing** | `numpy` | `import numpy as np` | Numerical operations |
| | `pandas` | `import pandas as pd` | Data manipulation, CSV/Excel I/O |
| **Chemistry** | `rdkit.Chem` | `from rdkit import Chem` | Molecular structure parsing |
| | `rdkit.Chem.Descriptors` | `from rdkit.Chem import Descriptors` | Molecular descriptor calculation |
| | `rdkit.Chem.Lipinski` | `from rdkit.Chem import Lipinski` | Lipinski's Rule of Five |
| | `chembl_webresource_client` | `from chembl_webresource_client.new_client import new_client` | ChEMBL database API access |
| **Visualization** | `seaborn` | `import seaborn as sns` | Statistical visualizations |
| | `matplotlib.pyplot` | `import matplotlib.pyplot as plt` | Plotting |
| | `plotly.express` | `import plotly.express as px` | Interactive plots |
| **Statistics** | `scipy.stats` | `from scipy.stats import mannwhitneyu` | Mann-Whitney U test |
| **Machine Learning** | `sklearn.model_selection` | `from sklearn.model_selection import train_test_split, cross_val_score, KFold, GridSearchCV` | Train/test splitting, cross-validation, hyperparameter tuning |
| | `sklearn.ensemble` | `from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor, GradientBoostingRegressor, VotingRegressor` | Ensemble models |
| | `sklearn.svm` | `from sklearn.svm import SVR` | Support Vector Regression |
| | `sklearn.metrics` | `from sklearn.metrics import mean_squared_error, r2_score` | Model evaluation |
| | `sklearn.feature_selection` | `from sklearn.feature_selection import VarianceThreshold` | Feature selection |
| | `lazypredict` | `from lazypredict.Supervised import LazyRegressor` | Automated model benchmarking |
| | `xgboost` | `from xgboost import XGBRegressor` | XGBoost model |
| | `lightgbm` | `from lightgbm import LGBMRegressor` | LightGBM model |
| **Environment** | `google.colab` | `from google.colab import drive` | Google Colab drive mounting |

## Summary by Package

| Package | Submodules Used | Count |
|---------|-----------------|-------|
| scikit-learn | model_selection, ensemble, svm, metrics, feature_selection | 5 |
| rdkit | Chem, Descriptors, Lipinski | 3 |
| Visualization | seaborn, matplotlib, plotly | 3 |
| Data | pandas, numpy | 2 |
| Boosting | xgboost, lightgbm | 2 |
| chembl_webresource_client | new_client | 1 |
| scipy | stats | 1 |
| lazypredict | Supervised | 1 |

## Installation (Recommended: uv)

Using `uv` for faster, smarter dependency resolution:

```bash
# Install uv
pip install uv

# Install all dependencies with compatible versions
uv pip install --system \
    scikit-learn numpy pandas matplotlib seaborn \
    lazypredict chembl_webresource_client rdkit \
    openpyxl padelpy xgboost lightgbm
```

### Alternative: Standard pip

```bash
pip install pandas numpy rdkit scikit-learn lazypredict seaborn matplotlib plotly scipy chembl_webresource_client openpyxl xgboost lightgbm
```

## Notes

- **uv** is recommended for Google Colab to avoid numpy/scikit-learn version conflicts
- `google.colab` is only needed when running on Google Colab (pre-installed there)
- `openpyxl` is required by pandas for reading Excel files
- RDKit installation: Use `pip install rdkit` (not conda for macOS)
- `lazypredict` runs 30+ ML algorithms automatically for benchmarking
- PaDEL-Descriptor (external tool, not Python library) used for PubChem fingerprint generation
- `xgboost` and `lightgbm` are used by LazyPredict and for model comparison
