# Libraries Used in Project

## Python Version
- Python 3.7+

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
| **Machine Learning** | `sklearn.model_selection` | `from sklearn.model_selection import train_test_split` | Train/test splitting |
| | `sklearn.ensemble` | `from sklearn.ensemble import RandomForestRegressor` | Random Forest model |
| | `sklearn.linear_model` | `from sklearn.linear_model import LinearRegression` | Linear Regression model |
| | `sklearn.metrics` | `from sklearn.metrics import mean_squared_error` | Model evaluation |
| | `sklearn.feature_selection` | `from sklearn.feature_selection import VarianceThreshold` | Feature selection |
| | `lazypredict` | `from lazypredict.Supervised import LazyRegressor` | Automated model benchmarking |
| **Environment** | `google.colab` | `from google.colab import drive` | Google Colab drive mounting |

## Summary by Package

| Package | Submodules Used | Count |
|---------|-----------------|-------|
| scikit-learn | model_selection, ensemble, linear_model, metrics, feature_selection | 5 |
| rdkit | Chem, Descriptors, Lipinski | 3 |
| Visualization | seaborn, matplotlib, plotly | 3 |
| Data | pandas, numpy | 2 |
| chembl_webresource_client | new_client | 1 |
| scipy | stats | 1 |
| lazypredict | Supervised | 1 |

## Installation

```bash
pip install pandas numpy rdkit scikit-learn lazypredict seaborn matplotlib plotly scipy chembl_webresource_client openpyxl
```

## Notes

- `google.colab` is only needed when running on Google Colab (pre-installed there)
- `openpyxl` is required by pandas for reading Excel files
- RDKit installation: Use `pip install rdkit` (not conda for macOS)
- `lazypredict` runs 30+ ML algorithms automatically for benchmarking
- PaDEL-Descriptor (external tool, not Python library) used for PubChem fingerprint generation
