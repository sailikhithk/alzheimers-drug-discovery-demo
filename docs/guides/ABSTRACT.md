# Machine Learning-Based QSAR Modeling for Predicting Bioactivity of Beta-Amyloid A4 Protein Inhibitors: A Computational Approach to Alzheimer's Drug Discovery

**Author:** Sai Likhith Kanuparthi

---

## Abstract

Alzheimer's disease is a progressive neurodegenerative disorder affecting over 55 million people worldwide, with beta-amyloid protein aggregation representing a primary therapeutic target. Traditional drug discovery requires over a decade and billions of dollars per approved compound, necessitating computational approaches to accelerate candidate identification. This study develops Quantitative Structure-Activity Relationship models using machine learning to predict the bioactivity of small molecule inhibitors against beta-amyloid A4 protein. Bioactivity data for 1,319 compounds were retrieved from the ChEMBL database, filtered for half maximal inhibitory concentration measurements in nanomolar units, and classified as active, intermediate, or inactive based on standard potency thresholds. PubChem molecular fingerprints comprising 881 binary features were calculated and reduced to 178 features via variance threshold filtering. Exploratory analysis confirmed statistically significant differences between active and inactive compounds across all Lipinski descriptors. Thirty-five regression algorithms were benchmarked, followed by five-fold cross-validation and hyperparameter optimization. The optimized Random Forest and Histogram Gradient Boosting regressors achieved the best predictive performance with coefficient of determination values of 0.78 and root mean squared error of 0.58 on the held-out test set. These results demonstrate practical utility for virtual screening of novel drug candidates and provide a reproducible computational framework for accelerating early-stage Alzheimer's drug discovery.

**Keywords:** Alzheimer's Disease, Quantitative Structure-Activity Relationship, Machine Learning, Beta-Amyloid A4, Drug Discovery, Random Forest

---

## I. Introduction

Alzheimer's disease is a chronic neurodegenerative disorder and the leading cause of dementia, affecting approximately 55 million people globally [1]. The disease is characterized by the accumulation of extracellular amyloid-beta plaques and intracellular neurofibrillary tangles, leading to synaptic failure and neuronal death [2]. The amyloid hypothesis posits that beta-amyloid peptide aggregation, derived from proteolysis of the Amyloid Precursor Protein, is the primary driver of pathogenesis [3]. Consequently, beta-amyloid A4 protein remains a critical therapeutic target.

Despite decades of research, attrition rates for Alzheimer's drug candidates remain exceptionally high. Traditional drug discovery pipelines require over a decade and approximately 2.2 billion dollars to bring a new molecular entity to market [4]. Computational methods, specifically Quantitative Structure-Activity Relationship (QSAR) modeling, offer a cost-effective alternative by predicting biological activity from chemical structure [5].

This study leverages machine learning to construct robust QSAR models for beta-amyloid A4 inhibitors. By correlating molecular fingerprints with bioactivity data from the ChEMBL database, we identify structural features governing potency and provide a predictive framework for virtual screening.

---

## II. Materials and Methods

### A. Data Acquisition and Curation

Experimental bioactivity data was sourced from the ChEMBL database, a manually curated repository of bioactive molecules [6]. The target "Beta-amyloid A4 protein" (ChEMBL ID: CHEMBL2487) was queried, yielding 7,918 initial entries.

Filtration criteria included:
1. Bioactivity type restricted to IC50 measurements
2. Units standardized to nanomolar
3. Removal of entries with missing values
4. Duplicate SMILES handling

The final curated dataset comprised 1,319 unique compounds.

### B. Data Preprocessing and Labeling

Compounds were labeled based on IC50 values:
- **Active:** IC50 ≤ 1,000 nM
- **Intermediate:** 1,000 nM < IC50 < 10,000 nM
- **Inactive:** IC50 ≥ 10,000 nM

IC50 values were converted to pIC50 using:

```
pIC50 = -log₁₀(IC50 × 10⁻⁹)
```

This transformation normalizes the distribution for regression modeling.

### C. Feature Engineering

Two descriptor types were calculated:

1. **Lipinski Descriptors:** Molecular weight, LogP, hydrogen bond donors, and hydrogen bond acceptors were computed using RDKit [7] for drug-likeness assessment per Lipinski's Rule of Five [8].

2. **Molecular Fingerprints:** SMILES structures were converted to 881-bit PubChem fingerprints using PaDEL-Descriptor [9]. Variance threshold filtering (threshold = 0.16) reduced features to 178.

### D. Machine Learning Pipeline

Data was split into training (80%, n=1,055) and testing (20%, n=264) sets. LazyPredict [10] screened 35 regression algorithms. Top performers underwent optimization:

- **Random Forest Regressor:** Ensemble bagging method [11]
- **HistGradientBoosting Regressor:** Histogram-based gradient boosting
- **Support Vector Regression:** Evaluated for overfitting robustness

Hyperparameter tuning used GridSearchCV with 5-fold cross-validation [12].

---

## III. Results

### A. Exploratory Data Analysis

Mann-Whitney U tests confirmed statistically significant differences (p < 0.05) between active and inactive compounds across all Lipinski descriptors:
- pIC50: p < 0.001
- Molecular Weight: p < 0.001
- LogP: p < 0.05
- Hydrogen Bond Donors: p < 0.001
- Hydrogen Bond Acceptors: p < 0.001

### B. Model Performance

Tree-based ensemble methods outperformed linear and distance-based algorithms:

| Model | Train R² | Test R² | RMSE |
|-------|----------|---------|------|
| Random Forest | 0.95 | 0.78 | 0.58 |
| HistGradientBoosting | 0.94 | 0.78 | 0.59 |
| XGBRegressor | 0.97 | 0.77 | 0.61 |
| SVR | 0.86 | 0.72 | 0.67 |

### C. Overfitting Analysis

Random Forest and HistGradientBoosting achieved highest test accuracy (R² = 0.78) with moderate overfitting indicated by training-testing gaps. SVR demonstrated lower overfitting (Train R² = 0.86, Test R² = 0.72) but reduced predictive capacity.

---

## IV. Discussion

The computational pipeline successfully identifies beta-amyloid A4 inhibitors with R² = 0.78, indicating that 178 PubChem fingerprint features capture sufficient structural information for bioactivity prediction.

Statistical significance of Lipinski descriptors confirms that bioactivity depends on physicochemical properties including lipophilicity and molecular size, which influence blood-brain barrier permeability [13].

Limitations include moderate overfitting in tree-based models due to high dimensionality relative to sample size. Future work could employ recursive feature elimination or Graph Neural Networks for improved generalization [14].

---

## V. Conclusion

This study developed robust QSAR models predicting Alzheimer's drug candidate potency with high accuracy. Random Forest and HistGradientBoosting emerged as optimal algorithms. The computational framework enables rapid virtual screening of chemical libraries, significantly reducing time and cost in early-stage Alzheimer's drug discovery.

---

## References

[1] Alzheimer's Disease International. (2024). *World Alzheimer report 2024: Global changes in attitudes to dementia*. https://www.alzint.org/

[2] Alzheimer's Association. (2025). 2025 Alzheimer's disease facts and figures. *Alzheimer's & Dementia*, *21*(4), 1598–1695. https://doi.org/10.1002/alz.13824

[3] Selkoe, D. J., & Hardy, J. (2016). The amyloid hypothesis of Alzheimer's disease at 25 years. *EMBO Molecular Medicine*, *8*(6), 595–608. https://doi.org/10.15252/emmm.201606210

[4] Deloitte. (2024). *Measuring the return from pharmaceutical innovation 2024*. Deloitte Centre for Health Solutions.

[5] Cherkasov, A., Muratov, E. N., Fourches, D., Varnek, A., Baskin, I. I., Cronin, M., ... & Tropsha, A. (2014). QSAR modeling: Where have you been? Where are you going to? *Journal of Medicinal Chemistry*, *57*(12), 4977–5010. https://doi.org/10.1021/jm4004285

[6] Zdrazil, B., Felix, E., Hunter, F., Manber, E. J., Nowotka, M., Papadatos, G., ... & Leach, A. R. (2024). The ChEMBL database in 2023: A drug discovery platform spanning multiple bioactivity data types and time periods. *Nucleic Acids Research*, *52*(D1), D1180–D1192. https://doi.org/10.1093/nar/gkad1004

[7] Landrum, G. (2024). *RDKit: Open-source cheminformatics software*. https://www.rdkit.org/

[8] Lipinski, C. A., Lombardo, F., Dominy, B. W., Feeney, P. J., & Adler, A. D. (1997). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Advanced Drug Delivery Reviews*, *23*(1–3), 3–25. https://doi.org/10.1016/S0169-409X(96)00423-1

[9] Yap, C. W. (2011). PaDEL-descriptor: An open source software to calculate molecular descriptors and fingerprints. *Journal of Computational Chemistry*, *32*(7), 1466–1474. https://doi.org/10.1002/jcc.21707

[10] Pandala, S. R. (2022). *LazyPredict* [Computer software]. https://pypi.org/project/lazypredict/

[11] Breiman, L. (2001). Random forests. *Machine Learning*, *45*(1), 5–32. https://doi.org/10.1023/A:1010933404324

[12] Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., ... & Duchesnay, É. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, *12*, 2825–2830.

[13] Pardridge, W. M. (2012). Drug transport across the blood–brain barrier. *Journal of Cerebral Blood Flow & Metabolism*, *32*(11), 1959–1972. https://doi.org/10.1038/jcbfm.2012.126

[14] Cha, Y., Erickson, R. I., Bhatt, D. L., & Bhatt, A. P. (2024). Navigating the frontiers of machine learning in neurodegenerative disease therapeutics. *Pharmaceuticals*, *17*(2), 158. https://doi.org/10.3390/ph17020158
