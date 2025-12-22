# Index: Machine Learning for Alzheimer's Drug Discovery

This index maps the entire project workflow—from biological theory to the final Python code—explaining the **justification** (Why), the **action** (What), and the **implementation** (How) for every step.

---

## Quick Reference

| Resource | Description |
|----------|-------------|
| [glossary.md](glossary.md) | Domain terminology definitions |
| [libraries.md](libraries.md) | Complete library reference with import statements |
| [linux.md](linux.md) | Shell commands used in project |
| [notebook-analysis.md](notebook-analysis.md) | Detailed cell-by-cell notebook evaluation |

---

## Part I: The Foundation (Setting the Stage)

### **Chapter 1: The Digital Laboratory Setup**
* **1.1. Environment Configuration**
    * **Topic:** Virtual Environment & Dependencies.
    * **Action (What):** Installing `chembl_webresource_client`, `rdkit`, `lazypredict` and mounting Google Drive.
    * **Justification (Why):** Standard data science environments do not include specialized chemistry tools. We need a persistent workspace (Drive) to prevent data loss and specialized libraries to talk to biological databases.
    * **Implementation (How):** `!pip install`, `drive.mount('/content/drive')`, `os.chdir()`.
    * **Notebook Cells:** 1-5

### **Chapter 2: The Biological Target (Domain Knowledge)**
* **2.1. Understanding the Enemy: Alzheimer's Disease**
    * **Topic:** The Amyloid Cascade Hypothesis.
    * **Action (What):** Identifying **Beta-amyloid A4 protein (APP)** as the specific drug target.
    * **Justification (Why):** You cannot build a model without a clear objective. This protein is responsible for creating the plaques that kill neurons in Alzheimer's patients; stopping it is the goal.
    * **Notebook Cells:** 6-10 (Section 1: "Alzheimer's disease: an ongoing challenge")
* **2.2. The Strategy: QSAR**
    * **Topic:** Quantitative Structure-Activity Relationship.
    * **Action (What):** Adopting the hypothesis that **Biological Activity = Function(Chemical Structure)**.
    * **Justification (Why):** This mathematical framework allows us to predict the efficacy of *untested* chemicals using computer models, saving billions in lab costs.

---

## Part II: Data Mining & Curation (The Target Variable $y$)

### **Chapter 3: Mining the ChEMBL Database**
* **3.1. Database Querying**
    * **Topic:** Fetching Bioactivity Data.
    * **Action (What):** Querying ChEMBL for ID `CHEMBL2487` (Beta-amyloid A4).
    * **Justification (Why):** Machine Learning requires "ground truth" data to learn from. ChEMBL is the gold-standard repository for experimental drug results.
    * **Implementation (How):** `new_client.target.filter(target_chembl_id='CHEMBL2487')`.
    * **Notebook Cells:** 11-20
    * **Output:** Raw bioactivity data (~7,918 rows from API or Excel import)
* **3.2. Data Standardization**
    * **Topic:** Filtering for **IC50**.
    * **Action (What):** Retaining only experiments that measured "Half Maximal Inhibitory Concentration."
    * **Justification (Why):** Experimental results come in many units (Ki, Kd, EC50). Mixing them would confuse the model. We must enforce a single, consistent unit of potency (IC50).
    * **Implementation (How):** `bioactivity[bioactivity.standard_type == 'IC50']`.
    * **Output:** `Beta_amyloid_A4_protein_bioactivity_data_curated.csv` (1,319 rows)

### **Chapter 4: Data Preprocessing & Labeling**
* **4.1. Cleaning & Handling Missing Data**
    * **Topic:** Data Sanitization.
    * **Action (What):** Removing rows with missing `standard_value` or `canonical_smiles`.
    * **Justification (Why):** Models cannot calculate properties for non-existent structures or learn from empty target values.
    * **Implementation (How):** `df.dropna()`.
    * **Notebook Cells:** 21-35
* **4.2. Bioactivity Classification**
    * **Topic:** Thresholding.
    * **Action (What):** Categorizing drugs into **Active** ($\le$ 1,000 nM), **Inactive** ($\ge$ 10,000 nM), and **Intermediate**.
    * **Justification (Why):** Raw numbers are noisy. Creating clear "buckets" allows for simple "Good vs. Bad" analysis and enables classification modeling.
    * **Implementation (How):** Custom iteration loop classifying based on `float(standard_value)`.
    * **Output:** `Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv`

---

## Part III: Feature Engineering (The Input Features $X$)

### **Chapter 5: Chemical Descriptors & Lipinski's Rules**
* **5.1. The "Rule of 5" Analysis**
    * **Topic:** Drug-Likeness.
    * **Action (What):** Calculating Molecular Weight, LogP, and H-Bond donors/acceptors.
    * **Justification (Why):** Potency isn't enough; a drug must be absorbable. This filter identifies if "Active" compounds are physically capable of being oral pills (not too heavy, not too greasy).
    * **Implementation (How):** `Descriptors.MolWt(mol)`, `Descriptors.MolLogP(mol)` using **RDKit**.
    * **Notebook Cells:** 36-70
    * **Libraries:** `rdkit.Chem`, `rdkit.Chem.Descriptors`, `rdkit.Chem.Lipinski`
* **5.2. Data Normalization (The pIC50 Transformation)**
    * **Topic:** Logarithmic Scaling.
    * **Action (What):** Converting raw **IC50** (nM) to **pIC50** ($-log_{10}(M)$).
    * **Justification (Why):**
        1.  **Distribution:** IC50 data is skewed (spanning 1 to 100,000). Log transformation fixes this for regression.
        2.  **Intuition:** It converts the confusing "Lower is Better" scale into a standard "Higher is Better" scale.
    * **Implementation (How):** Custom function `norm_value()` and `pIC50()`.
    * **Output:** `Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_3class_pIC50.csv`

### **Chapter 6: Molecular Fingerprinting (PaDEL)**
* **6.1. From Text to Math**
    * **Topic:** PubChem Fingerprints.
    * **Action (What):** Converting SMILES strings (text) into binary vectors of **881 bits** (0s and 1s).
    * **Justification (Why):** Algorithms (like Random Forest) cannot read text descriptions of molecules. They need numerical arrays representing chemical substructures (e.g., "Has Ring? Yes/No").
    * **Implementation (How):** `!bash padel.sh` (Using the Java-based PaDEL engine).
    * **Notebook Cells:** 71-85
    * **Input:** `molecule.smi` (SMILES strings)
    * **Output:** `descriptors_output.csv` (881 PubChem fingerprint columns)

---

## Part IV: Exploratory Data Analysis (EDA)

### **Chapter 7: Statistical Validation**
* **7.1. Hypothesis Testing**
    * **Topic:** Mann-Whitney U Test.
    * **Action (What):** Comparing the distributions of descriptors between Active and Inactive groups.
    * **Justification (Why):** To prove that observed differences (e.g., "Active drugs are smaller") are **statistically significant** and not just due to random chance.
    * **Implementation (How):** `mannwhitneyu()` from `scipy.stats`.
    * **Notebook Cells:** 86-110
    * **Findings:** All four Lipinski descriptors show statistically significant differences (p < 0.05):
        - MW: p = 4.44e-05
        - LogP: p = 3.51e-12
        - NumHDonors: p = 0.0005
        - NumHAcceptors: p = 6.74e-06

### **Chapter 8: Visualization (The Chemical Space)**
* **8.1. Chemical Space Plots**
    * **Topic:** Scatter & Frequency Plots.
    * **Action (What):** Plotting MW vs. LogP and pIC50 frequency.
    * **Justification (Why):** Visual confirmation that "Active" compounds cluster in specific, drug-like regions of chemical space (validating Lipinski's rules).
    * **Implementation (How):** `matplotlib` and `seaborn`.
    * **Notebook Cells:** 111-130
    * **Output Files:**
        - `plot_MW.pdf` — Molecular Weight distribution
        - `plot_LogP.pdf` — LogP distribution
        - `plot_NumHDonors.pdf` — H-bond donors distribution
        - `plot_NumHAcceptors.pdf` — H-bond acceptors distribution
        - `plot_MW_vs_LogP.pdf` — Chemical space scatter plot
        - `plot_bioactivity_class.pdf` — Class distribution
        - `plot_ic50.pdf` — IC50 distribution

---

## Part V: Machine Learning Model Building

### **Chapter 9: Model Preparation**
* **9.1. Feature Selection**
    * **Topic:** Removing Low Variance Features.
    * **Action (What):** Deleting fingerprint columns that are identical for almost all molecules.
    * **Justification (Why):** If a feature is "0" for 99% of the data, it provides no information. Removing it speeds up training and reduces the "Curse of Dimensionality."
    * **Implementation (How):** `VarianceThreshold(threshold=0.1)`.
    * **Notebook Cells:** 131-145
    * **Result:** 881 features → 178 features after variance filtering
* **9.2. Train/Test Split**
    * **Topic:** Data Partitioning.
    * **Action (What):** Splitting data into 80% Training and 20% Testing.
    * **Justification (Why):** We must evaluate the model on data it has *never seen* to ensure it isn't just memorizing the training set (Overfitting).
    * **Implementation (How):** `train_test_split(X, Y, test_size=0.2)`.

### **Chapter 10: The Algorithm Tournament**
* **10.1. LazyPredict Benchmarking**
    * **Topic:** Model Comparison.
    * **Action (What):** Training 30+ Regression Models simultaneously.
    * **Justification (Why):** There is no single "best" AI model. Instead of guessing, we run a tournament to empirically find which algorithm works best for this specific chemical dataset.
    * **Implementation (How):** `LazyRegressor.fit(X_train, X_test, Y_train, Y_test)`.
    * **Notebook Cells:** 146-160
    * **Libraries:** `lazypredict.Supervised.LazyRegressor`

---

## Part VI: Evaluation & Conclusion

### **Chapter 11: The Report Card**
* **11.1. Performance Metrics**
    * **Topic:** $R^2$ and RMSE.
    * **Action (What):** Calculating the Coefficient of Determination ($R^2$) and Root Mean Squared Error.
    * **Justification (Why):** These metrics quantify reliability. A high $R^2$ (>0.6) proves the model has learned the underlying chemical rules.
    * **Findings:** Random Forest and Decision Trees typically perform best (~0.71).
    * **Notebook Cells:** 161-166
* **11.2. Visualization of Accuracy**
    * **Topic:** Experimental vs. Predicted Plot.
    * **Action (What):** Plotting Real Lab Values ($X$) vs. AI Predictions ($Y$).
    * **Justification (Why):** A visual sanity check. Ideally, all points fall on a straight diagonal line, indicating perfect prediction.

### **Chapter 12: Final Product**
* **12.1. The QSAR Model**
    * **Outcome:** A validated computational tool.
    * **Justification (Why):** We can now screen millions of *new, unmade* chemical structures and predict their potential to treat Alzheimer's Disease, drastically accelerating the drug discovery pipeline.

---

## Data Pipeline Summary

```
Raw Data (7,918 rows)
    ↓ Filter IC50, remove nulls
Curated Data (1,319 rows)
    ↓ Add bioactivity labels
Labelled Data (1,319 rows)
    ↓ Remove intermediate class
3-Class Data (Active + Inactive)
    ↓ Calculate pIC50
pIC50 Data
    ↓ Generate PubChem fingerprints
Fingerprint Data (881 features)
    ↓ Variance threshold filtering
ML-Ready Data (178 features)
    ↓ Train/Test split (80/20)
Model Training → Predictions
```

## Output Files Generated

| File | Generated In | Description |
|------|--------------|-------------|
| `Beta_amyloid_A4_protein_bioactivity_data_curated.csv` | Section 2 | Filtered for IC50, nM units |
| `Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv` | Section 2 | With bioactivity class labels |
| `Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_2class_pIC50.csv` | Section 3 | Binary classification |
| `Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_3class_pIC50.csv` | Section 3 | 3-class with pIC50 |
| `molecule.smi` | Section 5 | SMILES for fingerprinting |
| `descriptors_output.csv` | Section 5 | PubChem fingerprints |
| `plot_*.pdf` | Section 4 | EDA visualizations |
