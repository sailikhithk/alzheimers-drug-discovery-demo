# Glossary: Machine Learning for Alzheimer's Drug Discovery

This glossary defines the technical terms, biological concepts, and data science metrics used throughout the Alzheimer's Drug Discovery QSAR project.

---

## Phase 1: Biological Foundation & Data Mining

### **Alzheimer's Disease**
A progressive neurodegenerative disease causing 60-70% of dementia cases. Characterized by memory loss, cognitive decline, and behavioral changes. The disease involves abnormal accumulation of beta-amyloid plaques in the brain.

### **ChEMBL**
A massive, open-access database maintained by the European Bioinformatics Institute (EBI). It contains information on bioactive molecules with drug-like properties.
* *Context:* Used as the primary data source to fetch experimental results for the target protein.

### **Target Protein (Beta-amyloid A4)**
The specific biological molecule (ID: `CHEMBL2487`) implicated in Alzheimer's Disease.
* *Mechanism:* This protein cleaves into sticky fragments (Aβ) that aggregate into plaques, leading to neuron death.
* *Goal:* To find chemical compounds that bind to this target and inhibit its activity.

### **Inhibitor**
A chemical compound (drug) that binds to a target protein and decreases its biological activity.
* *Context:* We are specifically filtering for inhibitors that stop the Beta-amyloid aggregation process.

### **IC50 (Half Maximal Inhibitory Concentration)**
A quantitative measure of a drug's potency. It indicates how much of a particular drug is needed to inhibit a biological process by 50%.
* *Scale:* **Lower is Better**. A low IC50 (e.g., 5 nM) means only a tiny amount is needed to work. A high IC50 (e.g., 50,000 nM) means a huge amount is needed (weak drug).

### **Standard Units (nM)**
The unit of measurement for IC50, typically **Nanomolar** ($10^{-9}$ mol/L).
* *Context:* We standardize all data to nM to ensure consistent comparison between different experiments.

---

## Phase 2: Data Cleaning & Labeling

### **SMILES**
**Simplified Molecular Input Line Entry System**. A specification for describing the structure of chemical molecules using short ASCII strings.
* *Example:* `CCO` represents Ethanol.
* *Function:* Allows 3D chemical structures to be stored as text for computer processing.

### **Canonical SMILES**
A "standardized" unique text representation of a molecule. Since a SMILES string can be written in multiple ways for the same molecule (depending on where you start counting atoms), the "Canonical" version ensures every unique molecule has exactly one unique identifier.

### **Bioactivity Class**
A categorical label created to bucket compounds based on their potency.
* **Active:** IC50 ≤ 1,000 nM (Strong candidates) — corresponds to pIC50 ≥ 6
* **Inactive:** IC50 ≥ 10,000 nM (Weak/Useless candidates) — corresponds to pIC50 ≤ 5
* **Intermediate:** 1,000 < IC50 < 10,000 nM (Borderline candidates)
* *Why these thresholds:* Industry standard cutoffs. 1,000 nM (1 μM) is a common "hit" threshold; 10,000 nM (10 μM) is a typical screening cutoff.

### **Missing Data (NA)**
Rows in the dataset where critical values (like `standard_value` or `canonical_smiles`) are undefined. These are dropped during preprocessing.

---

## Phase 3: Feature Engineering (Chemistry to Math)

### **Lipinski's Rule of 5**
A rule of thumb used to evaluate **drug-likeness**. It predicts if a chemical compound with a certain pharmacological or biological activity has properties that would make it a likely orally active drug in humans.
* **Molecular Weight (MW):** Should be < 500 Daltons (Size).
* **LogP:** Should be < 5 (Solubility).
* **H-Bond Donors:** Should be ≤ 5.
* **H-Bond Acceptors:** Should be ≤ 10.
* *Why "Rule of 5":* All thresholds are multiples of 5 — easy to remember.

### **LogP (Partition Coefficient)**
A measure of the differential solubility of a compound in two immiscible solvents (octanol and water).
* **High LogP (>5):** Lipophilic (Fat-loving). Can cross cell membranes but may accumulate in fat tissue.
* **Low LogP (<0):** Hydrophilic (Water-loving). Dissolves easily in blood but struggles to cross membranes.
* **Ideal range:** 1-3 for oral drugs.

### **pIC50**
The negative logarithm of the IC50 value in Molar units.
* **Formula:** pIC50 = -log₁₀(IC50 in Molar)
* **Conversion from nM:** pIC50 = -log₁₀(IC50_nM × 10⁻⁹) = 9 - log₁₀(IC50_nM)
* **Why used?**
    1. **Linearity:** Converts the exponential scale of IC50 into a linear scale suitable for regression.
    2. **Intuition:** Converts "Lower is Better" to "**Higher is Better**."
    3. **Distribution:** Transforms skewed IC50 data to more normal distribution.

---

## Phase 4: Molecular Descriptors (The "Fingerprints")

### **Molecular Descriptor**
A mathematical representation of a molecule's physical or chemical properties (e.g., molecular weight, number of rings, charge).

### **Molecular Fingerprint**
A binary vector (array of 0s and 1s) representing the presence or absence of specific substructures in a molecule.
* *Analogy:* Like a barcode for a molecule.

### **PubChem Fingerprint**
A specific fingerprinting system that checks for 881 specific chemical substructures.
* *Output:* A list of 881 binary bits for every molecule.
* *Categories:*
  - Bits 1-115: Element counts (Has ≥1 Carbon, Has ≥2 Oxygen, etc.)
  - Bits 116-263: Ring systems (Has benzene, Has 5-member ring, etc.)
  - Bits 264-459: Atom pairs (C-C bond, C-N bond, C=O bond, etc.)
  - Bits 460-579: Atom environments (Carbon with 3 neighbors, etc.)
  - Bits 580-881: SMARTS patterns (specific substructure patterns)

### **PaDEL-Descriptor**
A software tool used to calculate molecular descriptors and fingerprints from SMILES strings.
* *Requirement:* Java runtime environment.

---

## Phase 5: Exploratory Data Analysis (EDA)

### **Chemical Space**
A multidimensional space defined by the molecular descriptors.
* *Visualization:* We plot MW vs. LogP to visualize if "Active" compounds cluster in a specific "drug-like" region of this space.

### **Mann-Whitney U Test**
A non-parametric statistical test used to determine whether there is a significant difference between the distributions of two independent groups (Active vs. Inactive).
* *Why non-parametric:* Doesn't assume normal distribution of data.
* *Goal:* To prove that the difference in properties (like MW) between active and inactive drugs is not just due to random chance.

### **P-Value**
The probability that the observed results occurred by chance.
* **p < 0.05:** Statistically significant (The difference is real).
* **p ≥ 0.05:** Not statistically significant (Cannot reject null hypothesis).

### **Project Results (Mann-Whitney U Test)**
All five descriptors showed statistically significant differences between active and inactive compounds:
* pIC50: p = 1.70e-155 ✓
* Molecular Weight: p = 7.80e-18 ✓
* LogP: p = 4.24e-02 ✓
* NumHDonors: p = 1.45e-11 ✓
* NumHAcceptors: p = 3.49e-14 ✓

---

## Phase 6: Model Building (Machine Learning)

### **Feature Selection (Variance Threshold)**
A technique to reduce the number of input variables ($X$).
* *Process:* Removing columns (fingerprints) that have the same value (0 or 1) for almost all samples, as they provide no useful information to the model.
* *Threshold used:* 0.8 × (1 - 0.8) = 0.16
* *Result:* Reduced from 881 features to 177 features.

### **Train/Test Split**
Dividing the dataset into two subsets:
* **Training Set (80%):** Used to train the model (~1,055 samples).
* **Testing Set (20%):** Used to evaluate the model's performance on unseen data (~264 samples).
* *Purpose:* To detect and prevent **Overfitting**.

### **Overfitting**
When a model learns the training data too well, including noise and outliers, and fails to generalize to new data.
* *Symptom:* High training accuracy, low test accuracy.
* *Project observation:* Training R² = 0.97, Test R² = 0.37 — classic overfitting.

### **LazyPredict**
A Python library that helps build a baseline by running many machine learning models with default parameters on a dataset and ranking them by performance metrics.
* *Models tested:* 30+ regression algorithms.

### **Random Forest Regressor**
An ensemble learning method that constructs a multitude of decision trees at training time and outputs the mean prediction of the individual trees.
* *Project result:* Test R² = 0.29

### **HistGradientBoostingRegressor**
A gradient boosting algorithm optimized for speed using histogram-based splitting.
* *Project result:* Best test performance with R² = 0.37

---

## Phase 7: Evaluation (Grading the AI)

### **R² (Coefficient of Determination)**
A statistical measure that represents the proportion of the variance for a dependent variable ($y$) that's explained by an independent variable ($X$).
* *Range:* 0 to 1 (can be negative for poor models).
* *Interpretation:* R² = 0.37 means model explains 37% of pIC50 variation.
* *Goal for QSAR:* > 0.6 is generally considered good.

### **RMSE (Root Mean Squared Error)**
The standard deviation of the prediction errors (residuals). It tells you how concentrated the data is around the line of best fit.
* *Goal:* A lower value indicates better fit.

### **Scatter Plot (Experimental vs. Predicted)**
A visualization comparing the ground truth values from the lab (X-axis) against the values predicted by the AI (Y-axis).
* *Ideal:* All points falling on a straight diagonal line.

### **Project Results Summary**
* Best test R²: 0.37 (HistGradientBoostingRegressor)
* Significant overfitting observed (Training R² ~0.97 vs Test R² ~0.37)
* All Lipinski descriptors statistically significant for distinguishing active/inactive compounds
* **Next steps:** Feature engineering, hyperparameter tuning, cross-validation, ensemble methods
