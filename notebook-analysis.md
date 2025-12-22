# Notebook Analysis: ML for Alzheimer's Drug Discovery

This document provides a sequential evaluation of the Jupyter notebook `ml_for_alzheimer_s_drug_discovery.ipynb`, analyzing whether markdown descriptions and code cells are consistent and make sense together.

---

## Overview

| Metric | Value |
|--------|-------|
| Total Sections | 7 main sections |
| Data Source | ChEMBL Database (Beta amyloid A4 protein) |
| Initial Compounds | 7,918 |
| Final Dataset | 1,319 compounds (after filtering) |
| Features | 881 PubChem fingerprints → 178 after variance filtering |
| Target Variable | pIC50 (continuous) |

---

## Section 1: Alzheimer's Disease Background

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 1 | MD | Title: "Alzheimers Drug Discovery Project Demo with QSAR" | ✅ OK | Clear project title |
| 2 | MD | Separator (***) | ✅ OK | Visual separator |
| 3 | Code | Google Colab drive mount | ✅ OK | Environment setup |
| 4 | Code | Set working directory, import numpy/pandas | ✅ OK | Standard setup |
| 5 | MD | Alzheimer's disease explanation | ✅ OK | Comprehensive background on AD, Beta-amyloid A4 protein, therapeutic approaches |

**Section Assessment**: ✅ GOOD - Provides excellent domain context for data scientists unfamiliar with Alzheimer's research. Explains why Beta-amyloid A4 is a drug target.

---

## Section 2: QSAR Introduction

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 6 | MD | Section header "# 2. QSAR..." | ✅ OK | Clear section start |
| 7 | MD | QSAR explanation with Figure 1 reference | ✅ OK | Explains ML workflow for drug discovery |
| 8 | MD | Detailed QSAR workflow explanation | ✅ OK | Describes molecular descriptors, x/y variables, model training |
| 9 | MD | Project objective (yellow box) | ✅ OK | Clear goal statement |

**Section Assessment**: ✅ GOOD - Excellent pedagogical flow. Explains QSAR concept before diving into implementation.

---

## Section 3: Download Bioactivity Data

### 3.1 Installing and Importing Libraries

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 10 | MD | Section header "# 3. Download Bioactivity Data" | ✅ OK | |
| 11 | MD | ChEMBL database introduction | ✅ OK | Explains data source |
| 12 | MD | "## 3.1 Installing and importing libraries" | ✅ OK | |
| 13 | Code | pip install commands (commented out) | ✅ OK | Dependencies listed |
| 14 | Code | `from chembl_webresource_client.new_client import new_client` | ✅ OK | Import statement |

### 3.2 Search for Target Proteins

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 15 | MD | "## 3.2 Search for target proteins..." | ✅ OK | |
| 16 | MD | Explains ChEMBL data retrieval | ✅ OK | Notes API is slow |
| 17 | Code | Target search (commented out) | ✅ OK | Shows API approach but uses pre-downloaded data |
| 18 | MD | Explains using xlsx file instead | ✅ OK | Practical workaround |
| 19 | Code | `pd.read_excel('Beta_amyloid A4_protein_active_compounds.xlsx')` | ✅ OK | Loads 7,918 compounds |
| 20 | MD | "7,918 entries with 45 features" | ✅ OK | Matches output |

### 3.3 Data Cleaning and Labeling

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 21 | MD | "## 3.3 Data Cleaning and Labeling" | ✅ OK | |
| 22 | MD | Focus on specific features | ✅ OK | |
| 23 | Code | Check unique "Target Organism" | ✅ OK | Output: ['Homo sapiens'] |
| 24 | Code | Check unique "Target Type" | ✅ OK | Output: ['SINGLE PROTEIN'] |
| 25 | MD | Filter for IC50 only | ✅ OK | Explains uniformity need |
| 26 | Code | Filter `Standard Type == 'IC50'` | ✅ OK | Reduces to 1,497 rows |
| 27 | MD | Explains IC50 filtering | ✅ OK | |
| 28 | Code | Filter for nM units | ✅ OK | |
| 29 | Code | Drop null Standard Value | ✅ OK | |
| 30 | Code | Drop duplicate SMILES | ✅ OK | Final: 1,319 compounds |
| 31 | Code | Save curated data | ✅ OK | |

### 3.4 Labeling Compounds

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 32 | MD | "## 3.4 Labeling compounds..." | ✅ OK | |
| 33 | MD | Explains IC50 thresholds | ✅ OK | Active ≤1000nM, Inactive ≥10000nM |
| 34 | Code | Define `bioactivity_class()` function | ✅ OK | Implements 3-class labeling |
| 35 | Code | Apply labeling | ✅ OK | |
| 36 | Code | Value counts | ✅ OK | Shows class distribution |

**Section Assessment**: ✅ GOOD - Clear data pipeline. Each filtering step is explained before code execution.

---

## Section 4: Exploratory Data Analysis (EDA)

### 4.1 Calculating Lipinski Descriptors

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 37 | MD | "# 4. Exploratory Data Analysis" | ✅ OK | |
| 38 | MD | "## 4.1 Calculating Lipinski Descriptors" | ✅ OK | |
| 39 | MD | Lipinski Rule of 5 explanation | ✅ OK | Explains drug-likeness criteria |
| 40 | Code | Import RDKit | ✅ OK | |
| 41 | Code | Define `lipinski()` function | ✅ OK | Calculates MW, LogP, HBD, HBA |
| 42 | Code | Apply Lipinski descriptors | ✅ OK | |

### 4.2 Convert IC50 to pIC50

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 43 | MD | "## 4.2 Convert IC50 to pIC50" | ✅ OK | |
| 44 | MD | Explains pIC50 transformation | ✅ OK | pIC50 = -log10(IC50 in M) |
| 45 | Code | Define `pIC50()` function | ✅ OK | Handles unit conversion |
| 46 | Code | Apply pIC50 calculation | ✅ OK | |

### 4.3 Removing Intermediate Class

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 47 | MD | "## 4.3 Removing the 'intermediate' bioactivity class" | ✅ OK | |
| 48 | MD | Explains binary classification rationale | ✅ OK | |
| 49 | Code | Filter to active/inactive only | ✅ OK | Creates df_2class |

### 4.4 EDA via Lipinski Descriptors

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 50 | MD | "## 4.4 Exploratory Data Analysis via Lipinski Descriptors" | ✅ OK | |
| 51 | MD | "### 4.4.1 How many active and inactive compounds?" | ✅ OK | |
| 52 | Code | Bioactivity class bar plot | ✅ OK | Visualization |
| 53 | MD | "### 4.4.2 Molecular Weight and LogP chemical spaces" | ✅ OK | |
| 54 | Code | Scatter plot MW vs LogP | ✅ OK | |
| 55 | MD | "### 4.4.3 Distribution of pIC50 Values" | ✅ OK | |
| 56 | Code | pIC50 box plot | ✅ OK | |
| 57 | MD | "### 4.4.4 Distribution of Molecular Weights" | ✅ OK | |
| 58 | Code | MW box plot | ✅ OK | |
| 59 | MD | "### 4.4.5 Distribution of LogP Values" | ✅ OK | |
| 60 | Code | LogP box plot | ✅ OK | |
| 61 | MD | "### 4.4.6 Distribution of Hydrogen Bond Donors" | ✅ OK | |
| 62 | Code | NumHDonors box plot | ✅ OK | |
| 63 | MD | "### 4.4.7 Distribution of Hydrogen Bond Acceptors" | ✅ OK | |
| 64 | Code | NumHAcceptors box plot | ✅ OK | |

### 4.5 Statistical Validation: Mann-Whitney U Test

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 65 | MD | "## 4.5 Statistical Validation: Mann-Whitney U Test" | ✅ OK | |
| 66 | MD | Explains why Mann-Whitney U (non-parametric, two groups) | ✅ OK | Good pedagogical explanation |
| 67 | Code | Define `mann_whitney_test()` function | ✅ OK | |
| 68 | Code | Run test on 5 descriptors | ✅ OK | All p < 0.05 |
| 69 | MD | Interpretation of results | ✅ OK | Explains significance |

### 4.6 Key Findings: EDA

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 70 | MD | "## 4.6 Key Findings: Exploratory Data Analysis" | ✅ OK | |
| 71 | MD | pIC50 findings (green box) | ✅ OK | Correctly states significant |
| 72 | MD | Lipinski findings (green box) | ✅ OK | **FIXED** - Now correctly states all 4 descriptors significant |

**Section Assessment**: ✅ GOOD - Comprehensive EDA with statistical validation. Section 4.5 and 4.6 are now consistent (both confirm significance).

---

## Section 5: Descriptor Calculation and Dataset Preparation

### 5.1 Calculate Fingerprint Descriptors

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 73 | MD | "# 5. Descriptor Calculation and Dataset Preparation" | ✅ OK | |
| 74 | MD | "## 5.1 Calculate Fingerprint Descriptors" | ✅ OK | |
| 75 | MD | Select SMILES and ChEMBL ID | ✅ OK | |
| 76 | Code | Create molecule.smi file | ✅ OK | Tab-separated SMILES |
| 77 | Code | Verify molecule.smi content | ✅ OK | |
| 78 | Code | Count molecules (1,319) | ✅ OK | |

### 5.1.1 Calculate PaDEL Descriptors

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 79 | MD | "### 5.1.1 Calculate PaDEL Descriptors" | ✅ OK | |
| 80 | Code | PaDEL calculation (commented out) | ⚠️ NOTE | Uses pre-computed file |
| 81 | Code | List files to verify descriptors_output.csv | ✅ OK | |

### 5.2 Preparing X and Y Matrices

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 82 | MD | "## 5.2 Preparing the X and Y Data Matrices" | ✅ OK | |
| 83 | MD | "### 5.2.1 Prepare X matrix" | ✅ OK | |
| 84 | MD | Explains X = PaDEL features | ✅ OK | |
| 85 | Code | Load descriptors, create df3_X | ✅ OK | 881 features |
| 86 | MD | "### 5.2.2 Prepare Y matrix" | ✅ OK | |
| 87 | Code | Create df3_Y (pIC50 values) | ✅ OK | |
| 88 | MD | "### 5.2.3 Combining X and Y" | ✅ OK | |
| 89 | Code | Merge X and Y | ✅ OK | |

**Section Assessment**: ✅ GOOD - Clear explanation of feature engineering. PaDEL fingerprints (881 binary features) are well documented.

---

## Section 6: Model Building

### 6.1 Remove Low Variance Features

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 90 | MD | "## 6.1 Remove Low Variance Features" | ✅ OK | |
| 91 | MD | Explains variance threshold rationale | ✅ OK | |
| 92 | Code | VarianceThreshold(0.8 * 0.2) | ✅ OK | 881 → 178 features |
| 93 | MD | "177 features with variance ≥80%" | ⚠️ MINOR | Output shows 178, text says 177 |

### 6.2 Data Split

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 94 | MD | "## 6.2 Data split" | ✅ OK | |
| 95 | Code | Assign Y = df3_Y | ✅ OK | |
| 96 | MD | Explains 80/20 split | ✅ OK | |
| 97 | Code | train_test_split(test_size=0.2) | ✅ OK | |
| 98 | Code | X_train.shape: (1055, 178) | ✅ OK | |
| 99 | Code | X_test.shape: (264, 178) | ✅ OK | |
| 100 | MD | "1055 training, 264 testing" | ✅ OK | Matches output |
| 101 | Code | Print sample counts | ✅ OK | |

### 6.3 Compare ML Models

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 102 | MD | "## 6.3 Compare ML Models" | ✅ OK | |
| 103 | MD | Explains LazyRegressor approach | ✅ OK | |
| 104 | Code | Import LazyRegressor | ✅ OK | |
| 105 | MD | LazyRegressor explanation | ✅ OK | |
| 106 | Code | Run LazyRegressor | ✅ OK | Trains 30+ models |
| 107 | Code | Display results table | ✅ OK | Shows R², RMSE, Time |

**Section Assessment**: ✅ GOOD - Clear ML pipeline. Minor discrepancy in feature count (177 vs 178).

---

## Section 7: Provisional Conclusions

| Cell | Type | Content | Status | Notes |
|------|------|---------|--------|-------|
| 108 | MD | "# 7. Provisional Conclusions" | ✅ OK | |
| 109 | MD | Comprehensive conclusions | ✅ OK | **UPDATED** - Now accurately reflects results |

**Conclusions Content**:
- ✅ Statistical validation: All 5 descriptors significant (p < 0.05)
- ✅ Training R²: 0.94-0.97 (ExtraTrees, RandomForest, DecisionTree)
- ✅ Test R²: Max 0.76 (RandomForestRegressor & HistGradientBoostingRegressor)
- ✅ Identifies overfitting problem (tree models show moderate to severe overfitting)
- ✅ Lists next steps for improvement

**Section Assessment**: ✅ GOOD - Conclusions now accurately reflect the actual results.

---

## Summary: Issues Found and Fixed

| Issue | Location | Status |
|-------|----------|--------|
| Section 4.6 said "no significant difference" for Lipinski descriptors | Cell 72 | ✅ FIXED - Now says "all significant" |
| Conclusions didn't reflect actual model performance | Cell 109 | ✅ FIXED - Now shows correct R² values |
| Feature count discrepancy (177 vs 178) | Section 6.1 | ⚠️ MINOR - Text says 177, output shows 178 |

---

## Overall Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Pedagogical Flow** | ⭐⭐⭐⭐⭐ | Excellent progression from biology → chemistry → ML |
| **Code-Markdown Consistency** | ⭐⭐⭐⭐ | Good after fixes; minor feature count discrepancy |
| **Domain Explanations** | ⭐⭐⭐⭐⭐ | Thorough explanations of AD, QSAR, Lipinski, IC50 |
| **Statistical Rigor** | ⭐⭐⭐⭐⭐ | Mann-Whitney U test validates visual observations |
| **ML Pipeline** | ⭐⭐⭐⭐ | Good use of LazyPredict; identifies overfitting |
| **Conclusions** | ⭐⭐⭐⭐⭐ | Honest about limitations, suggests improvements |

**Overall**: The notebook is well-structured for teaching a data scientist about computational drug discovery. The flow from disease background → data collection → EDA → feature engineering → ML modeling is logical and well-documented.


---

## Appendix A: Detailed CSV File Generation Mapping

This section traces exactly where each CSV file is generated in the notebook.

### File 1: Beta_amyloid_A4_protein_bioactivity_data_curated.csv

| Attribute | Value |
|-----------|-------|
| Generated at | Section 3.3, Cell 82 |
| DataFrame | `df4` |
| Code | `df4.to_csv('Beta_amyloid_A4_protein_bioactivity_data_curated.csv', index=False)` |
| Rows | 1,319 |
| Columns | 4 |

**Columns:**
```
Molecule ChEMBL ID, Molecular Weight, Smiles, Standard Value
```

**Preceding Markdown:**
> "Let's saves the final dataframe to CSV file"

**How df4 was created:**
1. Load Excel → df (7,918 rows)
2. Filter `Standard Type == 'IC50'` → 1,497 rows
3. Filter `Standard Units == 'nM'` → 1,351 rows
4. Drop null `Standard Value` → 1,319 rows
5. Select 4 columns → df4

---

### File 2: Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv

| Attribute | Value |
|-----------|-------|
| Generated at | Section 3.4, Cell 85 |
| DataFrame | `df5` |
| Code | `df5.to_csv('Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv', index=False)` |
| Rows | 1,319 |
| Columns | 5 |

**Columns:**
```
Molecule ChEMBL ID, Molecular Weight, Smiles, Standard Value, class
```

**Preceding Markdown:**
> "## 3.4 Labeling compounds as either being active, inactive or intermediate"

**How df5 was created:**
```python
# Bioactivity classification thresholds
def bioactivity_class(value):
    if value <= 1000:
        return 'active'
    elif value >= 10000:
        return 'inactive'
    else:
        return 'intermediate'

df5 = pd.concat([df4, bioactivity_class], axis=1)
```

**Class Distribution:**
- active: 449
- intermediate: 375
- inactive: 495

---

### File 3: Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_3class_pIC50.csv

| Attribute | Value |
|-----------|-------|
| Generated at | Section 4.2, Cell 94 |
| DataFrame | `df_combined` |
| Code | `df_combined.to_csv('Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_3class_pIC50.csv')` |
| Rows | 1,319 |
| Columns | 9 |

**Columns:**
```
(index), Molecule ChEMBL ID, Molecular Weight, Smiles, pIC50, class, LogP, NumHDonors, NumHAcceptors
```

**Preceding Markdown:**
> "Let's write this to CSV file"

**How df_combined was created:**
1. Calculate Lipinski descriptors using RDKit
2. Calculate pIC50: `pIC50 = -log10(IC50 * 1e-9)`
3. Combine all into df_combined

**Note:** This file has an index column (unnamed first column) because `index=False` was not specified.

---

### File 4: Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_2class_pIC50.csv

| Attribute | Value |
|-----------|-------|
| Generated at | Section 4.3, Cell 96 |
| DataFrame | `df_2class` |
| Code | `df_2class.to_csv('Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_2class_pIC50.csv')` |
| Rows | 944 |
| Columns | 9 |

**Columns:**
```
(index), Molecule ChEMBL ID, Molecular Weight, Smiles, pIC50, class, LogP, NumHDonors, NumHAcceptors
```

**Preceding Markdown:**
> "## 4.3 Removing the 'intermediate' bioactivity class"
> "In order to easily compare between active and inactive compounds, we have to delete from our analysis intermediate compounds"

**How df_2class was created:**
```python
df_2class = df_combined[df_combined['class'] != 'intermediate']
```

**Class Distribution:**
- active: 449
- inactive: 495

**Purpose:** Used for EDA visualizations and Mann-Whitney U test (comparing active vs inactive only).

---

### File 5: molecule.smi

| Attribute | Value |
|-----------|-------|
| Generated at | Section 5.1, Cell 116 |
| DataFrame | `df3_selection` |
| Code | `df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)` |
| Rows | 1,319 |
| Format | Tab-separated, no header |

**Content:**
```
SMILES<tab>ChEMBL_ID
```

**Example:**
```
Nc1nc(NCc2ccccc2)c2ccccc2n1	CHEMBL3815078
CC(C)Nc1nc(NCCc2ccccc2)c2ccccc2n1	CHEMBL4070462
```

**Purpose:** Input file for PaDEL-Descriptor to calculate PubChem fingerprints.

---

### File 6: Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv

| Attribute | Value |
|-----------|-------|
| Generated at | Section 5.2.3, Cell 125 |
| DataFrame | `dataset3` |
| Code | `dataset3.to_csv('Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv', index=False)` |
| Rows | 1,319 |
| Columns | 883 |

**Columns:**
```
Molecular Weight, PubchemFP0, PubchemFP1, ..., PubchemFP880, pIC50
```

**Preceding Markdown:**
> "dataset3 dataframe contains 881 input features and 1 output variable (pIC50 values)"
> "Let's download the CSV file for the Part 4 (Model Building)"

**How dataset3 was created:**
```python
dataset3 = pd.concat([df3_X, df3_Y], axis=1)
# df3_X = 881 PubChem fingerprints + Molecular Weight
# df3_Y = pIC50 values
```

**Purpose:** Final ML-ready dataset with features (X) and target (Y).

---

### File 7: descriptors_output.csv (INCOMPLETE)

| Attribute | Value |
|-----------|-------|
| Generated by | PaDEL-Descriptor (external tool) |
| Rows | 418 (should be 1,319) |
| Columns | 882 |

**Columns:**
```
Name, PubchemFP0, PubchemFP1, ..., PubchemFP880
```

**Note:** This file is incomplete (only 418 of 1,319 molecules). The notebook uses pre-computed fingerprints from `Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv` instead.

---

### File 8: aligned_dataset.csv

| Attribute | Value |
|-----------|-------|
| Generated by | `create_aligned_data.py` (utility script) |
| Rows | 1,319 |
| Columns | 888 |

**Columns:**
```
Molecule ChEMBL ID, Molecular Weight, Smiles, Standard Value, class, Name, PubchemFP0, ..., PubchemFP880, pIC50
```

**Purpose:** Alternative merged dataset created by utility script. Contains bioactivity data + fingerprints + pIC50.

---

## Appendix B: Data Files Analysis

### Input File: Beta_amyloid A4_protein_active_compounds.xlsx

| Attribute | Value |
|-----------|-------|
| Rows | 7,918 |
| Columns | 45 |
| Source | ChEMBL Database |
| Target | Beta amyloid A4 protein (CHEMBL2487) |

**Key Columns:**
| Column | Description | Example |
|--------|-------------|---------|
| Molecule ChEMBL ID | Unique compound identifier | CHEMBL3815078 |
| Smiles | Molecular structure (SMILES notation) | Nc1nc(NCc2ccccc2)c2ccccc2n1 |
| Standard Type | Measurement type | IC50, Inhibition, Activity, Ki |
| Standard Value | Measured value | 4800.0 |
| Standard Units | Units of measurement | nM, %, uM |
| Molecular Weight | Compound mass (Da) | 250.31 |
| AlogP | Lipophilicity | 2.82 |
| Target Name | Protein target | Beta amyloid A4 protein |
| Assay Description | Experimental details | Full assay protocol |

**Standard Type Distribution:**
| Type | Count | Notes |
|------|-------|-------|
| Inhibition | 4,179 | % inhibition at fixed concentration |
| IC50 | 1,497 | Half-maximal inhibitory concentration ✓ Used |
| Activity | 1,160 | General activity measure |
| Ki | 531 | Inhibition constant |
| Other | 551 | Potency, EC50, Kd, etc. |

**Why IC50 was chosen:**
- Quantitative measure (not just % at single concentration)
- Standard drug discovery metric
- Allows pIC50 transformation for ML

---

### Data Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────┐
│ STAGE 1: Raw Input                                                  │
│ File: Beta_amyloid A4_protein_active_compounds.xlsx                 │
│ Rows: 7,918 | Cols: 45                                              │
│ Contains: All assay types (IC50, Inhibition, Ki, Activity, etc.)    │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼ Filter: Standard Type == 'IC50'
                                    │         Standard Units == 'nM'
                                    │         Standard Value not null
┌─────────────────────────────────────────────────────────────────────┐
│ STAGE 2: Curated                                                    │
│ File: Beta_amyloid_A4_protein_bioactivity_data_curated.csv          │
│ Rows: 1,319 | Cols: 4                                               │
│ Columns: Molecule ChEMBL ID, Molecular Weight, Smiles, Standard Value│
│ Note: Contains 196 duplicate SMILES (not removed at this stage)     │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼ Add bioactivity class labels
                                    │ Active: IC50 ≤ 1,000 nM
                                    │ Intermediate: 1,000 < IC50 < 10,000 nM
                                    │ Inactive: IC50 ≥ 10,000 nM
┌─────────────────────────────────────────────────────────────────────┐
│ STAGE 3: Labelled (3-class)                                         │
│ File: Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv│
│ Rows: 1,319 | Cols: 5                                               │
│ Class distribution: active=449, intermediate=375, inactive=495      │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                    ┌───────────────┴───────────────┐
                    ▼                               ▼
┌─────────────────────────────────┐ ┌─────────────────────────────────┐
│ STAGE 4a: 3-class + pIC50       │ │ STAGE 4b: 2-class (Binary)      │
│ + Lipinski descriptors          │ │ File: ..._2class_pIC50.csv      │
│ File: ..._3class_pIC50.csv      │ │ Rows: 944 | Cols: 9             │
│ Rows: 1,319 | Cols: 9           │ │ (Intermediate removed)          │
│ Added: pIC50, LogP, NumHDonors, │ │ Class: active=449, inactive=495 │
│        NumHAcceptors            │ │ Used for: EDA, Mann-Whitney test│
└─────────────────────────────────┘ └─────────────────────────────────┘
                    │
                    ▼ Calculate PubChem fingerprints (881 bits)
┌─────────────────────────────────────────────────────────────────────┐
│ STAGE 5: With Fingerprints                                          │
│ File: Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv│
│ Rows: 1,319 | Cols: 883                                             │
│ Features: Molecular Weight + 881 PubChem fingerprints (PubchemFP0-880)│
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼ Variance threshold (keep variance ≥ 0.16)
                                    │ 881 → 178 features
┌─────────────────────────────────────────────────────────────────────┐
│ STAGE 6: ML Ready                                                   │
│ X: 1,319 samples × 178 features                                     │
│ Y: pIC50 values (continuous target)                                 │
│ Split: 80% train (1,055) / 20% test (264)                           │
└─────────────────────────────────────────────────────────────────────┘
```

---

### Intermediate Files Summary

| File | Rows | Cols | Purpose |
|------|------|------|---------|
| raw_input.csv | 7,918 | 45 | Converted from Excel for analysis |
| Beta_amyloid_A4_protein_bioactivity_data_curated.csv | 1,319 | 4 | IC50 only, nM units, no nulls |
| Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv | 1,319 | 5 | Added 3-class labels |
| Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_2class_pIC50.csv | 944 | 9 | Binary (no intermediate), for EDA |
| Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled_3class_pIC50.csv | 1,319 | 9 | 3-class with pIC50 + Lipinski |
| Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv | 1,319 | 883 | Final ML dataset with fingerprints |
| molecule.smi | 1,319 | 2 | SMILES for PaDEL input |
| descriptors_output.csv | 418 | 882 | Partial PaDEL output (incomplete) |
| aligned_dataset.csv | 1,319 | 888 | Alternative merged dataset |

---

### Key Observations

1. **Duplicate SMILES**: 196 duplicate SMILES exist in the curated data (1,319 rows but only 1,123 unique SMILES). These are NOT removed, meaning some molecules appear multiple times with different IC50 values from different assays.

2. **descriptors_output.csv Mismatch**: This file only has 418 rows (vs 1,319 expected). The notebook uses the pre-computed `Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv` instead.

3. **IC50 Range**: 0.3 nM to 800,000 nM (spanning 6 orders of magnitude), justifying the log transformation to pIC50.

4. **Class Balance**: Reasonably balanced - active (449), intermediate (375), inactive (495).

5. **pIC50 Range**: 3.10 to 9.52 (corresponding to IC50 of ~800,000 nM to ~0.3 nM).
