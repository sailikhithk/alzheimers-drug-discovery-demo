# Alzheimer's Drug Discovery using Machine Learning

A computational drug discovery project that uses QSAR (Quantitative Structure-Activity Relationship) modeling to predict the bioactivity of compounds against Beta-amyloid A4 protein - a key target in Alzheimer's disease.

---

## Table of Contents
1. [Project Overview](#1-project-overview)
2. [The Science Behind It](#2-the-science-behind-it)
3. [Data Pipeline - Detailed](#3-data-pipeline---detailed)
4. [Code Walkthrough](#4-code-walkthrough)
5. [File Descriptions](#5-file-descriptions)
6. [Key Concepts Explained](#6-key-concepts-explained)
7. [How to Run](#7-how-to-run)

---

## 1. Project Overview

**Goal:** Build machine learning models that can predict whether a chemical compound will be effective (active) against Beta-amyloid A4 protein.

**Why it matters:** If we can predict which molecules are likely to be active, we can:
- Save years of lab testing
- Reduce drug development costs
- Accelerate finding treatments for Alzheimer's

**What the ML model does:**
- Input: A molecule's structure (represented as numbers/features)
- Output: Predicted bioactivity (pIC50 value - how potent the drug is)

---

## 2. The Science Behind It

### 2.1 Alzheimer's Disease & Beta-Amyloid
- Alzheimer's causes 60-70% of dementia cases
- Beta-amyloid A4 protein forms "plaques" in the brain
- These plaques disrupt neurons and cause brain cell death
- **Our target:** Find molecules that can inhibit/block this protein

### 2.2 What is IC50?
- **IC50** = "Half maximal Inhibitory Concentration"
- It measures how much drug is needed to inhibit the target by 50%
- **Lower IC50 = More potent drug** (needs less to work)
- Measured in nanomolar (nM)

### 2.3 Bioactivity Classification
The project classifies compounds into 3 categories based on IC50:

| Class | IC50 Value | Meaning |
|-------|------------|---------|
| **Active** | < 1,000 nM | Strong inhibitor (good drug candidate) |
| **Intermediate** | 1,000 - 10,000 nM | Moderate effect |
| **Inactive** | > 10,000 nM | Weak/no effect |

### 2.4 What is QSAR?
**QSAR = Quantitative Structure-Activity Relationship**

The core idea: A molecule's **structure** determines its **activity**.

```
Molecule Structure → Mathematical Features → ML Model → Predicted Activity
     (SMILES)         (Fingerprints)                      (pIC50)
```

---

## 3. Data Pipeline - Detailed

### High-Level Overview

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                        COMPLETE DATA PIPELINE OVERVIEW                                 ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║   [ChEMBL Database] ──→ [Raw Excel] ──→ [Cleaned CSV] ──→ [Labelled CSV]              ║
║                                                              │                         ║
║                                                              ▼                         ║
║   [ML Model] ◄── [Train/Test Split] ◄── [Feature Matrix] ◄── [Fingerprints]          ║
║       │                                                                                ║
║       ▼                                                                                ║
║   [Predictions] ──→ [Drug Candidates]                                                  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

---

### STEP 1: Raw Data Acquisition

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 1: RAW DATA FROM ChEMBL DATABASE                                                 ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  SOURCE: ChEMBL (https://www.ebi.ac.uk/chembl/)                                       ║
║  - World's largest open bioactivity database                                           ║
║  - Contains millions of compound-target interactions                                   ║
║  - Curated from scientific literature                                                  ║
║                                                                                        ║
║  FILE: Beta_amyloid A4_protein_active_compounds.xlsx                                   ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ COLUMNS IN RAW DATA (~30+ columns):                                              │  ║
║  │                                                                                  │  ║
║  │ • Molecule ChEMBL ID    - Unique identifier (e.g., CHEMBL3815078)               │  ║
║  │ • Smiles                - Molecular structure as text                            │  ║
║  │ • Standard Type         - Type of measurement (IC50, Ki, EC50, etc.)            │  ║
║  │ • Standard Value        - The measured value                                     │  ║
║  │ • Standard Units        - Units (nM, uM, etc.)                                  │  ║
║  │ • Molecular Weight      - Mass in Daltons                                        │  ║
║  │ • Target Organism       - Species (Homo sapiens)                                 │  ║
║  │ • Target Type           - SINGLE PROTEIN                                         │  ║
║  │ • Assay Description     - How the experiment was done                            │  ║
║  │ • Document Year         - Publication year                                       │  ║
║  │ • ... and many more                                                              │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  RECORDS: ~1,320 compounds tested against Beta-amyloid A4 protein                     ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 2: Data Cleaning & Filtering

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 2: DATA CLEANING & STANDARDIZATION                                               ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  OUTPUT FILE: Beta_amyloid_A4_protein_bioactivity_data_curated.csv                    ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ FILTERING OPERATIONS:                                                            │  ║
║  │                                                                                  │  ║
║  │ 1. FILTER BY MEASUREMENT TYPE                                                    │  ║
║  │    ┌──────────────────────────────────────────────────────────────────────────┐ │  ║
║  │    │ Raw data has multiple measurement types:                                  │ │  ║
║  │    │   • IC50 (Half maximal Inhibitory Concentration) ← WE KEEP THIS          │ │  ║
║  │    │   • Ki (Inhibition constant)                                              │ │  ║
║  │    │   • EC50 (Half maximal Effective Concentration)                           │ │  ║
║  │    │   • Kd (Dissociation constant)                                            │ │  ║
║  │    │                                                                           │ │  ║
║  │    │ WHY IC50? Most common, standardized, comparable across studies            │ │  ║
║  │    │                                                                           │ │  ║
║  │    │ CODE: df_ic50 = df[df['Standard Type'] == 'IC50']                        │ │  ║
║  │    └──────────────────────────────────────────────────────────────────────────┘ │  ║
║  │                                                                                  │  ║
║  │ 2. FILTER BY UNITS                                                               │  ║
║  │    ┌──────────────────────────────────────────────────────────────────────────┐ │  ║
║  │    │ Different labs report in different units:                                 │ │  ║
║  │    │   • nM (nanomolar) ← WE KEEP THIS (most common)                          │ │  ║
║  │    │   • uM (micromolar) = 1000 nM                                             │ │  ║
║  │    │   • mM (millimolar) = 1,000,000 nM                                        │ │  ║
║  │    │                                                                           │ │  ║
║  │    │ WHY nM? Standard unit, avoids conversion errors                           │ │  ║
║  │    │                                                                           │ │  ║
║  │    │ CODE: df2 = df_ic50[df_ic50['Standard Units'] == 'nM']                   │ │  ║
║  │    └──────────────────────────────────────────────────────────────────────────┘ │  ║
║  │                                                                                  │  ║
║  │ 3. REMOVE NULL VALUES                                                            │  ║
║  │    ┌──────────────────────────────────────────────────────────────────────────┐ │  ║
║  │    │ Some experiments failed or weren't recorded properly                      │ │  ║
║  │    │                                                                           │ │  ║
║  │    │ CODE: df3 = df2[~df2['Standard Value'].isnull()]                         │ │  ║
║  │    └──────────────────────────────────────────────────────────────────────────┘ │  ║
║  │                                                                                  │  ║
║  │ 4. SELECT RELEVANT COLUMNS                                                       │  ║
║  │    ┌──────────────────────────────────────────────────────────────────────────┐ │  ║
║  │    │ From 30+ columns, we keep only 4:                                         │ │  ║
║  │    │                                                                           │ │  ║
║  │    │ • Molecule ChEMBL ID  → Unique identifier                                 │ │  ║
║  │    │ • Molecular Weight    → For Lipinski analysis                             │ │  ║
║  │    │ • Smiles              → Molecular structure (for fingerprints)            │ │  ║
║  │    │ • Standard Value      → IC50 value (our target variable)                  │ │  ║
║  │    └──────────────────────────────────────────────────────────────────────────┘ │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  RESULT: Clean dataset with ~1,320 rows × 4 columns                                   ║
║                                                                                        ║
║  SAMPLE OUTPUT:                                                                        ║
║  ┌──────────────────┬──────────────────┬─────────────────────────┬──────────────┐     ║
║  │ Molecule ChEMBL  │ Molecular Weight │ Smiles                  │ Standard Val │     ║
║  ├──────────────────┼──────────────────┼─────────────────────────┼──────────────┤     ║
║  │ CHEMBL3815078    │ 250.31           │ Nc1nc(NCc2ccccc2)c2...  │ 4800.0       │     ║
║  │ CHEMBL4070462    │ 306.41           │ CC(C)Nc1nc(NCCc2ccc...  │ 3800.0       │     ║
║  │ CHEMBL3911126    │ 361.40           │ COc1cc(O)ccc1/C=C/C...  │ 520.0        │     ║
║  │ CHEMBL3943188    │ 434.49           │ CO[C@H]1CC[C@]2(CC1...  │ 4.8          │     ║
║  └──────────────────┴──────────────────┴─────────────────────────┴──────────────┘     ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 3: Bioactivity Labeling

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 3: BIOACTIVITY CLASSIFICATION (LABELING)                                         ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  OUTPUT FILE: Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv       ║
║                                                                                        ║
║  PURPOSE: Convert continuous IC50 values into categorical labels                       ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ CLASSIFICATION LOGIC:                                                            │  ║
║  │                                                                                  │  ║
║  │                        IC50 Scale (nM)                                           │  ║
║  │   ◄────────────────────────────────────────────────────────────────────────►    │  ║
║  │   0.1    1     10    100   1,000  10,000  100,000  1,000,000                     │  ║
║  │    │     │      │      │      │       │        │         │                       │  ║
║  │    └─────┴──────┴──────┴──────┘       └────────┴─────────┘                       │  ║
║  │              │                              │                                     │  ║
║  │         ACTIVE                         INACTIVE                                   │  ║
║  │       (< 1,000 nM)                   (> 10,000 nM)                                │  ║
║  │                                                                                  │  ║
║  │                    ┌─────────────────┐                                           │  ║
║  │                    │  INTERMEDIATE   │                                           │  ║
║  │                    │ (1,000-10,000)  │                                           │  ║
║  │                    └─────────────────┘                                           │  ║
║  │                                                                                  │  ║
║  │ THRESHOLDS:                                                                      │  ║
║  │   • IC50 < 1,000 nM      → "active"       (potent inhibitor)                    │  ║
║  │   • 1,000 ≤ IC50 ≤ 10,000 → "intermediate" (moderate effect)                    │  ║
║  │   • IC50 > 10,000 nM     → "inactive"     (weak/no effect)                      │  ║
║  │                                                                                  │  ║
║  │ WHY THESE THRESHOLDS?                                                            │  ║
║  │   • Industry standard for drug discovery                                         │  ║
║  │   • 1,000 nM = 1 μM (common cutoff for "hit" compounds)                         │  ║
║  │   • 10,000 nM = 10 μM (typical screening cutoff)                                │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  CODE:                                                                                 ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ for ic50 in df["Standard Value"]:                                                │  ║
║  │     if ic50 >= 10000:                                                            │  ║
║  │         label = "inactive"                                                       │  ║
║  │     elif ic50 <= 1000:                                                           │  ║
║  │         label = "active"                                                         │  ║
║  │     else:                                                                        │  ║
║  │         label = "intermediate"                                                   │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  RESULT: Added 'class' column                                                          ║
║                                                                                        ║
║  CLASS DISTRIBUTION (approximate):                                                     ║
║  ┌────────────────┬─────────┬────────────────────────────────────────────────────┐    ║
║  │ Class          │ Count   │ Bar                                                │    ║
║  ├────────────────┼─────────┼────────────────────────────────────────────────────┤    ║
║  │ active         │ ~450    │ ████████████████████                               │    ║
║  │ intermediate   │ ~350    │ ███████████████                                    │    ║
║  │ inactive       │ ~520    │ ██████████████████████                             │    ║
║  └────────────────┴─────────┴────────────────────────────────────────────────────┘    ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 4: Exploratory Data Analysis (EDA)

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 4: EXPLORATORY DATA ANALYSIS (EDA)                                               ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  PURPOSE: Understand the data before building ML models                                ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ 4A. CALCULATE LIPINSKI DESCRIPTORS                                               │  ║
║  │                                                                                  │  ║
║  │ Using RDKit library to compute drug-likeness properties:                         │  ║
║  │                                                                                  │  ║
║  │ SMILES String ──► RDKit Molecule Object ──► Descriptors                          │  ║
║  │                                                                                  │  ║
║  │ ┌────────────────────────────────────────────────────────────────────────────┐  │  ║
║  │ │ DESCRIPTOR        │ WHAT IT MEASURES              │ DRUG-LIKE RANGE        │  │  ║
║  │ ├────────────────────────────────────────────────────────────────────────────┤  │  ║
║  │ │ LogP              │ Fat vs water solubility       │ -0.4 to 5.6            │  │  ║
║  │ │ NumHDonors        │ Hydrogen bond donors          │ ≤ 5                    │  │  ║
║  │ │ NumHAcceptors     │ Hydrogen bond acceptors       │ ≤ 10                   │  │  ║
║  │ │ Molecular Weight  │ Size of molecule              │ ≤ 500 Da               │  │  ║
║  │ └────────────────────────────────────────────────────────────────────────────┘  │  ║
║  │                                                                                  │  ║
║  │ CODE:                                                                            │  ║
║  │   from rdkit import Chem                                                         │  ║
║  │   from rdkit.Chem import Descriptors, Lipinski                                   │  ║
║  │   mol = Chem.MolFromSmiles(smiles_string)                                        │  ║
║  │   logp = Descriptors.MolLogP(mol)                                                │  ║
║  │   hdonors = Lipinski.NumHDonors(mol)                                             │  ║
║  │   hacceptors = Lipinski.NumHAcceptors(mol)                                       │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ 4B. CONVERT IC50 TO pIC50                                                        │  ║
║  │                                                                                  │  ║
║  │ PROBLEM: IC50 values span a huge range (0.3 to 300,000 nM)                       │  ║
║  │ SOLUTION: Use logarithmic scale (pIC50)                                          │  ║
║  │                                                                                  │  ║
║  │ FORMULA:                                                                         │  ║
║  │ ┌────────────────────────────────────────────────────────────────────────────┐  │  ║
║  │ │                                                                            │  │  ║
║  │ │   pIC50 = -log₁₀(IC50 in Molar)                                           │  │  ║
║  │ │                                                                            │  │  ║
║  │ │   Step 1: Convert nM to M  →  IC50_M = IC50_nM / 1,000,000,000            │  │  ║
║  │ │   Step 2: Take negative log →  pIC50 = -log₁₀(IC50_M)                     │  │  ║
║  │ │                                                                            │  │  ║
║  │ └────────────────────────────────────────────────────────────────────────────┘  │  ║
║  │                                                                                  │  ║
║  │ CONVERSION TABLE:                                                                │  ║
║  │ ┌────────────────┬────────────────┬────────────────────────────────────────┐   │  ║
║  │ │ IC50 (nM)      │ pIC50          │ Interpretation                         │   │  ║
║  │ ├────────────────┼────────────────┼────────────────────────────────────────┤   │  ║
║  │ │ 1              │ 9.0            │ Extremely potent (rare)                │   │  ║
║  │ │ 10             │ 8.0            │ Very potent                            │   │  ║
║  │ │ 100            │ 7.0            │ Potent                                 │   │  ║
║  │ │ 1,000          │ 6.0            │ Active threshold ◄──────────────────── │   │  ║
║  │ │ 10,000         │ 5.0            │ Inactive threshold ◄────────────────── │   │  ║
║  │ │ 100,000        │ 4.0            │ Weak                                   │   │  ║
║  │ └────────────────┴────────────────┴────────────────────────────────────────┘   │  ║
║  │                                                                                  │  ║
║  │ BENEFIT: Higher pIC50 = More potent (intuitive!)                                │  ║
║  │          Better distribution for ML algorithms                                   │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ 4C. VISUALIZATIONS GENERATED                                                     │  ║
║  │                                                                                  │  ║
║  │ • plot_bioactivity_class.pdf  → Bar chart of active vs inactive counts          │  ║
║  │ • plot_MW_vs_LogP.pdf         → Scatter plot showing chemical space             │  ║
║  │ • plot_ic50.pdf               → Box plot of pIC50 by class                      │  ║
║  │ • plot_MW.pdf                 → Box plot of Molecular Weight by class           │  ║
║  │ • plot_LogP.pdf               → Box plot of LogP by class                       │  ║
║  │ • plot_NumHAcceptors.pdf      → Box plot of H-bond acceptors by class           │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 5: Molecular Fingerprint Calculation

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 5: MOLECULAR FINGERPRINT CALCULATION                                             ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  PURPOSE: Convert molecular structures into numerical features for ML                  ║
║                                                                                        ║
║  INPUT FILE:  molecule.smi (SMILES strings)                                           ║
║  OUTPUT FILE: descriptors_output.csv (881 fingerprint features)                       ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ WHAT ARE MOLECULAR FINGERPRINTS?                                                 │  ║
║  │                                                                                  │  ║
║  │ A way to represent molecular structure as a binary vector (0s and 1s)            │  ║
║  │                                                                                  │  ║
║  │ EXAMPLE:                                                                         │  ║
║  │                                                                                  │  ║
║  │   Molecule: Aspirin                                                              │  ║
║  │   SMILES: CC(=O)OC1=CC=CC=C1C(=O)O                                              │  ║
║  │                                                                                  │  ║
║  │   Fingerprint (simplified):                                                      │  ║
║  │   ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┐                │  ║
║  │   │ Bit │ Bit │ Bit │ Bit │ Bit │ Bit │ Bit │ Bit │ ... │ Bit │                │  ║
║  │   │  1  │  2  │  3  │  4  │  5  │  6  │  7  │  8  │     │ 881 │                │  ║
║  │   ├─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤                │  ║
║  │   │  1  │  0  │  1  │  1  │  0  │  0  │  1  │  0  │ ... │  1  │                │  ║
║  │   └─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┘                │  ║
║  │                                                                                  │  ║
║  │   Each bit represents presence (1) or absence (0) of a substructure:            │  ║
║  │   • Bit 1: Has benzene ring? → 1 (yes)                                          │  ║
║  │   • Bit 2: Has nitrogen? → 0 (no)                                               │  ║
║  │   • Bit 3: Has carboxylic acid? → 1 (yes)                                       │  ║
║  │   • ... and so on for 881 different substructures                               │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ TOOL USED: PaDEL-Descriptor                                                      │  ║
║  │                                                                                  │  ║
║  │ • Calculates PubChem fingerprints (881 bits)                                     │  ║
║  │ • PubChem = NIH's public chemistry database                                      │  ║
║  │ • Standardized fingerprint definition                                            │  ║
║  │                                                                                  │  ║
║  │ PROCESS:                                                                         │  ║
║  │                                                                                  │  ║
║  │   ┌─────────────┐      ┌─────────────┐      ┌─────────────────────────┐         │  ║
║  │   │ molecule.smi│ ───► │   PaDEL     │ ───► │ descriptors_output.csv  │         │  ║
║  │   │ (SMILES)    │      │ (Java tool) │      │ (881 columns)           │         │  ║
║  │   └─────────────┘      └─────────────┘      └─────────────────────────┘         │  ║
║  │                                                                                  │  ║
║  │ CODE:                                                                            │  ║
║  │   # Save SMILES to file                                                          │  ║
║  │   df[['Smiles','Molecule ChEMBL ID']].to_csv('molecule.smi', sep='\t')          │  ║
║  │                                                                                  │  ║
║  │   # Run PaDEL (bash command)                                                     │  ║
║  │   ! bash padel.sh                                                                │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ PUBCHEM FINGERPRINT CATEGORIES (881 bits total):                                 │  ║
║  │                                                                                  │  ║
║  │ ┌──────────────────────────────────────────────────────────────────────────┐    │  ║
║  │ │ CATEGORY                    │ BITS      │ EXAMPLES                       │    │  ║
║  │ ├──────────────────────────────────────────────────────────────────────────┤    │  ║
║  │ │ Element counts              │ 1-115     │ Has ≥1 Carbon, Has ≥2 Oxygen   │    │  ║
║  │ │ Ring systems                │ 116-263   │ Has benzene, Has 5-member ring │    │  ║
║  │ │ Atom pairs                  │ 264-459   │ C-C bond, C-N bond, C=O bond   │    │  ║
║  │ │ Atom nearest neighbors      │ 460-579   │ C with 3 neighbors             │    │  ║
║  │ │ SMARTS patterns             │ 580-881   │ Specific substructure patterns │    │  ║
║  │ └──────────────────────────────────────────────────────────────────────────┘    │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  OUTPUT SAMPLE (descriptors_output.csv):                                              ║
║  ┌──────────────────┬────────┬────────┬────────┬─────┬────────┐                       ║
║  │ Name             │ PubFP0 │ PubFP1 │ PubFP2 │ ... │ PubFP880│                      ║
║  ├──────────────────┼────────┼────────┼────────┼─────┼────────┤                       ║
║  │ CHEMBL3815078    │ 1      │ 1      │ 1      │ ... │ 0      │                       ║
║  │ CHEMBL4070462    │ 1      │ 1      │ 1      │ ... │ 0      │                       ║
║  │ CHEMBL3911126    │ 1      │ 1      │ 1      │ ... │ 1      │                       ║
║  └──────────────────┴────────┴────────┴────────┴─────┴────────┘                       ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 6: Prepare ML Dataset (X and Y Matrices)

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 6: PREPARE ML DATASET (X AND Y MATRICES)                                         ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  OUTPUT FILE: Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv    ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ COMBINING FEATURES (X) AND TARGET (Y)                                            │  ║
║  │                                                                                  │  ║
║  │ ┌─────────────────────────────────────────────────────────────────────────────┐ │  ║
║  │ │                                                                             │ │  ║
║  │ │   X (Features)                              Y (Target)                      │ │  ║
║  │ │   ┌────────────────────────────────────┐   ┌─────────┐                     │ │  ║
║  │ │   │ PubFP0 │ PubFP1 │ ... │ PubFP880   │   │ pIC50   │                     │ │  ║
║  │ │   ├────────┼────────┼─────┼────────────┤   ├─────────┤                     │ │  ║
║  │ │   │   1    │   1    │ ... │     0      │   │  5.32   │                     │ │  ║
║  │ │   │   1    │   1    │ ... │     0      │   │  5.42   │                     │ │  ║
║  │ │   │   1    │   1    │ ... │     1      │   │  6.28   │                     │ │  ║
║  │ │   │   1    │   0    │ ... │     0      │   │  8.32   │                     │ │  ║
║  │ │   │  ...   │  ...   │ ... │    ...     │   │  ...    │                     │ │  ║
║  │ │   └────────┴────────┴─────┴────────────┘   └─────────┘                     │ │  ║
║  │ │        881 columns                          1 column                        │ │  ║
║  │ │                                                                             │ │  ║
║  │ └─────────────────────────────────────────────────────────────────────────────┘ │  ║
║  │                                                                                  │  ║
║  │ CODE:                                                                            │  ║
║  │   # Load fingerprints                                                            │  ║
║  │   df3_X = pd.read_csv('descriptors_output.csv')                                  │  ║
║  │   df3_X = df3_X.drop(columns=['Name'])  # Remove ID column                       │  ║
║  │                                                                                  │  ║
║  │   # Get target variable                                                          │  ║
║  │   df3_Y = df_padel['pIC50']                                                      │  ║
║  │                                                                                  │  ║
║  │   # Combine                                                                      │  ║
║  │   dataset = pd.concat([df3_X, df3_Y], axis=1)                                    │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  FINAL DATASET SHAPE: ~1,320 rows × 882 columns (881 features + 1 target)             ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 7: Feature Selection

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 7: FEATURE SELECTION (VARIANCE THRESHOLD)                                        ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  PURPOSE: Remove uninformative features to improve model performance                   ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ PROBLEM: LOW VARIANCE FEATURES                                                   │  ║
║  │                                                                                  │  ║
║  │ Some fingerprint bits are almost always 0 or almost always 1:                    │  ║
║  │                                                                                  │  ║
║  │   Feature A (LOW VARIANCE - USELESS):                                            │  ║
║  │   ┌───┬───┬───┬───┬───┬───┬───┬───┬───┬───┐                                     │  ║
║  │   │ 0 │ 0 │ 0 │ 0 │ 0 │ 0 │ 0 │ 0 │ 0 │ 0 │  ← All same value!                 │  ║
║  │   └───┴───┴───┴───┴───┴───┴───┴───┴───┴───┘    Can't distinguish molecules      │  ║
║  │                                                                                  │  ║
║  │   Feature B (HIGH VARIANCE - USEFUL):                                            │  ║
║  │   ┌───┬───┬───┬───┬───┬───┬───┬───┬───┬───┐                                     │  ║
║  │   │ 1 │ 0 │ 1 │ 1 │ 0 │ 0 │ 1 │ 0 │ 1 │ 0 │  ← Good variation!                 │  ║
║  │   └───┴───┴───┴───┴───┴───┴───┴───┴───┴───┘    Helps distinguish molecules      │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ SOLUTION: VARIANCE THRESHOLD                                                     │  ║
║  │                                                                                  │  ║
║  │ Remove features where variance < threshold                                       │  ║
║  │                                                                                  │  ║
║  │ THRESHOLD CALCULATION:                                                           │  ║
║  │   threshold = 0.8 × (1 - 0.8) = 0.16                                            │  ║
║  │                                                                                  │  ║
║  │ This removes features that are:                                                  │  ║
║  │   • >80% zeros (mostly absent substructure)                                      │  ║
║  │   • >80% ones (mostly present substructure)                                      │  ║
║  │                                                                                  │  ║
║  │ CODE:                                                                            │  ║
║  │   from sklearn.feature_selection import VarianceThreshold                        │  ║
║  │   selector = VarianceThreshold(threshold=0.16)                                   │  ║
║  │   X = selector.fit_transform(df3_X)                                              │  ║
║  │                                                                                  │  ║
║  │ RESULT:                                                                          │  ║
║  │   Before: 881 features                                                           │  ║
║  │   After:  ~200-300 features (varies based on data)                               │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 8: Train/Test Split

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 8: TRAIN/TEST SPLIT                                                              ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  PURPOSE: Separate data for training and evaluation                                    ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ WHY SPLIT?                                                                       │  ║
║  │                                                                                  │  ║
║  │ • Training set: Model learns patterns from this data                             │  ║
║  │ • Test set: Evaluate how well model generalizes to NEW data                      │  ║
║  │ • Prevents overfitting (memorizing instead of learning)                          │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ SPLIT VISUALIZATION:                                                             │  ║
║  │                                                                                  │  ║
║  │   FULL DATASET (~1,320 samples)                                                  │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │████████████████████████████████████████████████████████████████████████│   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  │                              │                                                   │  ║
║  │                              ▼                                                   │  ║
║  │   ┌────────────────────────────────────────────────────────┐ ┌──────────────┐   │  ║
║  │   │████████████████████████████████████████████████████████│ │██████████████│   │  ║
║  │   │              TRAINING SET (80%)                        │ │  TEST (20%)  │   │  ║
║  │   │              ~1,056 samples                            │ │ ~264 samples │   │  ║
║  │   └────────────────────────────────────────────────────────┘ └──────────────┘   │  ║
║  │                              │                                      │            │  ║
║  │                              ▼                                      ▼            │  ║
║  │                    Model learns here              Model evaluated here           │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  CODE:                                                                                 ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ from sklearn.model_selection import train_test_split                             │  ║
║  │                                                                                  │  ║
║  │ X_train, X_test, Y_train, Y_test = train_test_split(                            │  ║
║  │     X,                    # Features                                             │  ║
║  │     Y,                    # Target (pIC50)                                       │  ║
║  │     test_size=0.2,        # 20% for testing                                      │  ║
║  │     random_state=42       # For reproducibility (optional)                       │  ║
║  │ )                                                                                │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 9: Model Training & Comparison

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 9: MODEL TRAINING & COMPARISON                                                   ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  TOOL USED: LazyPredict (automatically tests many ML models)                          ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ WHAT LAZYPREDICT DOES:                                                           │  ║
║  │                                                                                  │  ║
║  │   ┌─────────────┐                                                                │  ║
║  │   │ Your Data   │                                                                │  ║
║  │   └──────┬──────┘                                                                │  ║
║  │          │                                                                       │  ║
║  │          ▼                                                                       │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │                        LazyPredict                                       │   │  ║
║  │   │  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌─────────────┐        │   │  ║
║  │   │  │   Random    │ │  Gradient   │ │   Ridge     │ │    SVR      │  ...   │   │  ║
║  │   │  │   Forest    │ │  Boosting   │ │ Regression  │ │             │        │   │  ║
║  │   │  └─────────────┘ └─────────────┘ └─────────────┘ └─────────────┘        │   │  ║
║  │   │                                                                          │   │  ║
║  │   │  Trains 30+ models automatically!                                        │   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  │          │                                                                       │  ║
║  │          ▼                                                                       │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │ COMPARISON TABLE (sorted by R² score)                                    │   │  ║
║  │   │                                                                          │   │  ║
║  │   │ Model                    │ R²     │ RMSE   │ Time (s)                    │   │  ║
║  │   │ ─────────────────────────┼────────┼────────┼─────────                    │   │  ║
║  │   │ ExtraTreesRegressor      │ 0.85   │ 0.42   │ 1.2                         │   │  ║
║  │   │ RandomForestRegressor    │ 0.83   │ 0.45   │ 2.1                         │   │  ║
║  │   │ GradientBoostingRegressor│ 0.81   │ 0.48   │ 3.5                         │   │  ║
║  │   │ XGBRegressor             │ 0.80   │ 0.49   │ 1.8                         │   │  ║
║  │   │ ...                      │ ...    │ ...    │ ...                         │   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  CODE:                                                                                 ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ from lazypredict.Supervised import LazyRegressor                                 │  ║
║  │                                                                                  │  ║
║  │ # Create LazyRegressor instance                                                  │  ║
║  │ clf = LazyRegressor(verbose=0, ignore_warnings=True)                             │  ║
║  │                                                                                  │  ║
║  │ # Fit on training data, evaluate on test data                                    │  ║
║  │ models, predictions = clf.fit(X_train, X_test, Y_train, Y_test)                  │  ║
║  │                                                                                  │  ║
║  │ # View results                                                                   │  ║
║  │ print(predictions)  # Sorted by R² score                                         │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ EVALUATION METRICS:                                                              │  ║
║  │                                                                                  │  ║
║  │ • R² (R-squared): How much variance is explained (0-1, higher = better)         │  ║
║  │ • RMSE: Root Mean Squared Error (lower = better)                                 │  ║
║  │ • MAE: Mean Absolute Error (lower = better)                                      │  ║
║  │                                                                                  │  ║
║  │ INTERPRETATION:                                                                  │  ║
║  │   R² = 0.85 means model explains 85% of pIC50 variation                         │  ║
║  │   RMSE = 0.5 means predictions are off by ~0.5 pIC50 units on average           │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
                                          │
                                          ▼
```

### STEP 10: Final Output & Usage

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║  STEP 10: FINAL OUTPUT & PRACTICAL USAGE                                               ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ WHAT YOU GET:                                                                    │  ║
║  │                                                                                  │  ║
║  │ 1. TRAINED MODEL                                                                 │  ║
║  │    • Can predict pIC50 for new molecules                                         │  ║
║  │    • Input: 881 PubChem fingerprint features                                     │  ║
║  │    • Output: Predicted pIC50 value                                               │  ║
║  │                                                                                  │  ║
║  │ 2. MODEL COMPARISON                                                              │  ║
║  │    • Know which algorithm works best for this data                               │  ║
║  │    • Baseline performance metrics                                                │  ║
║  │                                                                                  │  ║
║  │ 3. PROCESSED DATASETS                                                            │  ║
║  │    • Ready for further analysis or different ML approaches                       │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ HOW TO USE FOR NEW MOLECULES:                                                    │  ║
║  │                                                                                  │  ║
║  │   NEW MOLECULE                                                                   │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │ SMILES: CCO...                                                           │   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  │                              │                                                   │  ║
║  │                              ▼                                                   │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │ STEP A: Calculate PubChem Fingerprints (881 features)                    │   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  │                              │                                                   │  ║
║  │                              ▼                                                   │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │ STEP B: Feed to trained model                                            │   │  ║
║  │   │         predicted_pIC50 = model.predict(fingerprints)                    │   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  │                              │                                                   │  ║
║  │                              ▼                                                   │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │ STEP C: Interpret result                                                 │   │  ║
║  │   │                                                                          │   │  ║
║  │   │   pIC50 > 6.0  →  Likely ACTIVE (good drug candidate!)                  │   │  ║
║  │   │   pIC50 5-6    →  INTERMEDIATE (needs optimization)                     │   │  ║
║  │   │   pIC50 < 5.0  →  Likely INACTIVE (probably not useful)                 │   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

---

### Complete Pipeline Summary

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                           COMPLETE PIPELINE SUMMARY                                    ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  [1] ChEMBL Database                                                             │  ║
║  │       │                                                                          │  ║
║  │       ▼                                                                          │  ║
║  │  [2] Raw Excel (1320 compounds, 30+ columns)                                     │  ║
║  │       │                                                                          │  ║
║  │       ▼ Filter IC50, nM units, remove nulls                                      │  ║
║  │  [3] Curated CSV (1320 compounds, 4 columns)                                     │  ║
║  │       │                                                                          │  ║
║  │       ▼ Add bioactivity labels                                                   │  ║
║  │  [4] Labelled CSV (+ class column)                                               │  ║
║  │       │                                                                          │  ║
║  │       ├──────────────────────────────────────┐                                   │  ║
║  │       ▼                                      ▼                                   │  ║
║  │  [5] EDA & Lipinski                    [6] SMILES → Fingerprints                 │  ║
║  │       │                                      │                                   │  ║
║  │       ▼                                      ▼                                   │  ║
║  │  [7] pIC50 conversion                  [8] 881 PubChem features                  │  ║
║  │       │                                      │                                   │  ║
║  │       └──────────────┬───────────────────────┘                                   │  ║
║  │                      ▼                                                           │  ║
║  │  [9] ML-Ready Dataset (881 features + pIC50)                                     │  ║
║  │       │                                                                          │  ║
║  │       ▼ Variance threshold                                                       │  ║
║  │  [10] Feature Selection (~200-300 features)                                      │  ║
║  │       │                                                                          │  ║
║  │       ▼ 80/20 split                                                              │  ║
║  │  [11] Train/Test Split                                                           │  ║
║  │       │                                                                          │  ║
║  │       ▼ LazyPredict                                                              │  ║
║  │  [12] Model Comparison (30+ algorithms)                                          │  ║
║  │       │                                                                          │  ║
║  │       ▼                                                                          │  ║
║  │  [13] Best Model → Predict new molecules!                                        │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

---

## 4. Code Walkthrough

### 4.1 Loading & Cleaning Data

```python
# Load raw data from ChEMBL
df = pd.read_excel('Beta_amyloid A4_protein_active_compounds.xlsx')

# Filter for IC50 measurements only
df_ic50 = df[df['Standard Type'] == 'IC50']

# Keep only nanomolar (nM) units for consistency
df2 = df_ic50[df_ic50['Standard Units'] == 'nM']

# Remove rows with missing values
df3 = df2[~df2['Standard Value'].isnull()]

# Select relevant columns
selection = ['Molecule ChEMBL ID', 'Molecular Weight', 'Smiles', 'Standard Value']
df4 = df3[selection]
```

### 4.2 Labeling Compounds

```python
# Classify based on IC50 thresholds
bioactivity_threshold = []
for ic50 in df4["Standard Value"]:
    if float(ic50) >= 10000:
        bioactivity_threshold.append("inactive")    # Weak drug
    elif float(ic50) <= 1000:
        bioactivity_threshold.append("active")      # Strong drug
    else:
        bioactivity_threshold.append("intermediate") # Medium

# Add class column
df5 = pd.concat([df4, pd.Series(bioactivity_threshold, name='class')], axis=1)
```

### 4.3 Calculating Lipinski Descriptors

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def lipinski(smiles):
    """Calculate drug-likeness descriptors"""
    descriptors = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        descriptors.append({
            'LogP': Descriptors.MolLogP(mol),           # Fat solubility
            'NumHDonors': Lipinski.NumHDonors(mol),     # H-bond donors
            'NumHAcceptors': Lipinski.NumHAcceptors(mol) # H-bond acceptors
        })
    return pd.DataFrame(descriptors)
```

### 4.4 Converting IC50 to pIC50

```python
# Convert nM to M, then take negative log
df["pIC50"] = -np.log10(df["Standard Value"] / 1e9)
```

### 4.5 Feature Selection & Model Training

```python
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from lazypredict.Supervised import LazyRegressor

# Remove low variance features
selector = VarianceThreshold(threshold=0.16)
X = selector.fit_transform(df3_X)

# Split data
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)

# Compare models
clf = LazyRegressor(verbose=0, ignore_warnings=True)
models, predictions = clf.fit(X_train, X_test, Y_train, Y_test)
```

---

## 5. File Descriptions

| File | Description |
|------|-------------|
| `Beta_amyloid A4_protein_active_compounds.xlsx` | Raw data from ChEMBL database |
| `Beta_amyloid_A4_protein_bioactivity_data_curated.csv` | Cleaned data (IC50 only, nM units) |
| `Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv` | + bioactivity class labels |
| `Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv` | ML-ready dataset with fingerprints |
| `molecule.smi` | SMILES strings for fingerprint calculation |
| `descriptors_output.csv` | 881 PubChem fingerprint features |
| `ml-for-alzheimer-s-drug-discovery-in-progress.ipynb` | Main Jupyter notebook |
| `plot_*.pdf` | EDA visualizations |

---

## 6. Key Concepts Explained

### 6.1 SMILES Notation
```
CCO                     → Ethanol (C-C-O)
c1ccccc1                → Benzene (aromatic ring)
CC(=O)O                 → Acetic acid
```

### 6.2 Lipinski's Rule of 5
A drug-like molecule should have:
- MW ≤ 500 Da
- LogP ≤ 5
- H-bond donors ≤ 5
- H-bond acceptors ≤ 10

### 6.3 Machine Learning Models
LazyPredict tests: Random Forest, Gradient Boosting, Ridge, Lasso, SVR, KNN, and many more.

---

## 7. How to Run

```bash
# Install dependencies
pip install pandas numpy rdkit scikit-learn lazypredict seaborn matplotlib

# Run notebook
jupyter notebook ml-for-alzheimer-s-drug-discovery-in-progress.ipynb
```

### Quick Start
```python
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

# Load data
df = pd.read_csv('Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv')
X = df.drop('pIC50', axis=1)
Y = df['pIC50']

# Train
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)
model = RandomForestRegressor(n_estimators=100)
model.fit(X_train, Y_train)

# Predict
predictions = model.predict(X_test)
```
