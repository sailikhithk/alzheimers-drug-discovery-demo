# Glossary: Machine Learning for Alzheimer’s Drug Discovery

This glossary defines the technical terms, biological concepts, and data science metrics used throughout the Alzheimer's Drug Discovery QSAR project.

---

## Phase 1: Biological Foundation & Data Mining

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
* **Active:** IC50 $\le$ 1,000 nM (Strong candidates).
* **Inactive:** IC50 $\ge$ 10,000 nM (Weak/Useless candidates).
* **Intermediate:** 1,000 < IC50 < 10,000 nM (Borderline candidates).

### **Missing Data (NA)**
Rows in the dataset where critical values (like `standard_value` or `canonical_smiles`) are undefined. These are dropped during preprocessing.

---

## Phase 3: Feature Engineering (Chemistry to Math)

### **Lipinski's Rule of 5**
A rule of thumb used to evaluate **drug-likeness**. It predicts if a chemical compound with a certain pharmacological or biological activity has properties that would make it a likely orally active drug in humans.
* **Molecular Weight (MW):** Should be < 500 Daltons (Size).
* **LogP:** Should be < 5 (Solubility).
* **H-Bond Donors:** Should be < 5.
* **H-Bond Acceptors:** Should be < 10.

### **LogP (Partition Coefficient)**
A measure of the differential solubility of a compound in two immiscible solvents (octanol and water).
* **High LogP (>5):** Lipophilic (Fat-loving/Greasy). Can cross cell membranes but may get trapped in fat tissue or be toxic.
* **Low LogP (<5):** Hydrophilic (Water-loving). Dissolves easily in blood.

### **pIC50**
The negative logarithm of the IC50 value: $pIC50 = -log_{10}(IC50)$.
* **Why used?**
    1.  **Linearity:** Converts the exponential scale of IC50 into a linear scale suitable for regression.
    2.  **Intuition:** Converts "Lower is Better" to "**Higher is Better**."

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

### **PaDEL-Descriptor**
A software tool (Passive ADME Learning) used to calculate molecular descriptors and fingerprints from SMILES strings.

---

## Phase 5: Exploratory Data Analysis (EDA)

### **Chemical Space**
A multidimensional space defined by the molecular descriptors.
* *Visualization:* We often plot MW vs. LogP to visualize if "Active" compounds cluster in a specific "drug-like" region of this space.

### **Mann-Whitney U Test**
A non-parametric statistical test used to determine whether there is a significant difference between the distributions of two independent groups (Active vs. Inactive).
* *Goal:* To prove that the difference in properties (like MW) between active and inactive drugs is not just due to random chance.

### **P-Value**
The probability that the observed results occurred by chance.
* **p < 0.05:** Statistically significant (The difference is real).

---

## Phase 6: Model Building (Machine Learning)

### **Feature Selection (Variance Threshold)**
A technique to reduce the number of input variables ($X$).
* *Process:* Removing columns (fingerprints) that have the same value (0 or 1) for almost all samples, as they provide no useful information to the model.

### **Train/Test Split**
Dividing the dataset into two subsets:
* **Training Set (80%):** Used to train the model.
* **Testing Set (20%):** Used to evaluate the model's performance on unseen data.
* *Purpose:* To detect and prevent **Overfitting**.

### **LazyPredict**
A Python library that helps build a baseline by running many machine learning models with default parameters on a dataset and ranking them by performance metrics.

### **Random Forest Regressor**
An ensemble learning method that constructs a multitude of decision trees at training time and outputs the mean prediction of the individual trees.
* *Status:* Often the top-performing algorithm for chemical QSAR data.

---

## Phase 7: Evaluation (Grading the AI)

### **$R^2$ (Coefficient of Determination)**
A statistical measure that represents the proportion of the variance for a dependent variable ($y$) that's explained by an independent variable ($X$).
* *Range:* 0 to 1.
* *Goal:* > 0.6 is generally considered good for QSAR models.

### **RMSE (Root Mean Squared Error)**
The standard deviation of the prediction errors (residuals). It tells you how concentrated the data is around the line of best fit.
* *Goal:* A lower value indicates better fit.

### **Scatter Plot (Experimental vs. Predicted)**
A visualization comparing the ground truth values from the lab (X-axis) against the values predicted by the AI (Y-axis).
* *Ideal:* All points falling on a straight diagonal line.