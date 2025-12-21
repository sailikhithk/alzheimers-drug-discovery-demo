# Alzheimer's Drug Discovery using Machine Learning
## A Complete Guide for Data Scientists

A computational drug discovery project that uses QSAR modeling to predict the bioactivity of compounds against Beta-amyloid A4 protein - a key target in Alzheimer's disease.

---

## Table of Contents
1. [Project Overview](#1-project-overview)
2. [The Biology: Understanding Alzheimer's Disease](#2-the-biology-understanding-alzheimers-disease)
3. [The Chemistry: Understanding Drug Molecules](#3-the-chemistry-understanding-drug-molecules)
4. [The Data Science: QSAR and Machine Learning](#4-the-data-science-qsar-and-machine-learning)
5. [Data Pipeline - Detailed](#5-data-pipeline---detailed)
6. [Code Walkthrough](#6-code-walkthrough)
7. [File Descriptions](#7-file-descriptions)
8. [How to Run](#8-how-to-run)

---

## 1. Project Overview

**Goal:** Build ML models to predict if a chemical compound will be effective against Beta-amyloid A4 protein.

**Why it matters:**
- Save years of lab testing
- Reduce drug development costs (typically $2.6 billion per drug!)
- Accelerate finding treatments for Alzheimer's

---

## 2. The Biology: Understanding Alzheimer's Disease

### 2.1 What is Alzheimer's Disease?

Alzheimer's is a progressive neurodegenerative disease that destroys memory and cognitive function.

**Statistics:**
- Causes 60-70% of all dementia cases worldwide
- Affects ~50 million people globally
- 6th leading cause of death in the US
- No cure currently exists

### 2.2 What Causes Alzheimer's? The Amyloid Hypothesis

**The Key Player: Beta-Amyloid A4 Protein (Aβ)**

WHAT IS IT?
- A small protein fragment (36-43 amino acids long)
- Produced when a larger protein (APP) is cut by enzymes
- Normally cleared from the brain
- In Alzheimer's: accumulates and forms toxic clumps


**THE DISEASE PROCESS:**

```
STEP 1: APP (Amyloid Precursor Protein) sits in cell membrane
        ════════════════════════════════════════
                       │
                       ▼ Cut by β-secretase enzyme

STEP 2: Intermediate fragment
        ══════════════════════
                       │
                       ▼ Cut by γ-secretase enzyme

STEP 3: Beta-amyloid (Aβ) fragments released
        ▪▪▪▪▪▪▪▪  ▪▪▪▪▪▪▪▪  ▪▪▪▪▪▪▪▪
                       │
                       ▼ Fragments stick together

STEP 4: Oligomers form (small toxic clusters)
        ╔═══════════════════╗
        ║ ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ║  ← MOST TOXIC FORM!
        ╚═══════════════════╝
                       │
                       ▼ Continue aggregating

STEP 5: Amyloid plaques (large deposits in brain)
        ╔═══════════════════════════════════════╗
        ║ ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ ║
        ╚═══════════════════════════════════════╝
```

**WHY IS THIS BAD?**

1. **SYNAPTIC DYSFUNCTION**
   - Aβ oligomers bind to synapses (connections between neurons)
   - Disrupts signal transmission
   - Memory formation requires synaptic activity

2. **INFLAMMATION**
   - Immune cells (microglia) try to clear Aβ
   - Chronic inflammation damages surrounding tissue

3. **OXIDATIVE STRESS**
   - Aβ generates reactive oxygen species (free radicals)
   - Damages cell membranes, proteins, DNA

4. **NEURONAL DEATH**
   - Neurons eventually die
   - Brain shrinks (atrophy)
   - Cognitive function declines

### 2.3 Why Target Beta-Amyloid A4?

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    WHY BETA-AMYLOID IS A DRUG TARGET                                   ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  THERAPEUTIC STRATEGIES:                                                               ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  STRATEGY 1: REDUCE Aβ PRODUCTION                                               │  ║
║  │  ─────────────────────────────────                                              │  ║
║  │  • Block β-secretase (BACE inhibitors)                                          │  ║
║  │  • Block γ-secretase (γ-secretase inhibitors/modulators)                        │  ║
║  │  • Prevent APP from being cut                                                    │  ║
║  │                                                                                  │  ║
║  │  STRATEGY 2: PREVENT Aβ AGGREGATION                                             │  ║
║  │  ─────────────────────────────────────                                          │  ║
║  │  • Small molecules that bind to Aβ                                              │  ║
║  │  • Prevent oligomer/plaque formation                                            │  ║
║  │  • THIS IS WHAT OUR PROJECT TARGETS! ◄────────────────────────────────────────  │  ║
║  │                                                                                  │  ║
║  │  STRATEGY 3: ENHANCE Aβ CLEARANCE                                               │  ║
║  │  ────────────────────────────────                                               │  ║
║  │  • Antibodies that bind and remove Aβ (immunotherapy)                           │  ║
║  │  • Activate microglia to clear plaques                                          │  ║
║  │  • Recent FDA-approved drugs: Aducanumab, Lecanemab                             │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  OUR PROJECT'S APPROACH:                                                               ║
║  ────────────────────────                                                              ║
║  We're looking for small molecules that can:                                           ║
║  • Bind to Beta-amyloid A4 protein                                                    ║
║  • Inhibit its toxic effects                                                           ║
║  • Potentially prevent aggregation                                                     ║
║                                                                                        ║
║  The ChEMBL database contains compounds that have been TESTED in labs                 ║
║  against Beta-amyloid A4. We use this data to train ML models.                        ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

### 2.4 Understanding IC50: The Key Measurement

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    IC50: THE GOLD STANDARD MEASUREMENT                                 ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  WHAT IS IC50?                                                                         ║
║  ─────────────                                                                         ║
║  IC50 = "Half maximal Inhibitory Concentration"                                        ║
║                                                                                        ║
║  DEFINITION:                                                                           ║
║  The concentration of drug needed to inhibit the target by 50%                        ║
║                                                                                        ║
║  VISUAL EXPLANATION:                                                                   ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  Target Activity (%)                                                             │  ║
║  │  100% ┤ ████████████                                                            │  ║
║  │       │             ████                                                         │  ║
║  │       │                 ████                                                     │  ║
║  │   50% ┤ - - - - - - - - - -████- - - - - - - - ← 50% inhibition                 │  ║
║  │       │                        ████                                              │  ║
║  │       │                            ████                                          │  ║
║  │    0% ┤                                ████████████                              │  ║
║  │       └────────────────────────────────────────────────────────────►            │  ║
║  │         Low                    │                              High               │  ║
║  │                                │                                                 │  ║
║  │                              IC50                                                │  ║
║  │                     (Drug Concentration)                                         │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  KEY INSIGHT: LOWER IC50 = MORE POTENT DRUG                                           ║
║  ───────────────────────────────────────────                                          ║
║                                                                                        ║
║  WHY? Because you need LESS drug to achieve the same effect!                          ║
║                                                                                        ║
║  EXAMPLE:                                                                              ║
║  ┌────────────────────────────────────────────────────────────────────────────────┐   ║
║  │                                                                                 │   ║
║  │  Drug A: IC50 = 10 nM     → Need only 10 nanomolar to inhibit 50%              │   ║
║  │  Drug B: IC50 = 10,000 nM → Need 10,000 nanomolar to inhibit 50%               │   ║
║  │                                                                                 │   ║
║  │  Drug A is 1000x MORE POTENT than Drug B!                                       │   ║
║  │                                                                                 │   ║
║  └────────────────────────────────────────────────────────────────────────────────┘   ║
║                                                                                        ║
║  UNITS EXPLAINED:                                                                      ║
║  ────────────────                                                                      ║
║  • nM = nanomolar = 10⁻⁹ moles per liter (very small!)                               ║
║  • μM = micromolar = 10⁻⁶ moles per liter                                            ║
║  • mM = millimolar = 10⁻³ moles per liter                                            ║
║                                                                                        ║
║  1 μM = 1,000 nM                                                                      ║
║  1 mM = 1,000,000 nM                                                                  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```


### 2.5 Bioactivity Classification: Active vs Inactive

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    BIOACTIVITY CLASSIFICATION EXPLAINED                                ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  WHY CLASSIFY?                                                                         ║
║  ─────────────                                                                         ║
║  • IC50 is a continuous value (0.1 to 1,000,000 nM)                                   ║
║  • For drug discovery, we need to make decisions: "Is this compound worth pursuing?"  ║
║  • Classification gives us clear categories                                            ║
║                                                                                        ║
║  THE THRESHOLDS (Industry Standard):                                                   ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │                        IC50 Scale (logarithmic)                                  │  ║
║  │                                                                                  │  ║
║  │   0.1 nM    1 nM    10 nM   100 nM   1 μM    10 μM   100 μM   1 mM              │  ║
║  │     │        │        │        │       │        │        │       │               │  ║
║  │     ▼        ▼        ▼        ▼       ▼        ▼        ▼       ▼               │  ║
║  │   ══════════════════════════════════════════════════════════════════            │  ║
║  │   │                              │                │                │             │  ║
║  │   │◄─────── ACTIVE ─────────────►│◄─INTERMEDIATE─►│◄── INACTIVE ──►│             │  ║
║  │   │        (< 1,000 nM)          │  (1-10 μM)     │   (> 10 μM)    │             │  ║
║  │   │                              │                │                │             │  ║
║  │                                                                                  │  ║
║  │   1,000 nM = 1 μM               10,000 nM = 10 μM                               │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  WHY THESE SPECIFIC THRESHOLDS?                                                        ║
║  ──────────────────────────────                                                        ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ ACTIVE (IC50 < 1,000 nM = 1 μM)                                                 │  ║
║  │ ───────────────────────────────                                                 │  ║
║  │ REASONING:                                                                       │  ║
║  │ • 1 μM is a common "hit" threshold in drug screening                            │  ║
║  │ • Compounds below this are considered "drug-like"                               │  ║
║  │ • Achievable drug concentrations in the body                                     │  ║
║  │ • Good starting point for optimization                                           │  ║
║  │                                                                                  │  ║
║  │ WHAT IT MEANS:                                                                   │  ║
║  │ • Strong binding to target                                                       │  ║
║  │ • Likely to have therapeutic effect                                              │  ║
║  │ • Worth investing in further development                                         │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ INTERMEDIATE (1,000 - 10,000 nM = 1-10 μM)                                      │  ║
║  │ ──────────────────────────────────────────                                      │  ║
║  │ REASONING:                                                                       │  ║
║  │ • Shows some activity but not strong enough                                      │  ║
║  │ • Might be optimized through chemical modification                               │  ║
║  │ • Could be a "lead" compound for further development                             │  ║
║  │                                                                                  │  ║
║  │ WHAT IT MEANS:                                                                   │  ║
║  │ • Moderate binding                                                               │  ║
║  │ • May need structural optimization                                               │  ║
║  │ • Not immediately useful but has potential                                       │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ INACTIVE (IC50 > 10,000 nM = 10 μM)                                             │  ║
║  │ ─────────────────────────────────────                                           │  ║
║  │ REASONING:                                                                       │  ║
║  │ • 10 μM is a common screening cutoff                                            │  ║
║  │ • Very high concentrations needed for effect                                     │  ║
║  │ • Unlikely to achieve therapeutic levels in body                                 │  ║
║  │ • Would likely cause toxicity before being effective                             │  ║
║  │                                                                                  │  ║
║  │ WHAT IT MEANS:                                                                   │  ║
║  │ • Weak or no binding to target                                                   │  ║
║  │ • Not worth pursuing                                                             │  ║
║  │ • Discard from drug development pipeline                                         │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

---

## 3. The Chemistry: Understanding Drug Molecules

### 3.1 SMILES Notation: The Language of Molecules

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    SMILES: SIMPLIFIED MOLECULAR INPUT LINE ENTRY SYSTEM                ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  WHAT IS SMILES?                                                                       ║
║  ───────────────                                                                       ║
║  A way to represent molecular structure as a TEXT STRING                               ║
║  Like a "language" that computers can read to understand molecules                     ║
║                                                                                        ║
║  WHY USE SMILES?                                                                       ║
║  ───────────────                                                                       ║
║  • Compact: A complex molecule in just a few characters                                ║
║  • Searchable: Can search databases for molecules                                      ║
║  • Computable: Computers can process and analyze                                       ║
║  • Universal: Standard format across chemistry software                                ║
║                                                                                        ║
║  BASIC RULES:                                                                          ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  ATOMS:                                                                          │  ║
║  │  • C = Carbon       • N = Nitrogen     • O = Oxygen                             │  ║
║  │  • S = Sulfur       • P = Phosphorus   • F, Cl, Br, I = Halogens               │  ║
║  │  • Lowercase = aromatic (c, n, o, s)                                            │  ║
║  │                                                                                  │  ║
║  │  BONDS:                                                                          │  ║
║  │  • Single bond: implied (CC = C-C)                                              │  ║
║  │  • Double bond: = (C=O)                                                          │  ║
║  │  • Triple bond: # (C#N)                                                          │  ║
║  │  • Aromatic: lowercase letters                                                   │  ║
║  │                                                                                  │  ║
║  │  BRANCHES:                                                                       │  ║
║  │  • Parentheses for side chains: CC(C)C = isobutane                              │  ║
║  │                                                                                  │  ║
║  │  RINGS:                                                                          │  ║
║  │  • Numbers mark ring closures: c1ccccc1 = benzene                               │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  EXAMPLES:                                                                             ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  SMILES              MOLECULE           STRUCTURE                                │  ║
║  │  ──────────────────────────────────────────────────────────────────────────     │  ║
║  │                                                                                  │  ║
║  │  C                   Methane            H                                        │  ║
║  │                                          │                                       │  ║
║  │                                      H ─ C ─ H                                   │  ║
║  │                                          │                                       │  ║
║  │                                          H                                       │  ║
║  │                                                                                  │  ║
║  │  CCO                 Ethanol            H   H                                    │  ║
║  │                                          │   │                                   │  ║
║  │                                      H ─ C ─ C ─ O ─ H                           │  ║
║  │                                          │   │                                   │  ║
║  │                                          H   H                                   │  ║
║  │                                                                                  │  ║
║  │  c1ccccc1            Benzene                 H                                   │  ║
║  │                                              │                                   │  ║
║  │                                          H ─ C ═ C ─ H                           │  ║
║  │                                              ║   ║                               │  ║
║  │                                          H ─ C ─ C ─ H                           │  ║
║  │                                              │                                   │  ║
║  │                                              H                                   │  ║
║  │                                         (aromatic ring)                          │  ║
║  │                                                                                  │  ║
║  │  CC(=O)O             Acetic acid            O                                    │  ║
║  │                                             ║                                    │  ║
║  │                                      H₃C ─ C ─ O ─ H                             │  ║
║  │                                                                                  │  ║
║  │  CC(=O)Nc1ccc...     Acetaminophen    (complex structure)                       │  ║
║  │                      (Tylenol)                                                   │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  REAL EXAMPLE FROM OUR DATA:                                                           ║
║  ────────────────────────────                                                          ║
║  SMILES: Nc1nc(NCc2ccccc2)c2ccccc2n1                                                  ║
║                                                                                        ║
║  BREAKDOWN:                                                                            ║
║  • Nc1nc(...)c2ccccc2n1 = Two fused aromatic rings with nitrogens                     ║
║  • NCc2ccccc2 = Amino group connected to benzyl group                                 ║
║  • This is a quinazoline derivative (common drug scaffold)                            ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```


### 3.2 Molecular Descriptors: Quantifying Molecules

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    MOLECULAR DESCRIPTORS EXPLAINED                                     ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  WHAT ARE MOLECULAR DESCRIPTORS?                                                       ║
║  ───────────────────────────────                                                       ║
║  Numerical values that describe properties of a molecule                               ║
║  They convert molecular structure into NUMBERS that ML can use                         ║
║                                                                                        ║
║  TYPES OF DESCRIPTORS:                                                                 ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  1. CONSTITUTIONAL DESCRIPTORS                                                   │  ║
║  │     • Count atoms, bonds, rings                                                  │  ║
║  │     • Example: Number of carbons, number of rings                                │  ║
║  │                                                                                  │  ║
║  │  2. PHYSICOCHEMICAL DESCRIPTORS                                                  │  ║
║  │     • Physical/chemical properties                                               │  ║
║  │     • Example: Molecular weight, LogP, polar surface area                        │  ║
║  │                                                                                  │  ║
║  │  3. TOPOLOGICAL DESCRIPTORS                                                      │  ║
║  │     • Based on molecular graph (atoms as nodes, bonds as edges)                  │  ║
║  │     • Example: Wiener index, connectivity indices                                │  ║
║  │                                                                                  │  ║
║  │  4. FINGERPRINTS (Binary descriptors)                                            │  ║
║  │     • Presence/absence of substructures                                          │  ║
║  │     • Example: PubChem fingerprints (881 bits)                                   │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

### 3.3 Lipinski's Rule of 5: Drug-Likeness

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    LIPINSKI'S RULE OF 5: PREDICTING DRUG-LIKENESS                      ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  BACKGROUND:                                                                           ║
║  ───────────                                                                           ║
║  In 1997, Christopher Lipinski at Pfizer analyzed thousands of drugs                  ║
║  He found patterns that predict if a molecule can be an ORAL drug                     ║
║                                                                                        ║
║  THE RULE OF 5:                                                                        ║
║  ──────────────                                                                        ║
║  A molecule is likely to be orally active if it has:                                  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  PROPERTY              │ THRESHOLD │ WHY IT MATTERS                             │  ║
║  │  ──────────────────────┼───────────┼────────────────────────────────────────── │  ║
║  │                        │           │                                            │  ║
║  │  Molecular Weight      │ ≤ 500 Da  │ Larger molecules can't cross membranes    │  ║
║  │                        │           │ easily. They're too big to be absorbed    │  ║
║  │                        │           │ from the gut into the bloodstream.        │  ║
║  │                        │           │                                            │  ║
║  │  LogP                  │ ≤ 5       │ Measures fat-solubility. Too high = drug  │  ║
║  │  (Partition coeff.)    │           │ accumulates in fat tissue, poor           │  ║
║  │                        │           │ distribution. Too low = can't cross       │  ║
║  │                        │           │ cell membranes.                           │  ║
║  │                        │           │                                            │  ║
║  │  H-bond Donors         │ ≤ 5       │ Groups like -OH, -NH. Too many = molecule │  ║
║  │  (NumHDonors)          │           │ is too polar, can't cross membranes.      │  ║
║  │                        │           │                                            │  ║
║  │  H-bond Acceptors      │ ≤ 10      │ Groups like =O, -O-, -N<. Same reason     │  ║
║  │  (NumHAcceptors)       │           │ as donors - affects membrane crossing.    │  ║
║  │                        │           │                                            │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  WHY "RULE OF 5"?                                                                      ║
║  ────────────────                                                                      ║
║  All thresholds are multiples of 5! (500, 5, 5, 10)                                   ║
║  Easy to remember.                                                                     ║
║                                                                                        ║
║  DETAILED EXPLANATIONS:                                                                ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ MOLECULAR WEIGHT (MW)                                                            │  ║
║  │ ─────────────────────                                                            │  ║
║  │                                                                                  │  ║
║  │ WHAT: Total mass of all atoms in the molecule (in Daltons)                      │  ║
║  │                                                                                  │  ║
║  │ HOW CALCULATED:                                                                  │  ║
║  │   MW = Σ (atomic mass × count)                                                  │  ║
║  │   Example: Water (H₂O) = 2(1.008) + 16.00 = 18.02 Da                           │  ║
║  │                                                                                  │  ║
║  │ WHY IT MATTERS FOR DRUGS:                                                        │  ║
║  │                                                                                  │  ║
║  │   Small molecule (MW < 500)          Large molecule (MW > 500)                  │  ║
║  │   ┌─────────────────────┐            ┌─────────────────────────────┐            │  ║
║  │   │         ○           │            │                             │            │  ║
║  │   │    Can easily       │            │    ████████████████████     │            │  ║
║  │   │    pass through     │            │    Too big to fit through   │            │  ║
║  │   │    cell membrane    │            │    membrane pores           │            │  ║
║  │   └─────────────────────┘            └─────────────────────────────┘            │  ║
║  │         ↓ ↓ ↓                                    ✗                              │  ║
║  │   ═══════════════════              ═══════════════════════════════              │  ║
║  │      Cell membrane                        Cell membrane                          │  ║
║  │                                                                                  │  ║
║  │ TYPICAL VALUES IN OUR DATA: 150 - 850 Da                                        │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ LogP (PARTITION COEFFICIENT)                                                     │  ║
║  │ ────────────────────────────                                                     │  ║
║  │                                                                                  │  ║
║  │ WHAT: Measures how much a molecule prefers fat vs water                         │  ║
║  │                                                                                  │  ║
║  │ HOW MEASURED:                                                                    │  ║
║  │   Shake molecule in a mixture of octanol (fat-like) and water                   │  ║
║  │   Measure concentration in each layer                                            │  ║
║  │                                                                                  │  ║
║  │   LogP = log₁₀(concentration in octanol / concentration in water)               │  ║
║  │                                                                                  │  ║
║  │ INTERPRETATION:                                                                  │  ║
║  │                                                                                  │  ║
║  │   LogP < 0     │ Prefers water (hydrophilic)                                    │  ║
║  │   LogP = 0     │ Equal preference                                               │  ║
║  │   LogP > 0     │ Prefers fat (lipophilic)                                       │  ║
║  │   LogP > 5     │ Too fat-soluble (problematic)                                  │  ║
║  │                                                                                  │  ║
║  │ WHY IT MATTERS:                                                                  │  ║
║  │                                                                                  │  ║
║  │   ┌─────────────────────────────────────────────────────────────────────────┐   │  ║
║  │   │                                                                         │   │  ║
║  │   │   Too LOW LogP (< 0)           │   Too HIGH LogP (> 5)                 │   │  ║
║  │   │   ─────────────────            │   ──────────────────                  │   │  ║
║  │   │   • Can't cross membranes      │   • Accumulates in fat tissue        │   │  ║
║  │   │   • Poor absorption            │   • Poor water solubility            │   │  ║
║  │   │   • Stays in blood/water       │   • Hard to formulate                │   │  ║
║  │   │                                │   • Potential toxicity               │   │  ║
║  │   │                                │                                       │   │  ║
║  │   │              IDEAL RANGE: 1 - 3                                        │   │  ║
║  │   │                                                                         │   │  ║
║  │   └─────────────────────────────────────────────────────────────────────────┘   │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ HYDROGEN BOND DONORS & ACCEPTORS                                                 │  ║
║  │ ────────────────────────────────                                                 │  ║
║  │                                                                                  │  ║
║  │ WHAT ARE HYDROGEN BONDS?                                                         │  ║
║  │   Weak attractions between H attached to O/N and another O/N                    │  ║
║  │   Important for molecular interactions                                           │  ║
║  │                                                                                  │  ║
║  │ H-BOND DONOR: Atom that "gives" hydrogen                                        │  ║
║  │   • -OH (hydroxyl)                                                               │  ║
║  │   • -NH₂ (amine)                                                                │  ║
║  │   • -NH- (secondary amine)                                                       │  ║
║  │                                                                                  │  ║
║  │ H-BOND ACCEPTOR: Atom that "receives" hydrogen                                  │  ║
║  │   • =O (carbonyl)                                                                │  ║
║  │   • -O- (ether)                                                                  │  ║
║  │   • -N< (tertiary amine)                                                         │  ║
║  │                                                                                  │  ║
║  │ WHY THEY MATTER:                                                                 │  ║
║  │   • Too many = molecule is too polar                                             │  ║
║  │   • Polar molecules can't cross the lipid (fat) cell membrane                   │  ║
║  │   • They get "stuck" in water                                                    │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```


### 3.4 Molecular Fingerprints: The ML Features

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    MOLECULAR FINGERPRINTS: DEEP DIVE                                   ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  WHAT ARE FINGERPRINTS?                                                                ║
║  ──────────────────────                                                                ║
║  A way to encode molecular structure as a BINARY VECTOR (0s and 1s)                   ║
║  Each position (bit) represents a specific structural feature                          ║
║                                                                                        ║
║  ANALOGY:                                                                              ║
║  ─────────                                                                             ║
║  Think of it like a checklist:                                                         ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  □ Has benzene ring?                    → 1 (yes) or 0 (no)                     │  ║
║  │  □ Has hydroxyl group (-OH)?            → 1 (yes) or 0 (no)                     │  ║
║  │  □ Has nitrogen?                        → 1 (yes) or 0 (no)                     │  ║
║  │  □ Has 5-membered ring?                 → 1 (yes) or 0 (no)                     │  ║
║  │  □ Has carbon-carbon double bond?       → 1 (yes) or 0 (no)                     │  ║
║  │  ... (881 questions for PubChem fingerprints)                                   │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  WHY USE FINGERPRINTS FOR ML?                                                          ║
║  ────────────────────────────                                                          ║
║                                                                                        ║
║  1. FIXED LENGTH                                                                       ║
║     • Every molecule → same number of features (881)                                  ║
║     • ML algorithms need consistent input size                                         ║
║                                                                                        ║
║  2. CAPTURES STRUCTURE                                                                 ║
║     • Similar molecules → similar fingerprints                                        ║
║     • Different molecules → different fingerprints                                    ║
║                                                                                        ║
║  3. INTERPRETABLE                                                                      ║
║     • Each bit has a meaning                                                           ║
║     • Can understand what features matter                                              ║
║                                                                                        ║
║  PUBCHEM FINGERPRINTS (881 bits):                                                      ║
║  ─────────────────────────────────                                                     ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  SECTION 1: ELEMENT COUNTS (Bits 1-115)                                         │  ║
║  │  ──────────────────────────────────────                                         │  ║
║  │  • Bit 1: ≥1 Hydrogen                                                           │  ║
║  │  • Bit 2: ≥2 Hydrogen                                                           │  ║
║  │  • Bit 3: ≥4 Hydrogen                                                           │  ║
║  │  • ...                                                                           │  ║
║  │  • Bit 10: ≥1 Carbon                                                            │  ║
║  │  • Bit 11: ≥2 Carbon                                                            │  ║
║  │  • ...                                                                           │  ║
║  │  • Bit 50: ≥1 Nitrogen                                                          │  ║
║  │  • Bit 60: ≥1 Oxygen                                                            │  ║
║  │  • Bit 70: ≥1 Sulfur                                                            │  ║
║  │  • Bit 80: ≥1 Fluorine                                                          │  ║
║  │  • Bit 90: ≥1 Chlorine                                                          │  ║
║  │  • ...                                                                           │  ║
║  │                                                                                  │  ║
║  │  SECTION 2: RING SYSTEMS (Bits 116-263)                                         │  ║
║  │  ──────────────────────────────────────                                         │  ║
║  │  • Bit 116: Has any ring                                                        │  ║
║  │  • Bit 120: Has 3-membered ring                                                 │  ║
║  │  • Bit 125: Has 4-membered ring                                                 │  ║
║  │  • Bit 130: Has 5-membered ring                                                 │  ║
║  │  • Bit 135: Has 6-membered ring (like benzene)                                  │  ║
║  │  • Bit 140: Has aromatic ring                                                   │  ║
║  │  • Bit 150: Has fused rings                                                     │  ║
║  │  • ...                                                                           │  ║
║  │                                                                                  │  ║
║  │  SECTION 3: ATOM PAIRS (Bits 264-459)                                           │  ║
║  │  ────────────────────────────────────                                           │  ║
║  │  • Bit 264: C-C single bond                                                     │  ║
║  │  • Bit 270: C=C double bond                                                     │  ║
║  │  • Bit 280: C-N bond                                                            │  ║
║  │  • Bit 290: C-O bond                                                            │  ║
║  │  • Bit 300: C=O bond (carbonyl)                                                 │  ║
║  │  • Bit 310: N-H bond                                                            │  ║
║  │  • Bit 320: O-H bond                                                            │  ║
║  │  • ...                                                                           │  ║
║  │                                                                                  │  ║
║  │  SECTION 4: ATOM ENVIRONMENTS (Bits 460-579)                                    │  ║
║  │  ───────────────────────────────────────────                                    │  ║
║  │  • Bit 460: Carbon with 1 neighbor                                              │  ║
║  │  • Bit 465: Carbon with 2 neighbors                                             │  ║
║  │  • Bit 470: Carbon with 3 neighbors                                             │  ║
║  │  • Bit 475: Carbon with 4 neighbors                                             │  ║
║  │  • Bit 480: Nitrogen with 1 neighbor                                            │  ║
║  │  • ...                                                                           │  ║
║  │                                                                                  │  ║
║  │  SECTION 5: SMARTS PATTERNS (Bits 580-881)                                      │  ║
║  │  ─────────────────────────────────────────                                      │  ║
║  │  • Specific substructure patterns                                               │  ║
║  │  • Functional groups                                                             │  ║
║  │  • Common drug scaffolds                                                         │  ║
║  │  • Pharmacophore features                                                        │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  EXAMPLE:                                                                              ║
║  ─────────                                                                             ║
║  Aspirin (CC(=O)Oc1ccccc1C(=O)O)                                                      ║
║                                                                                        ║
║  Fingerprint (first 20 bits):                                                          ║
║  [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, ...]                    ║
║                                                                                        ║
║  Interpretation:                                                                       ║
║  • Bit 1 = 1: Has hydrogen                                                            ║
║  • Bit 10 = 1: Has carbon                                                             ║
║  • Bit 60 = 1: Has oxygen                                                             ║
║  • Bit 135 = 1: Has 6-membered ring                                                   ║
║  • Bit 140 = 1: Has aromatic ring                                                     ║
║  • Bit 300 = 1: Has C=O bond                                                          ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

### 3.5 pIC50: Why We Transform IC50

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    pIC50 TRANSFORMATION: THE MATH AND REASONING                        ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  THE PROBLEM WITH RAW IC50:                                                            ║
║  ──────────────────────────                                                            ║
║                                                                                        ║
║  IC50 values in our data range from 0.3 nM to 300,000 nM                              ║
║  That's a 1,000,000x difference!                                                       ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  RAW IC50 DISTRIBUTION:                                                          │  ║
║  │                                                                                  │  ║
║  │  Count                                                                           │  ║
║  │    │                                                                             │  ║
║  │  █████                                                                           │  ║
║  │  █████                                                                           │  ║
║  │  █████                                                                           │  ║
║  │  █████ █                                                                         │  ║
║  │  █████ █                                                                         │  ║
║  │  █████ █ █                                                                       │  ║
║  │  █████ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █ █   │  ║
║  │  └─────────────────────────────────────────────────────────────────────────►    │  ║
║  │  0.1   1    10   100  1K   10K  100K  300K                                      │  ║
║  │                    IC50 (nM)                                                     │  ║
║  │                                                                                  │  ║
║  │  PROBLEM: Highly skewed! Most values clustered at low end.                      │  ║
║  │           ML algorithms struggle with this distribution.                         │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  THE SOLUTION: LOGARITHMIC TRANSFORMATION                                              ║
║  ────────────────────────────────────────                                              ║
║                                                                                        ║
║  FORMULA:                                                                              ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │   pIC50 = -log₁₀(IC50 in Molar)                                                 │  ║
║  │                                                                                  │  ║
║  │   Step 1: Convert nM to M                                                        │  ║
║  │           IC50_M = IC50_nM ÷ 1,000,000,000                                      │  ║
║  │           IC50_M = IC50_nM × 10⁻⁹                                               │  ║
║  │                                                                                  │  ║
║  │   Step 2: Take negative log                                                      │  ║
║  │           pIC50 = -log₁₀(IC50_M)                                                │  ║
║  │                                                                                  │  ║
║  │   COMBINED:                                                                      │  ║
║  │           pIC50 = -log₁₀(IC50_nM × 10⁻⁹)                                        │  ║
║  │           pIC50 = -log₁₀(IC50_nM) - log₁₀(10⁻⁹)                                │  ║
║  │           pIC50 = -log₁₀(IC50_nM) + 9                                           │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  CONVERSION TABLE:                                                                     ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  IC50 (nM)    │  IC50 (M)      │  pIC50   │  Interpretation                     │  ║
║  │  ─────────────┼────────────────┼──────────┼───────────────────────────────────  │  ║
║  │  0.1          │  10⁻¹⁰         │  10.0    │  Exceptionally potent               │  ║
║  │  1            │  10⁻⁹          │  9.0     │  Extremely potent                   │  ║
║  │  10           │  10⁻⁸          │  8.0     │  Very potent                        │  ║
║  │  100          │  10⁻⁷          │  7.0     │  Potent                             │  ║
║  │  1,000        │  10⁻⁶ (1 μM)   │  6.0     │  Active threshold ◄────────────    │  ║
║  │  10,000       │  10⁻⁵ (10 μM)  │  5.0     │  Inactive threshold ◄──────────    │  ║
║  │  100,000      │  10⁻⁴          │  4.0     │  Weak                               │  ║
║  │  1,000,000    │  10⁻³ (1 mM)   │  3.0     │  Very weak                          │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  BENEFITS OF pIC50:                                                                    ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  1. BETTER DISTRIBUTION                                                          │  ║
║  │     • Transforms skewed data to more normal distribution                         │  ║
║  │     • ML algorithms work better with normal distributions                        │  ║
║  │                                                                                  │  ║
║  │  pIC50 DISTRIBUTION:                                                             │  ║
║  │                                                                                  │  ║
║  │  Count                                                                           │  ║
║  │    │           ████                                                              │  ║
║  │    │         ████████                                                            │  ║
║  │    │       ████████████                                                          │  ║
║  │    │     ████████████████                                                        │  ║
║  │    │   ████████████████████                                                      │  ║
║  │    │ ████████████████████████                                                    │  ║
║  │    └─────────────────────────────────────────────────────────────────────►      │  ║
║  │      3     4     5     6     7     8     9    10                                 │  ║
║  │                      pIC50                                                       │  ║
║  │                                                                                  │  ║
║  │  2. INTUITIVE INTERPRETATION                                                     │  ║
║  │     • Higher pIC50 = More potent (easier to understand!)                        │  ║
║  │     • Each unit increase = 10x more potent                                       │  ║
║  │                                                                                  │  ║
║  │  3. ADDITIVE RELATIONSHIPS                                                       │  ║
║  │     • Differences in pIC50 are meaningful                                        │  ║
║  │     • pIC50 of 7 vs 6 = 10x difference in potency                               │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  CODE:                                                                                 ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  # Convert IC50 (nM) to pIC50                                                   │  ║
║  │  df["IC50_M"] = df["Standard Value"] / 1e9  # nM to M                           │  ║
║  │  df["pIC50"] = -np.log10(df["IC50_M"])                                          │  ║
║  │                                                                                  │  ║
║  │  # Or in one line:                                                               │  ║
║  │  df["pIC50"] = -np.log10(df["Standard Value"] * 1e-9)                           │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

---

## 4. The Data Science: QSAR and Machine Learning

### 4.1 What is QSAR?

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    QSAR: QUANTITATIVE STRUCTURE-ACTIVITY RELATIONSHIP                  ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  THE CORE IDEA:                                                                        ║
║  ──────────────                                                                        ║
║  A molecule's STRUCTURE determines its ACTIVITY                                        ║
║                                                                                        ║
║  If we can mathematically describe structure, we can PREDICT activity!                ║
║                                                                                        ║
║  THE QSAR EQUATION:                                                                    ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │   Activity = f(Structural Features)                                              │  ║
║  │                                                                                  │  ║
║  │   pIC50 = f(Fingerprint₁, Fingerprint₂, ..., Fingerprint₈₈₁)                    │  ║
║  │                                                                                  │  ║
║  │   Where f() is learned by machine learning                                       │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  THE WORKFLOW:                                                                         ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │   ┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐   │  ║
║  │   │  Molecule   │ ──► │  Features   │ ──► │  ML Model   │ ──► │  Predicted  │   │  ║
║  │   │  (SMILES)   │     │(Fingerprint)│     │  (trained)  │     │   pIC50     │   │  ║
║  │   └─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘   │  ║
║  │                                                                                  │  ║
║  │   "CCO..."           [1,0,1,1,0...]      RandomForest         6.5               │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  WHY QSAR WORKS:                                                                       ║
║  ───────────────                                                                       ║
║                                                                                        ║
║  1. SIMILAR STRUCTURES → SIMILAR ACTIVITIES                                           ║
║     • Molecules with similar shapes bind similarly to targets                         ║
║     • Small changes in structure → small changes in activity                          ║
║                                                                                        ║
║  2. STRUCTURE ENCODES PROPERTIES                                                       ║
║     • Functional groups determine chemical behavior                                    ║
║     • Shape determines how molecule fits into target                                   ║
║                                                                                        ║
║  3. PATTERNS ARE LEARNABLE                                                             ║
║     • ML can find complex patterns in high-dimensional data                           ║
║     • Fingerprints capture relevant structural information                            ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```


### 4.2 The Machine Learning Pipeline

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    THE COMPLETE ML PIPELINE                                            ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  STEP 1: DATA COLLECTION                                                         │  ║
║  │  ════════════════════════                                                        │  ║
║  │  • Source: ChEMBL database                                                       │  ║
║  │  • ~1,320 compounds tested against Beta-amyloid A4                              │  ║
║  │  • Each compound has: SMILES, IC50, molecular properties                        │  ║
║  │                                                                                  │  ║
║  │                              ↓                                                   │  ║
║  │                                                                                  │  ║
║  │  STEP 2: DATA CLEANING                                                           │  ║
║  │  ═══════════════════════                                                         │  ║
║  │  • Filter: Only IC50 measurements                                                │  ║
║  │  • Standardize: Only nM units                                                    │  ║
║  │  • Remove: Missing values                                                        │  ║
║  │                                                                                  │  ║
║  │                              ↓                                                   │  ║
║  │                                                                                  │  ║
║  │  STEP 3: FEATURE ENGINEERING                                                     │  ║
║  │  ═══════════════════════════                                                     │  ║
║  │  • Calculate: PubChem fingerprints (881 features)                               │  ║
║  │  • Transform: IC50 → pIC50                                                      │  ║
║  │  • Add: Lipinski descriptors (for EDA)                                          │  ║
║  │                                                                                  │  ║
║  │                              ↓                                                   │  ║
║  │                                                                                  │  ║
║  │  STEP 4: FEATURE SELECTION                                                       │  ║
║  │  ═════════════════════════                                                       │  ║
║  │  • Remove: Low variance features                                                 │  ║
║  │  • Keep: ~200-300 informative features                                          │  ║
║  │                                                                                  │  ║
║  │                              ↓                                                   │  ║
║  │                                                                                  │  ║
║  │  STEP 5: TRAIN/TEST SPLIT                                                        │  ║
║  │  ════════════════════════                                                        │  ║
║  │  • Training: 80% (~1,056 compounds)                                             │  ║
║  │  • Testing: 20% (~264 compounds)                                                │  ║
║  │                                                                                  │  ║
║  │                              ↓                                                   │  ║
║  │                                                                                  │  ║
║  │  STEP 6: MODEL TRAINING                                                          │  ║
║  │  ══════════════════════                                                          │  ║
║  │  • LazyPredict: Tests 30+ algorithms automatically                              │  ║
║  │  • Each model learns: fingerprints → pIC50                                      │  ║
║  │                                                                                  │  ║
║  │                              ↓                                                   │  ║
║  │                                                                                  │  ║
║  │  STEP 7: EVALUATION                                                              │  ║
║  │  ══════════════════                                                              │  ║
║  │  • Metrics: R², RMSE, MAE                                                       │  ║
║  │  • Compare: All models on test set                                              │  ║
║  │  • Select: Best performing model                                                 │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

### 4.3 Understanding the ML Models

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    ML MODELS EXPLAINED                                                 ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  This is a REGRESSION problem (predicting continuous pIC50 values)                    ║
║                                                                                        ║
║  MODELS TESTED BY LAZYPREDICT:                                                         ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ RANDOM FOREST REGRESSOR                                                          │  ║
║  │ ───────────────────────                                                          │  ║
║  │                                                                                  │  ║
║  │ HOW IT WORKS:                                                                    │  ║
║  │ • Creates many decision trees (a "forest")                                       │  ║
║  │ • Each tree trained on random subset of data                                     │  ║
║  │ • Final prediction = average of all trees                                        │  ║
║  │                                                                                  │  ║
║  │ WHY IT'S GOOD FOR DRUG DISCOVERY:                                               │  ║
║  │ • Handles high-dimensional data (881 features)                                   │  ║
║  │ • Captures non-linear relationships                                              │  ║
║  │ • Robust to noise                                                                │  ║
║  │ • Can identify important features                                                │  ║
║  │                                                                                  │  ║
║  │ VISUAL:                                                                          │  ║
║  │                                                                                  │  ║
║  │   Tree 1        Tree 2        Tree 3       ...      Tree N                      │  ║
║  │     │             │             │                     │                          │  ║
║  │     ▼             ▼             ▼                     ▼                          │  ║
║  │   pred=6.2      pred=6.5      pred=6.1             pred=6.4                     │  ║
║  │     │             │             │                     │                          │  ║
║  │     └─────────────┴─────────────┴─────────────────────┘                          │  ║
║  │                           │                                                      │  ║
║  │                           ▼                                                      │  ║
║  │                    Average = 6.3                                                 │  ║
║  │                  (Final prediction)                                              │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ GRADIENT BOOSTING REGRESSOR                                                      │  ║
║  │ ───────────────────────────                                                      │  ║
║  │                                                                                  │  ║
║  │ HOW IT WORKS:                                                                    │  ║
║  │ • Builds trees SEQUENTIALLY                                                      │  ║
║  │ • Each tree corrects errors of previous trees                                    │  ║
║  │ • Final prediction = sum of all tree predictions                                 │  ║
║  │                                                                                  │  ║
║  │ WHY IT'S GOOD:                                                                   │  ║
║  │ • Often achieves best accuracy                                                   │  ║
║  │ • Learns complex patterns                                                        │  ║
║  │ • Handles feature interactions                                                   │  ║
║  │                                                                                  │  ║
║  │ VISUAL:                                                                          │  ║
║  │                                                                                  │  ║
║  │   Tree 1 ──► Residuals ──► Tree 2 ──► Residuals ──► Tree 3 ──► ...             │  ║
║  │     │                        │                        │                          │  ║
║  │   pred=5.0                 pred=+1.0                pred=+0.3                   │  ║
║  │     │                        │                        │                          │  ║
║  │     └────────────────────────┴────────────────────────┘                          │  ║
║  │                           │                                                      │  ║
║  │                           ▼                                                      │  ║
║  │                    Sum = 6.3                                                     │  ║
║  │                  (Final prediction)                                              │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ OTHER MODELS TESTED:                                                             │  ║
║  │                                                                                  │  ║
║  │ • Ridge Regression: Linear model with L2 regularization                         │  ║
║  │ • Lasso Regression: Linear model with L1 regularization (feature selection)     │  ║
║  │ • ElasticNet: Combination of Ridge and Lasso                                    │  ║
║  │ • SVR: Support Vector Regression (finds optimal hyperplane)                     │  ║
║  │ • KNN: K-Nearest Neighbors (predicts based on similar molecules)                │  ║
║  │ • XGBoost: Optimized gradient boosting                                          │  ║
║  │ • LightGBM: Fast gradient boosting                                              │  ║
║  │ • AdaBoost: Adaptive boosting                                                   │  ║
║  │ • ExtraTrees: Extremely randomized trees                                        │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

### 4.4 Evaluation Metrics

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                    EVALUATION METRICS EXPLAINED                                        ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ R² (R-SQUARED) - COEFFICIENT OF DETERMINATION                                    │  ║
║  │ ─────────────────────────────────────────────                                    │  ║
║  │                                                                                  │  ║
║  │ WHAT IT MEASURES:                                                                │  ║
║  │ How much of the variance in pIC50 is explained by the model                     │  ║
║  │                                                                                  │  ║
║  │ FORMULA:                                                                         │  ║
║  │ R² = 1 - (Sum of Squared Residuals / Total Sum of Squares)                      │  ║
║  │ R² = 1 - Σ(actual - predicted)² / Σ(actual - mean)²                            │  ║
║  │                                                                                  │  ║
║  │ INTERPRETATION:                                                                  │  ║
║  │ • R² = 1.0: Perfect predictions                                                 │  ║
║  │ • R² = 0.8: Model explains 80% of variance (good!)                              │  ║
║  │ • R² = 0.5: Model explains 50% of variance (moderate)                           │  ║
║  │ • R² = 0.0: Model is no better than predicting the mean                         │  ║
║  │ • R² < 0.0: Model is worse than predicting the mean                             │  ║
║  │                                                                                  │  ║
║  │ FOR DRUG DISCOVERY:                                                              │  ║
║  │ • R² > 0.7 is generally considered good                                         │  ║
║  │ • R² > 0.8 is very good                                                         │  ║
║  │ • R² > 0.9 is excellent (rare in real-world data)                               │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ RMSE (ROOT MEAN SQUARED ERROR)                                                   │  ║
║  │ ──────────────────────────────                                                   │  ║
║  │                                                                                  │  ║
║  │ WHAT IT MEASURES:                                                                │  ║
║  │ Average magnitude of prediction errors (in pIC50 units)                         │  ║
║  │                                                                                  │  ║
║  │ FORMULA:                                                                         │  ║
║  │ RMSE = √(Σ(actual - predicted)² / n)                                            │  ║
║  │                                                                                  │  ║
║  │ INTERPRETATION:                                                                  │  ║
║  │ • RMSE = 0.5: Predictions are off by ~0.5 pIC50 units on average               │  ║
║  │ • RMSE = 1.0: Predictions are off by ~1.0 pIC50 units on average               │  ║
║  │                                                                                  │  ║
║  │ WHAT DOES 0.5 pIC50 ERROR MEAN?                                                 │  ║
║  │ • 0.5 pIC50 ≈ 3x error in IC50                                                  │  ║
║  │ • If true IC50 = 100 nM, prediction might be 30-300 nM                          │  ║
║  │ • Still useful for ranking compounds!                                            │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │ MAE (MEAN ABSOLUTE ERROR)                                                        │  ║
║  │ ─────────────────────────                                                        │  ║
║  │                                                                                  │  ║
║  │ WHAT IT MEASURES:                                                                │  ║
║  │ Average absolute difference between predicted and actual values                 │  ║
║  │                                                                                  │  ║
║  │ FORMULA:                                                                         │  ║
║  │ MAE = Σ|actual - predicted| / n                                                 │  ║
║  │                                                                                  │  ║
║  │ DIFFERENCE FROM RMSE:                                                            │  ║
║  │ • MAE treats all errors equally                                                  │  ║
║  │ • RMSE penalizes large errors more (due to squaring)                            │  ║
║  │ • MAE is more interpretable                                                      │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

---

## 5. Data Pipeline - Detailed

### Complete Pipeline Diagram

```
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                           COMPLETE DATA PIPELINE                                       ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║                                                                                        ║
║  ┌─────────────────────────────────────────────────────────────────────────────────┐  ║
║  │                                                                                  │  ║
║  │  [1] ChEMBL Database                                                             │  ║
║  │       │  • World's largest bioactivity database                                  │  ║
║  │       │  • Contains millions of compound-target interactions                     │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [2] Raw Excel (~1320 compounds, 30+ columns)                                    │  ║
║  │       │  • Beta_amyloid A4_protein_active_compounds.xlsx                        │  ║
║  │       │  • Contains: SMILES, IC50, Ki, EC50, MW, etc.                           │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [3] Filter & Clean                                                              │  ║
║  │       │  • Keep only IC50 measurements (most common)                            │  ║
║  │       │  • Keep only nM units (standardize)                                     │  ║
║  │       │  • Remove null values                                                    │  ║
║  │       │  • Select 4 columns: ID, MW, SMILES, IC50                               │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [4] Curated CSV                                                                 │  ║
║  │       │  • Beta_amyloid_A4_protein_bioactivity_data_curated.csv                 │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [5] Add Bioactivity Labels                                                      │  ║
║  │       │  • IC50 < 1000 nM → "active"                                            │  ║
║  │       │  • 1000-10000 nM → "intermediate"                                       │  ║
║  │       │  • IC50 > 10000 nM → "inactive"                                         │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [6] Labelled CSV                                                                │  ║
║  │       │  • Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv    │  ║
║  │       │                                                                          │  ║
║  │       ├──────────────────────────────────────┐                                   │  ║
║  │       ▼                                      ▼                                   │  ║
║  │                                                                                  │  ║
║  │  [7] EDA & Lipinski                    [8] Fingerprint Calculation               │  ║
║  │       │  • Calculate LogP, HBD, HBA         │  • SMILES → molecule.smi          │  ║
║  │       │  • Generate plots                   │  • PaDEL → 881 PubChem FPs        │  ║
║  │       │  • Analyze distributions            │  • Output: descriptors_output.csv │  ║
║  │       ▼                                      ▼                                   │  ║
║  │                                                                                  │  ║
║  │  [9] Convert IC50 → pIC50                                                       │  ║
║  │       │  • pIC50 = -log10(IC50 in M)                                            │  ║
║  │       │  • Better distribution for ML                                            │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [10] Combine Features + Target                                                  │  ║
║  │       │  • X = 881 fingerprint columns                                          │  ║
║  │       │  • Y = pIC50 column                                                     │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [11] ML-Ready Dataset                                                           │  ║
║  │       │  • Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [12] Variance Threshold Feature Selection                                       │  ║
║  │       │  • Remove features with variance < 0.16                                 │  ║
║  │       │  • ~200-300 features remain                                             │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [13] Train/Test Split (80/20)                                                   │  ║
║  │       │  • Training: ~1,056 compounds                                           │  ║
║  │       │  • Testing: ~264 compounds                                              │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [14] LazyPredict Model Comparison                                               │  ║
║  │       │  • Tests 30+ regression algorithms                                      │  ║
║  │       │  • Ranks by R², RMSE                                                    │  ║
║  │       ▼                                                                          │  ║
║  │                                                                                  │  ║
║  │  [15] Best Model → Ready for Predictions!                                        │  ║
║  │                                                                                  │  ║
║  └─────────────────────────────────────────────────────────────────────────────────┘  ║
║                                                                                        ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
```

---

## 6. Code Walkthrough

### 6.1 Data Loading and Cleaning

```python
import pandas as pd
import numpy as np

# Load raw data
df = pd.read_excel('Beta_amyloid A4_protein_active_compounds.xlsx')

# Filter for IC50 only (most common measurement type)
df_ic50 = df[df['Standard Type'] == 'IC50']

# Standardize units to nM
df2 = df_ic50[df_ic50['Standard Units'] == 'nM']

# Remove missing values
df3 = df2[~df2['Standard Value'].isnull()]

# Select relevant columns
df4 = df3[['Molecule ChEMBL ID', 'Molecular Weight', 'Smiles', 'Standard Value']]
```

### 6.2 Bioactivity Labeling

```python
# Classify compounds based on IC50
labels = []
for ic50 in df4["Standard Value"]:
    if ic50 >= 10000:
        labels.append("inactive")
    elif ic50 <= 1000:
        labels.append("active")
    else:
        labels.append("intermediate")

df5 = pd.concat([df4, pd.Series(labels, name='class')], axis=1)
```

### 6.3 pIC50 Conversion

```python
# Convert IC50 (nM) to pIC50
df5["pIC50"] = -np.log10(df5["Standard Value"] * 1e-9)
```

### 6.4 Feature Selection and Model Training

```python
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from lazypredict.Supervised import LazyRegressor

# Load fingerprints
X = pd.read_csv('descriptors_output.csv').drop('Name', axis=1)
Y = df5['pIC50']

# Remove low variance features
selector = VarianceThreshold(threshold=0.16)
X_selected = selector.fit_transform(X)

# Split data
X_train, X_test, Y_train, Y_test = train_test_split(X_selected, Y, test_size=0.2)

# Compare models
clf = LazyRegressor(verbose=0, ignore_warnings=True)
models, predictions = clf.fit(X_train, X_test, Y_train, Y_test)
print(predictions)
```

---

## 7. File Descriptions

| File | Description |
|------|-------------|
| `Beta_amyloid A4_protein_active_compounds.xlsx` | Raw ChEMBL data |
| `Beta_amyloid_A4_protein_bioactivity_data_curated.csv` | Cleaned data |
| `Beta_amyloid_A4_protein_bioactivity_data_curated_and_labelled.csv` | + class labels |
| `Beta_amyloid_A4_protein_bioactivity_data_3class_pIC50_pubchem_fp.csv` | ML-ready dataset |
| `molecule.smi` | SMILES for fingerprint calculation |
| `descriptors_output.csv` | 881 PubChem fingerprints |
| `ml-for-alzheimer-s-drug-discovery-in-progress.ipynb` | Main notebook |
| `plot_*.pdf` | EDA visualizations |

---

## 8. How to Run

```bash
# Install dependencies
pip install pandas numpy rdkit scikit-learn lazypredict seaborn matplotlib

# Run notebook
jupyter notebook ml-for-alzheimer-s-drug-discovery-in-progress.ipynb
```

---

## Summary

This project demonstrates computational drug discovery:

1. **Biology**: Target Beta-amyloid A4 protein involved in Alzheimer's
2. **Chemistry**: Represent molecules as SMILES → fingerprints
3. **Data Science**: Use QSAR + ML to predict bioactivity (pIC50)
4. **Result**: Model that can screen new molecules for potential Alzheimer's drugs
