# ICM-MD  
### Integrating TM-Specific Features and MD-Derived Structures for Accurate Prediction of Inter-Chain Contacts in α-Helical Transmembrane Homodimers

---

## Overview

ICM-MD is a deep learning framework for predicting **inter-chain residue contacts** in **α-helical transmembrane (TM) homodimers** by integrating:

- TM-specific sequence and physicochemical features
- Evolutionary coupling information
- Protein language model embeddings
- Molecular dynamics (MD)-derived structural features

The framework combines biologically informed descriptors with modern machine learning to improve residue-pair contact prediction in membrane protein complexes.

---

# Repository Structure

```text
ICM-MD/
│
├── data/
│   ├── raw/
│   ├── processed/
│   ├── msa/
│   ├── embeddings/
│   ├── ccmpred/
│   ├── esmfold/
│   └── md_features/
│
├── models/
│   ├── checkpoints/
│   └── trained_models/
│
├── scripts/
│   ├── feature_extraction/
│   ├── preprocessing/
│   ├── training/
│   ├── evaluation/
│   └── utilities/
│
├── notebooks/
│
├── configs/
│
├── environment/
│   └── icm_md.yml
│
├── results/
│
├── README.md
└── LICENSE


