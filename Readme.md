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

## External Dependencies

The following tools/databases must be installed separately.

| Tool            | Purpose                                        |
| --------------- | ---------------------------------------------- |
| HHblits         | Multiple sequence alignment generation         |
| UniRef/UniClust | Sequence database for HHblits                  |
| MMseqs2         | Sequence clustering and fast similarity search |
| CCMpred         | Coevolutionary coupling prediction             |
| ESMFold         | Coevolutionary coupling prediction             |



