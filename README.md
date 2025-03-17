# MEA_stacking
An RNA prediction approach from alignments of homologous RNAs in the spirit of RNAalifold from the Vienna RNA package. Here, predict structures using Maximum Expected Accuracy (MEA) with stacking


## Overview
This project implements various algorithms for RNA secondary structure prediction and alignment. It includes implementations of:
- **Nussinov algorithm**  for single RNA folding and traceback.
- **Maximizing Expected Accuracy** for finding the consensus structure of RNA alignments, maximizing Base Pairs and Unpaired probabilities.
- **MEA xpected Accuracy with stacking** for taking into account the stacking process

## Getting Started
### 1. Running the Example Notebook
To quickly explore the main functionalities of the project, run the **Jupyter notebook** `notebook.ipynb`.

## Project Structure
```
/project_root/
    ├── modules/                    # Core Python modules for sequence alignment
    │   ├── simple_Mea.py                # Implements MEA
    │   ├── rna_structure.py             # Implements RNA structure with dot-bracket or base pairs list
    │   ├── nussinov.py                  # Simple nussinov algorithm .
    │   ├── RNA_alifold.py               #Implements RNA alifold
    │   ├── MEA_stacking.py              # Implements MEA with stacking

    ├── notebooks/ 
    │   ├── notebook.ipynb # The Jupyter Notebook

    ├── tests/                  # Some test functions
    │   ├── importdata.py    #imports the test data
    │   ├── compare.py          #compares the different algorithms
    │   ├── Testdata #folder containing  RNA alignements

    ├── README.md               # Project documentation
```

## Dependencies
Ensure the following Python libraries are installed:
```
Biopython
numpy
```
Install them via pip:
```
pip install biopython numpy
```

## Evaluation
The accuracy of the alignment is evaluated on the database in the Testdata folder


