# MEA_stacking
An RNA prediction approach from alignments of homologous RNAs in the spirit of RNAalifold from the Vienna RNA package. Here, predict structures using Maximum Expected Accuracy (MEA) with stacking


## Overview
This project implements various algorithms for RNA secondary structure prediction and alignment. It includes implementations of:
- **Nussinov algorithm**  for single RNA folding and traceback.
- **Maximum Total Accuracy** for finding the consensus structure of RNA alignments, maximizing Common Base Pairs and Unpaired Positions.
- **Maximizing Expected Accuracy** for finding the consensus structure of RNA alignments, maximizing Base Pairs and Unpaired probabilities.
- **MEA with Conservation Score.** for taking into account covariations and gaps in the consensus structure.

## Getting Started
### 1. Running the Example Notebook
To quickly explore the main functionalities of the project, run the **Jupyter notebook** `notebookfinal.ipynb`.

## Project Structure
```
/project_root/
    ├── modules/                    # Core Python modules for sequence alignment
    │   ├── map_structure_sequence.py    # Map an RNA structure to a sequence based on an alignment string
    │   ├── MEA_conservation.py          # Implements MEA with conservation score
    │   ├── MTA.py                       # Implements MEA
    │   ├── Mea.py                       # Implements MTA
    │   ├── rna_sequence2.py             # Simple functions to deal with sequence identity of an alignment
    │   ├── rna_structure2.py            # Implements RNA structure with dot-bracket or base pairs list
    │   ├── simplecount.py               # Simple count function
    │   ├── single_rna_structure.py      # Nussinov algorithm for a single rna strand.

    ├── notebooks/ 
    │   ├── finalnotebook.ipynb # The Jupyter Notebook
    │   ├── Testdata #folder containing  RNA alignements

    ├── tests/                  # Some test functions
    │   ├── test_nussinov.py    #imports the test data
    │   ├── compare.py          #compares the different consensus structures

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


