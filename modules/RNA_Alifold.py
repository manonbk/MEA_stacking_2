"""This module uses Alifold to create a consensus alignment"""

import RNA
# The RNA sequence

import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

def alifold_alignment(alignment):
    """
    Predicts the structure of a consensus RNA alignment using ViennaRNA's Alifold
    Arguments:
        alignment : list of string. A list of aligned sequences
    Returns :
        ss : str : the string of the dot bracket structure
    """
    print("Building consensus structure of : ")
    print(alignment)
    print("Using RNA Alifold")
    # create a new model details structure
    md = RNA.md()
    # optionally one could change some parameters
    # md.temperature = 25.0 # 25 Deg Celcius
    # md.dangles = 1 # keep default 2 for compatibility with partition folding
    # create a fold compound
    fc = RNA.fold_compound(alignment, md)
    # predict the  "Alifold" Minmum Free Energy and the corresponding secondary structure
    (ss, mfe) = fc.mfe()
    conservation_score = fc.eval_covar_structure(ss)
    print("%s\n%s [ %6.2f, %6.2f ]\n" % ('\n'.join(alignment), ss, mfe, conservation_score))
    return ss