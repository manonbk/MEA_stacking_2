import RNA

import sys
import os
parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.simple_MEA import *
from modules.MEA_stacking import *


##################TESTS###########################################

single_rna_example=["acguacggccauauccgagacacgcguaccggaacccauuccgaauuccga-agucaagcguccgcgag-uuggguuaguaaucuggugaaagaucacaggcgaacccccaa--u-gcuguacguc"]


list_of_bpp, list_of_struc=(list_of_bp_proba_and_struc(single_rna_example))

print(list_of_struc)


"""
func1= initcostMEA_stacking(list_of_bpp)
print(func1(3))

func2= costfunctionMEA_stacking(list_of_bpp)
print(func2((37,43)))
"""


MEA_stacking(single_rna_example)