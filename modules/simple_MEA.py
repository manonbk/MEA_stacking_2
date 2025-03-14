"""
This module implements the MEA of an alignement of RNA
"""
import RNA

import sys
import os
parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.nussinov import *
from modules.RNA_structure import *

############################## DEFINING COST FUNCTIONS ###############################################


def list_of_bp_proba_and_struc(list_of_sequences):
    '''Computes, for each structure, the matrix of base pair probabilities. 
    Will be useful then to compute unpaired_expected_accuracy and paired_expected accuracy, the cost function
    in the MEA nussinov algorithm.
    Arguments : 
            list_of_sequences : a list of sequences
    Returns : 
            list of the bpp : list of matrixes of the base pair probabilities
                    list_bpp[i] is the base pair proba matrix of sequence list_of_seq[i]
            list_of_struc a list of the structues associated to the sequences
    '''
    #list of the matrixes
    bpp_list=[]
    #list of the structures
    list_of_struc=[]
    # create model details
    md = RNA.md()
    # activate unique multibranch loop decomposition
    md.uniq_ML = 1
    for struc in list_of_sequences:
        # create fold compound object
        fc = RNA.fold_compound(struc, md)
        # compute MFE
        (ss, mfe) = fc.mfe()
        # rescale Boltzmann factors according to MFE; rescaling avoids numerical problems for long sequences
        fc.exp_params_rescale(mfe)
        # compute partition function to fill DP matrices
        fc.pf()
        bpp = fc.bpp()
        bpp_list.append(bpp)
        list_of_struc.append(ss)
    return bpp_list, list_of_struc



def unpaired_expected_acc(i, list_of_bpp):
    ''' Computes the cost function unpaired_expected_acc for the MEA
    Arg :
        i the int we are computing the cost of
        list_of_bpp : list of the bpp matrixes associated to the structure
    Returns : float : unpaired_expected_acc(i) 
    '''
    #print(i)
    somme=0
    for k in range (len(list_of_bpp)):
        paired_proba=0
        for j in range (len(list_of_bpp[0])):
            #print(i,j,k)
            #print(list_of_bpp[k])
            paired_proba+=list_of_bpp[k][j][i]
            paired_proba+=list_of_bpp[k][i][j]
        somme+= (1-paired_proba)
    return somme
    


def paired_expected_acc(indices,list_of_bpp, gamma):
    ''' Computes the cost function paired_acc
    Arg : list_of_struc : a list of structures, defined as a list of tuples of paired bases
        indices : the tuple of indices we are looking for the cost of
    Returns : int : paired_act(i)
    '''
    i,j=indices
    sum=0
    for k in range(len(list_of_bpp)):
        sum+=list_of_bpp[k][i][j]
        sum+=list_of_bpp[k][j][i] #only one of them is non null
    return 2*gamma*sum


def constfunctionMEA(gamma, list_of_bpp):
    """ 
    Defines the right cost function to use in the nussinov Algorithm, 
    taking into account gamma, and stacking the list of bpp tp not compute it every time
    """
    return lambda indices, list_of_struc : paired_expected_acc(indices, list_of_bpp, gamma)

def initcostMEA(list_of_bpp):
    """
    Defines the right initialisation function to use in the nussinov Algorithm, 
    stacking the list of bpp tp not compute it every time
    """
    return lambda i, list_of_struc : unpaired_expected_acc(i, list_of_bpp)



############################## RUN NUSSINOV FOR MEA ###############################################


def MEA(alignment, gamma=2,m=2):
    print("Building consensus structure of : ")
    print(alignment)
    print("Using the simple MEA with parameter gamma = " + str(gamma))
    print("...")
    list_of_bpp, list_of_struc=list_of_bp_proba_and_struc(alignment)
    print(list_of_struc)
    cost_fun=constfunctionMEA(gamma, list_of_bpp)
    init_cost=initcostMEA(list_of_bpp)
    nussMat=nussinov_matrix(list_of_struc, cost_fun, init_cost, m)
    """
    print("Nussinov Matrix:")
    for lign in nussMat:
        print(lign)
    """
    base_pairs= nussinov_traceback(list_of_struc,nussMat, cost_fun, m)
    return (dot_bracket_string(len(list_of_struc[0]),base_pairs,0))

