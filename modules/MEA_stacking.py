""" This module completes a nussinov algorithm adapted to the MEA with stacking
Predicts the structure of a single rna"""

import RNA

import sys
import os
parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.simple_MEA import *


############################THE NEW NUSSNOV FILLING AND TRACE BACK############################################

def stacking_matrices(list_of_structures, proba_function, init_cost, gamma, delta, m):
    """
    Constructs a two Matrices, F : minimizing the probability function and C : 
    Args : list_of_structures : list of the aligned structures (strings)
            proba_function : function : (i,j) (indices) -> float computes the sum of proba pij for each structure
            init_cost : function : int (indice) -> float computes the prize of being unpaired to initialize the matrix
                        (same as in simple MEA)
            gamma : float, the base pair contribution
            delta : float, the stacking contribution
            m : int, minimal loop length
    Returns : tuple of matrices (list of list) of nussinov:
            where F[i][j] = Best expected accuracy score for the subsequence [i,j][i,j].
            and C[i][j] =Best expected accuracy score for a closed structure between i and j.
    """
    n=len(list_of_structures[0])
    F = [[0 for i in range (n)]for j in range(n)]
    C = [[0 for i in range (n)]for j in range(n)]

    #initialize for all i-j<=m
    for k in range(m+1):
        i=0
        j=i+k
        while (i<n and j<n):
            if j==i:
                F[i][j]=init_cost(i)
            else :
                F[i][j]=F[i][j-1]+init_cost(j)
            i+=1
            j+=1

    #complete matrix
    for k in range(m, n):
        #indexes managment to fill matrix in diagonals
        i=0
        j=i+k
        while (i<n and j<n):

            #First fill C :
            if (i+1<n and j-1>=0):
                C[i][j]=gamma*proba_function((i,j)) + F[i+1][j-1] + delta*proba_function((i,j))*proba_function((i+1,j-1))
            else :
                C[i][j]=gamma*proba_function((i,j))

            if abs(i-j)>=2:
                alpha1 = max(C[i][k] + F[k+1][j]
                             for k in range(i+1, j))
            else :
                alpha1=0
            
            if i+1<n :
                alpha0=F[i+1][j] + init_cost(i)
            else :
                alpha0=init_cost(i)

    
            F[i][j]=max(alpha0, alpha1, C[i][j])
            i+=1
            j+=1


    return F,C



def stacking_traceback_rec(indices, F, C, init_cost, m):
    """
    The recursion for nussinov traceback: 
    Args : 
            F,C the matrices of the stacking matrices algorithm
            indices : the indices we are looking at
            m : int, minimal loop length
    Returns : the list of base pairs that took us to the optimal solution at step i,j = indices
    """
    previous_bases=[]
    i,j=indices
    found=False
    if j-i<=m:
        return []
    else :
        n=len(F[i])
        cost=F[i][j]

        #Find the previous bases

        #case i is paired to some k betweek i+1 and j-1
        for k in range(i+1,j):
            if cost== (C[i][k] + F[k+1][j]) and not found:
                found=True
                previous_bases=stacking_traceback_rec((i,k),F,C, init_cost,m) + stacking_traceback_rec((k+1,j),F,C, init_cost,m)
                previous_bases.append((i,k))

        
        #case i is unpaired
        if i+1<n and not found :
            if cost==F[i+1][j] + init_cost(i):
                previous_bases=stacking_traceback_rec((i+1,j),F,C, init_cost,m)
                found=True
        
        #case i and j are paired
        if not found:
            if cost== C[i][j]:
                previous_bases=stacking_traceback_rec((i+1,j-1),F,C, init_cost,m)
                previous_bases.append((i,j))
                found=True
        
        if not found:
            print("not found", i, j)
        
        return previous_bases
    

def stacking_traceback(F,C,init_cost,m):
    """"
    The traceback algorithm for the nussinov matrix
    Args : list_of_structures : list of the aligned structures
        const_function function -> tuple of indices, list of structures -> int returns the cost of adding a base pair ij with regard to the list of aligned structures
        init_cost : a cost fonction to initialize the matrix
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    n=len(F)
    return stacking_traceback_rec((0,n-1), F, C, init_cost,m)




############################THE ADAPTED COST FONCTIONS#################################################

#Imports from MEA some functions already comptuted, and adapt it.

def initcostMEA_stacking(list_of_bpp):
    """
    Defines the right initialisation function to use in the nussinov Algorithm, 
    stacking the list of bpp tp not compute it every time
    """
    return lambda i : unpaired_expected_acc(i, list_of_bpp)
    

def paired_expected_acc2(indices,list_of_bpp):
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
    return 2*sum


def costfunctionMEA_stacking(list_of_bpp):
    """ 
    Defines the right cost function to use in the nussinov Algorithm, 
    taking into account stacking the list of bpp tp not compute it every time
    """
    return lambda indices : paired_expected_acc2(indices, list_of_bpp)



#########################THE FINAL MEA WITH STACKING ALGORITHM############################################


def MEA_stacking(alignment, gamma=2, delta= 1/4,m=2):
    print("Building consensus structure of : ")
    print(alignment)
    print("Using the MEA_stacking with parameter gamma =" + str(gamma) + ", delta=" + str(delta))
    print("...")
    list_of_bpp, list_of_struc=list_of_bp_proba_and_struc(alignment)
    print(list_of_struc)
    cost_fun=costfunctionMEA_stacking(list_of_bpp)
    init_cost=initcostMEA_stacking(list_of_bpp)
    F,C=stacking_matrices(list_of_struc, cost_fun, init_cost, gamma, delta, m)
    """
    print("F Matrix")
    for lign in F:
        print(lign)
    print("C Matrix")
    for lign in C:
        print(lign)
    """
    base_pairs= stacking_traceback(F,C,init_cost,m)
    return (dot_bracket_string(len(list_of_struc[0]),base_pairs,0))