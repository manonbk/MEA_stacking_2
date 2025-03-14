import RNA

# The RNA sequence
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.RNA_structure import *
from modules.simple_MEA import *
from modules.MEA_stacking import *
from modules.RNA_Alifold import *


def consensus_structures(alignment, gamma, delta,m):
    result=[]
    result.append(MEA_stacking(alignment, gamma, delta,m))
    result.append(MEA(alignment, gamma,m))
    result.append(alifold_alignment(alignment, m))
    return result

def similarity (struc1, struc2):
    count=0
    for i in range(len(struc1)):
        if struc1[i]=='.' and struc2[i]=='.':
            count+=1
    bp_list1=list_of_base_pair(struc1)
    bp_list2=list_of_base_pair(struc2)
    for tuple in bp_list1:
        if tuple in bp_list2:
            count+=2
    return count/len(struc1)


def levenshtein_distance(s1, s2):
    len_s1, len_s2 = len(s1), len(s2)
    dp = [[0] * (len_s2 + 1) for _ in range(len_s1 + 1)]

    for i in range(len_s1 + 1):
        for j in range(len_s2 + 1):
            if i == 0:
                dp[i][j] = j  # Cost of inserting j characters
            elif j == 0:
                dp[i][j] = i  # Cost of deleting i characters
            else:
                cost = 0 if s1[i - 1] == s2[j - 1] else 1
                dp[i][j] = min(dp[i - 1][j] + 1,      # Deletion
                            dp[i][j - 1] + 1,      # Insertion
                            dp[i - 1][j - 1] + cost)  # Substitution

    return (len_s2-dp[len_s1][len_s2])/len_s2


#print(similarity("(((((((..((((........))))((((((.......))))))...............(((((.......))))))))))))",
#"(((((((..((((........))))(((((.........)))))...............(((((.......))))))))))))"))

def similarity_matrix(alignment, dist_function, gamma=2, delta=1,m=2):
    
    structures= consensus_structures(alignment, gamma, delta ,m)
    print(structures)
    matrix=np.zeros((3,3))
    for i in range (3):
        for j in range(i,3):
            sim_score=dist_function(structures[i], structures[j])
            matrix[i][j]=sim_score
            matrix[j][i]=sim_score
    return( structures,matrix)

def plot_similarity_matrix(matrix, method_names):
    plt.figure(figsize=(8, 6))
    sns.heatmap(matrix, annot=True, cmap="coolwarm", xticklabels=method_names, yticklabels=method_names)
    plt.title("Matrice de Similarité des Structures Consensus")
    plt.xlabel("Méthodes")
    plt.ylabel("Méthodes")
    plt.show()


method_names = ["MEA_stacking", "MEA", "Alifold"]
#sim_matrix = similarity_matrix(example_alignement)
#plot_similarity_matrix(sim_matrix, method_names)