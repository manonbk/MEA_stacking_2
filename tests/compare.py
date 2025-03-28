import RNA

# The RNA sequence
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd()))
parent_dir = os.path.join(parent_dir, "MEA_stacking")
sys.path.append(parent_dir)
#print(parent_dir)

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

#####start Madi's modif#######################
#Matthew's correlation coefficient (MCC)
from sklearn.metrics import matthews_corrcoef

def dot_bracket_to_pairs(dot_bracket):
    """Converts dot-bracket notation to a set of base pair tuples (i, j)."""
    stack = []
    pairs = set()
    
    for i, c in enumerate(dot_bracket):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                j = stack.pop()
                pairs.add((j, i))
    return pairs


def adjusted_mcc(ref_structure, pred_structure):
    """Computes an adjusted MCC where extra compatible pairs are ignored."""
    ref_pairs = dot_bracket_to_pairs(ref_structure)
    pred_pairs = dot_bracket_to_pairs(pred_structure)

    # True Positives (TP): Correctly predicted base pairs
    TP = len(ref_pairs & pred_pairs)

    # False Negatives (FN): Missing reference pairs
    FN = len(ref_pairs - pred_pairs)

    # False Positives (FP): Predicted pairs that contradict the reference
    FP = sum(1 for pair in pred_pairs if pair not in ref_pairs and not is_compatible(pair, ref_pairs))
    print(TP,FN,FP)

    # True Negatives (TN): Not explicitly counted in base pair evaluations
    TN = len(ref_structure) * (len(ref_structure) - 1) / 2 - len(ref_pairs) - sum(1 for pair in pred_pairs if pair not in ref_pairs and not is_compatible(pair, ref_pairs))


    # Compute MCC
    if FN+FP==0 :
        return 1
    if TP+TN==0:
        return -1
    if TP+FP==0:
        return 0
    if TP+FN == 0:
        return "encore not folded"
    if TN+FP == 0:
        return "trop folded"
    if TN+FN==0:
        return "vraiment vraiment trop la"
    
    return (TP*TN-FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))




def is_compatible(pair, ref_pairs):
    """Checks if a predicted pair can be added to the reference without contradiction."""
    i, j = pair
    for (ri, rj) in ref_pairs:
        if (i < ri < j < rj) or (ri < i < rj < j):  # Crossing pairs = contradiction
            return False
    return True  # If it doesn't contradict, we ignore it





#######end Madi's modif#####################################

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

import numpy as np

def evaluate_structure(reference, predicted):
    """
    Compute Sensitivity, Specificity, PPV, and F-score
    between a reference and predicted RNA secondary structure.
    """

    ref_pairs = dot_bracket_to_pairs(reference)
    pred_pairs = dot_bracket_to_pairs(predicted)

    TP = len(ref_pairs & pred_pairs)  # True Positives
    FP = sum(1 for pair in pred_pairs if pair not in ref_pairs and not is_compatible(pair, ref_pairs))
    FN = len(ref_pairs - pred_pairs)  # False Negatives
    TN = len(ref_pairs) * (len(ref_pairs) - 1) / 2 - len(ref_pairs) - sum(1 for pair in pred_pairs if pair not in ref_pairs and not is_compatible(pair, ref_pairs))

    print(TP,FP,FN, TN)
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TP + FP) > 0 else 0
    PPV = TP / (TP + FP) if (TP + FP) > 0 else 0
    F_score = 2 * (PPV * sensitivity) / (PPV + sensitivity) if (PPV + sensitivity) > 0 else 0

    return sensitivity, specificity, PPV, F_score


