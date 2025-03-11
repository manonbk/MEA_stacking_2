""" This module completes a nussinov algorithm using a cost function"""


def nussinov_matrix(list_of_structures, cost_function, init_cost, m=2):
    """
    Constructs a nussinov Matric, minimizing the cost function
    Args : list_of_structures : list of the aligned structures
        const_function function -> tuple of int (indices), list of structures -> int returns the cost of adding a base pair ij with regard to the list of aligned structures
        init_cost : function -> int (indice), list of structures -> int computes the prize of being unpaired to initialize the matrix
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    n=len(list_of_structures[0])
    nussinovMatrix = [[0 for i in range (n)]for j in range(n)]

    #initialize for all i-j<=m
    for k in range(m+1):
        i=0
        j=i+k
        while (i<n and j<n):
            if j==i:
                nussinovMatrix[i][j]=init_cost(i, list_of_structures)
            else :
                nussinovMatrix[i][j]=nussinovMatrix[i][j-1]+init_cost(j, list_of_structures)
            i+=1
            j+=1

    #complete matrix
    for k in range(m, n):
        i=0
        j=i+k
        while (i<n and j<n):
            if (i+1<n and j-1>=0):
                alpha0= nussinovMatrix[i+1][j-1] + cost_function((i,j), list_of_structures)
            else :
                alpha0=0
            alpha1 = max(nussinovMatrix[i][k-1] + nussinovMatrix[k][j]
                         for k in range(i, j+1))
            nussinovMatrix[i][j]=max(alpha0, alpha1)
            i+=1
            j+=1
    return nussinovMatrix



def nussinov_traceback_rec(list_of_structures,nussMat, cost_function, indices,m):
    """
    The recursive function for the Nussinov Traceback
    Args : list_of_structures : list of the aligned structures
        const_function function -> tuple of indices, list of structures -> int returns the cost of adding a base pair ij with regard to the list of aligned structures
        init_cost : a cost fonction to initialize the matrix
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    previous_bases=[]
    i,j=indices
    found=False
    if j-i<=m:
        return []
    else :
        n=len(nussMat[i])
        cost=nussMat[i][j]

        for k in range(i+1,j+1):
            if cost== (nussMat[i][k-1] + nussMat[k][j]) and not found:
                found=True
                previous_bases=nussinov_traceback_rec(list_of_structures,nussMat, cost_function,(i,k-1),m) + nussinov_traceback_rec(list_of_structures,nussMat, cost_function,(k,j),m)

        if (i+1<n and j-1>=0) and not found:
            alpha0= nussMat[i+1][j-1] + cost_function(indices, list_of_structures)
            if cost==alpha0:
                previous_bases=nussinov_traceback_rec(list_of_structures,nussMat, cost_function,(i+1,j-1),m)
                previous_bases.append((i,j))
                found=True
                #print("happens")
        
        if not found:
            print("not found", i, j)
        
        return previous_bases
    

def nussinov_traceback(list_of_structures,nussMat, costFun,m):
    """"
    The traceback algorithm for the nussinov matrix
    Args : list_of_structures : list of the aligned structures
        const_function function -> tuple of indices, list of structures -> int returns the cost of adding a base pair ij with regard to the list of aligned structures
        init_cost : a cost fonction to initialize the matrix
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    n=len(nussMat)
    return nussinov_traceback_rec(list_of_structures,nussMat, costFun, (0,n-1),m)