'''
Created on Oct 28, 2014

Contact: pachlioptas@gmail.com

Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import numpy as np
from time import time
from graph_analysis import IO
import scipy.linalg as linalg


def load_correlation_matrix(file_path):
    return np.loadtxt(file_path)
    

def decouple_correlation_matrix(corr_matrix, flip_negatives=True):
    pos_correlations = np.copy(corr_matrix)        
    pos_correlations[corr_matrix < 0] = 0
        
    neg_correlations = np.copy(corr_matrix)
    neg_correlations[corr_matrix > 0] = 0
    if flip_negatives:
        neg_correlations = np.abs(neg_correlations)
#    assert (np.all(pos_correlations - neg_correlations == corr_matrix))

    return pos_correlations, neg_correlations
        

def correlation_matrix_to_laplacian(corr_matrix):
    '''
    Converts a correlation matrix to a laplacian associated to the graph reflecting the correlations as edges weights.
    Input:
        corr_matrix:    mXm numpy.array corresponding to a dense, square symmetric matrix capturing
                        the correlations of random variables i and j at its (i,j) entry.
    Output:
        laplacian:      mXm numpy.array corresponding to the laplacian of a graph that links nodes i,j via an edge 
                        with weight equal to corr_matrix(i,j)     
    '''     
    laplacian = np.copy(corr_matrix)
    np.fill_diagonal(laplacian, 0)    
    degrees   = np.sum(laplacian, axis=1)
    laplacian = -laplacian    
    np.fill_diagonal(laplacian, degrees)
    return laplacian


def verify_correlation_matrix(corr_matrix):
    m, n = corr_matrix.shape
    # The correlation matrix is a square symmetric matrix with 1s in the diagonal.
    d    = corr_matrix.diagonal()
    if sum(d) != m:
        return False    
    if np.all(corr_matrix == corr_matrix.T) == False:
        return False
    #TODO assert entries in [-1, 1]
    return True

     
if __name__ == '__main__':
    corr_matrix_path = '/Users/optas/Documents/Financial_Networks/Corr_M_2013.txt'
    corr_matrix = load_correlation_matrix(corr_matrix_path)
    assert(verify_correlation_matrix(corr_matrix))
    
    pos_corr, neg_corr = decouple_correlation_matrix(corr_matrix)
    pos_laplacian      = correlation_matrix_to_laplacian(pos_corr)    
    neg_laplacian      = correlation_matrix_to_laplacian(neg_corr)
        
    evals, evecs = linalg.eigh(pos_laplacian, eigvals=(0,200))
    evecs        = np.transpose(evecs)                       
    graph_analysis.IO.save_data('pos_200_evals_evecs', evals, evecs)
    
    evals, evecs = linalg.eigh(neg_laplacian, eigvals=(0,200))
    evecs        = np.transpose(evecs)                       
    graph_analysis.IO.save_data('_200_evals_evecs', evals, evecs)
    
    
