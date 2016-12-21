'''
Created on Feb 1, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import numpy as np
from scipy.stats.mstats import rankdata

def jaccard_vectorial(v1, v2):
    '''
    Computes the Jaccard similarity between two real vectors. That is, sum_i{ min(v1[i], v2[i])}  / sum_i{ max(v1[i], v2[i])}

    Parameters
    ----------
    v1 : array-like
        Any array of arbitrary size.
    v2 : array-like
        Any other array of identical size.

    Returns
    -------
    jaccard : float
        Jaccard metric returned is a float on range [0,1].
        Maximum similarity = 1
        No similarity = 0
    
    Notes
    -----
    The order of inputs for `jaccard` is irrelevant. The result will be
    identical if `v1` and `v2` are switched.
    '''
    
    if type(v1) == np.ndarray and v1.shape != v2.shape:
        raise ValueError("Dimension mismatch: v1 and v2 must have the same dimensions.")     
    joined = np.vstack((v1, v2))
    return np.sum(joined.min(axis=0)) / np.sum(joined.max(axis=0))


def pdist_to_p_rank(pdist_matrix):
    ''' Converts a square distance matrix to a ranking matrix.
    That is, for each row of the input matrix it converts each distance of it, to an integer reflecting 
    the order/ranking of the distance. The smallest distance is encoded with zero. ties are broken by 'average'.    
    '''
    return rankdata(pdist_matrix, axis = 1)


def exclude_black_listed(distanceM, blackList, trueClusters):    
    
    '''
    Retains from the input square distance matrix, only the columns and rows listed in the whiteList.
    '''
    
    whiteList = set(np.arange(len(trueClusters)))    
    for i in xrange(len(blackList)):
        whiteList.difference_update(set(np.where(trueClusters == blackList[i])[0]))

    whiteList = list(whiteList)
    whiteList.sort()    

    return distanceM[whiteList, :][:, whiteList], trueClusters[whiteList]
    
    
    
    
    