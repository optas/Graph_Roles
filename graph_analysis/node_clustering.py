'''
Created on Feb 21, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import IO
import numpy as np
from collections import defaultdict
from scipy.misc import comb

def count_edges_between_clusters(g, clusters, edgeWeights=None, blackList=None, verbose=True):
    '''
    Counts the number of edges among different node clusters. 
     
    Input:
            g               -    (graph_tool) A graph.
            clusters         -   (np.array) identifying the cluster that each node of g belongs.
            blackList       -    (list or ndarray) The IDs of nodes that should be excluded from the output statistics.    
            edgeWeights     -     Edge property map, If provided the weight of edges will be reported separately.
            verbose         - 
    TODO:                    Enable blackListing and edgeWeights and directionality.
    '''
    

    edgesCount = defaultdict(lambda : 0)
    
    if edgeWeights != None:    
        weightOfCategory  = defaultdict(lambda : 0)
    
    for node1, node2 in g.edges():        
        key = tuple(sorted([clusters[node1], clusters[node2]]))  #Change if you want to count edge weight
        edgesCount[key] += 1        
            
    if verbose:    
        for categ1, categ2 in edgesCount.keys():
            if categ1 != categ2:
                edgeDensity = edgesCount[(categ1, categ2)] / float(np.sum(clusters == categ1) * np.sum(clusters == categ2))
            else:
                edgeDensity = edgesCount[(categ1, categ2)] / comb(np.sum(clusters == categ1), 2)

            print categ1,  categ2, edgeDensity
    

def charachteristics_of_node_clusters(g, clusters, characteristicProp="clusterSize", verbose=True):
    '''    Input:
                g               -    (graph_tool) A graph.
                clusters         -    Vertex property map, identifying the class of each node of g. Alternatively, a np.array carrying the same information.
    '''
    n = float(g.num_vertices())
                            
    if characteristicProp == "clusterSize":
        if verbose:
            print "Sizes of node clusters: \n<Class ID> <Size> <Relative Size>"
        clusterSizes = []        
        for i in xrange(min(clusters), max(clusters)+1):
            cluster_i_size = np.sum(clusters==i)
            clusterSizes.append((i, cluster_i_size))
            if verbose:
                print i, cluster_i_size, cluster_i_size/n
                                            
        return clusterSizes
    
    
def average_signature_of_clusters(signatures, clusters, blackList=None):    
    '''
        signatures     -    num_nodes x features
    '''
    assert(len(clusters) == signatures.shape[0])
    res = {}
    for i in xrange(min(clusters), max(clusters)+1):
        if blackList != None and i in blackList: continue
        clusterMembers = np.where(clusters == i)[0]
        clusterMembers = [[member] for member in clusterMembers] # Represent each member as a list [] for indexing the numpy.array.
        res[i] = np.average(signatures[clusterMembers,:], axis=0)    
    return res