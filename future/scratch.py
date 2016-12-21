'''
Keeping snipets of code that might be of some use or not.

Created on Feb 22, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

from graph_tool import Graph, generation
import numpy as np

def gnp(n, e):        
    '''
    Return an Erdos-renyi graph with n vertices and e edges
    '''
    g = Graph(directed=False)
    g.add_vertex(n)
   
    edgesAdded  = 0
    edgesNeeded = True
    while edgesNeeded:
        fromNode = range(n); toNode = range(n);
        np.random.shuffle(fromNode)
        np.random.shuffle(toNode) 
        for i in xrange(n):
            if fromNode[i] == toNode[i] or g.edge(fromNode[i], toNode[i]) != None:  continue #No self-loops or parallel edges
            g.add_edge(fromNode[i], toNode[i])
            edgesAdded += 1
            if  edgesAdded >= e: 
                edgesNeeded = False
                break
    generation.random_rewire(g, model="erdos")
    return g


def example_isClosestNeighborMatch():
    '''    
    Quick and Dirty: Is nearest neighbor of a node the same kind as she (i.e., part of a clique)     
    '''
    
    gFinal= erdos_circles_cliques([1000, 5000], nCliques=20, cliqueSize=15, nCircles=20, circleSize=30, save=False)
        
    cliques = set(np.where(gFinal.vp["cliq"].a == 1)[0])    
    circles = set(np.where(gFinal.vp["circ"].a == 1)[0])
    
    evals, evecs =  gk.laplacian_spectrum(gFinal, graphName=" ", eigsToCompute="all")

    eigsToUse = [5, 10, 20, 50, 100]
    res = {"cliq":[], "circ":[]}
    for k in eigsToUse:
        timeSample1 = gk.waveTimeSample2(minEval, maxEval, timePoints)
        sig1 = gk.WKS(evals, evecs, k, timeSample1)
 
        KDTree=spatial.cKDTree(sig1, leafsize=100)
        n = sig1.shape[0]
        print n
         
        clicF=circF=0
        for i in xrange(n):        
            if i in cliques and KDTree.query(sig1[i,:], k=2)[1][1] in cliques: clicF +=1                
            if i in circles and KDTree.query(sig1[i,:], k=2)[1][1] in circles: circF +=1
         
        clicAcu = clicF/float(len(cliques))
        circAcu = circF/float(len(circles))
        res["cliq"].append(clicAcu)
        res["circ"].append(circAcu)
        print k, clicAcu, circAcu
         
    label  = ["cliques", "circles"]
    inText = "Heat Kernel\nLogTimeSample |t|="+str(10)
    xlabel = "Largest Eigenvalues"
    myPLT.xQuantityVsAccuracy("syntheticMix", eigsToUse, [res["cliq"], res["circ"]], xlabel, label, inText, "-", True)
    


import numpy as np
from rpy2.robjects import r
import pandas.rpy.common as com
from pandas import DataFrame

from rpy2.robjects.packages import importr
utils = importr("utils")
importr('clusterCrit')
r('getCriteriaNames(TRUE)')


def internal_cluster_evaluation(feature_points, cluster_labels, metrics):
    fp = DataFrame(feature_points)
    fp = com.convert_to_r_dataframe(fp)
    r.assign("feature_points", fp)
    
    cl = DataFrame(cluster_labels)
    cl = com.convert_to_r_dataframe(cl)
    r.assign("cluster_labels", cl)
    
    r('cluster_labels')
#     r('intCriteria(feature_points, cluster_labels, c("C_index","Calinski_Harabasz","Dunn"))')


def to_r_frame(data, variable_name,save_path):        
#     data = {"feature_matrix" : data, "cluster_labels", labels }
    df = DataFrame(data)
    df = com.convert_to_r_dataframe(df)
    r.assign(variable_name, df)
    r("save(%s, file='%s.gzip', compress=TRUE)" % (variable_name, save_path))






def trimGroups(groups, bound):
    above = []
    below = []
    for i in np.unique(groups):
        if len(np.where(groups == i)[0]) > bound:    
            above.append(i)
        else:
            below.append(i)
    return above, below


def weighted_L1_Dist(vec1, vec2):
    '''
    A distance between 2 vectors. It was proposed in the Wave Kernel Signature paper as 
    a means to compute the similarity between two WK signature vectors. 
    weighted_L1_Dist(v1,v2) = Sum{i}  ( |v1_i - v2_i| \ |v1_i + v2_i| )
    '''        
    return np.dot(np.absolute(vec1-vec2), (float(1) / (np.absolute(vec1 + vec2))))



def pairwiseDists(X, Y, distanceF):
    res = np.empty( (X.shape[0], Y.shape[0]) )            
    for i, ptx in enumerate(X):
        for j, pty in enumerate(Y):
            res[i,j] = distanceF(ptx, pty)    
    return res



#Next 4 functions: very problem specific
def propagateLabels(d):
    res = {}
    #group must be list
    for key, group in d.iteritems():
        fullGroup = group+[key]
        res[key] = fullGroup
        for val in group:
            assert (val not in d.keys())
            res[val] = fullGroup
    return res
        
def listGroupingToDict(l):
    '''
    Input: l = List of Lists (each nested list represents an unordered grouping (e.g. node equivalence classes)
    Output: A dictionary, where every item of each nested list is becoming a 'key' with corresponding value all the items of its corresponding grouping.
    The values are stored internally as sets.
    '''
    res = dict()
    for group in l:
        for member in group:
            res[member] = group
    return res
    
def mergeEquivClasses(currentGroups):
    '''
    Input: currentGroups, a list of lists or a list of tuples. Each -inner- group corresponds to elements that belong to one Eq. class.
    e.g.
    a = [(1, 2), (6, 94, 3), (3, 5), (2, 4), (4, 8), (9, 11)]
    b = mergeEquivClasses(a)
    b = [{1, 2, 4, 8}, {3, 5, 6, 94}, {9, 11}]
    '''    

    equivClasses = list()
    equivClasses.append(set(currentGroups[0]))
        
    for i in xrange(1, len(currentGroups)):
        group_i = set(currentGroups[i])
        merged = False
        for j in xrange(len(equivClasses)):                        
            if len(group_i.intersection(equivClasses[j])) != 0:
                equivClasses[j] = equivClasses[j].union(group_i)
                merged = True
                break    
        if not merged: #make a new equiv. class
            equivClasses.append(group_i)     
            
    return equivClasses
    
def mergeTwoDictionaries(dict1, dict2, noOverlap=True):
    res = copy.deepcopy(dict1)
    if noOverlap:
        for key in dict2.keys():
            assert(key not in dict1)
            dict1[key] = dict2[key]
    else:
        for key in dict2.keys():            
            dict1[key] = dict1[key] + dict2[key]
    return res    
        
    
    
   
def breakTies(vector):
    '''
    Recalculates the ranking of the already sorted input observations 
    and asserts that if two observations have the same numeric value, 
    then they get the same rank.

    vector = [(location, rank), (location, rank)..., ]    
    '''
#     res = copy.deepcopy(vector)    
    n    = len(vector) # Todo see again
    allRanks = np.empty(n)
    allLocs  = np.empty(n)
    
    res = []
    rank = 0     
    for i in xrange(n):
#         res[i] = (rank, vector[i][0])
        allRanks[i] = rank
        allLocs[i]  = vector[i][0]
#         res.append((rank, vector[i][0]))
        if i < n-1  and vector[i][1] < vector[i+1][1]: rank +=1 #1st condition safeguards from an out_of_range Exception    
    
#     assert(is_increasing([r[0] for r in res]))
#     assert(range(n)== sorted([r[1] for r in res])) #TODO remove from code
    return allRanks, allLocs
        


def pdist_to_pRank(pdistM, tiesBreak=False):
    '''    
    tiesBreak = True, then all the objects that have the same distance from a given object will get the same rank.
    '''    
    prankM = np.empty(pdistM.shape)
    if tiesBreak: #slow
        for i, row in enumerate(pdistM):
            ranked = [t for t in sorted(enumerate(row), key=lambda x:x[1])] #[(rank, matrix_value), (rank, matrix_value)...]            
            allRanks, allLocs = breakTies(ranked)                        
#             print allRanks, allLocs
            for rank, loc in zip(allRanks, allLocs):
                prankM[i,loc] = rank
#             print prankM[i,]
               
    else:
        for i, row in enumerate(pdistM):
            temp = [j[0] for j in sorted(enumerate(row), key=lambda x:x[1])]
            for j, s in enumerate(temp):
                prankM[i,s] = j 
                   
    return prankM

            
def nmf_prediction():
    signatures = execute_method(methodName, inGraph, graphName, distances=False, ranks=False, distFunction="default", methodParams=methodParams)
    nmf        = decomposition.NMF(n_components=numOfTrueClusters, nls_max_iter=3000)
    nmfModel   = nmf.fit_transform(signatures)                
    predictedLabels = np.argmax(nmfModel, axis=1)