'''
Created on Mar 24, 2014

Contact: pachlioptas@gmail.com

Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''
 

import math, sys, itertools, os
from time import time
from scipy import spatial
from scipy.sparse import linalg
from heapq import nsmallest
from graph_analysis import IO
import numpy as np
from math import factorial
from graph_tool import spectral
from basic_utils import is_increasing


def heatKernel(tp, value):
    return math.exp(-1 * tp * value)    

def waveKernel(tp, value, sigma=1):    
    return math.exp(-1 *  ((tp - math.log(value))**2 ) / (2 * (sigma**2)))

def diffusionKernel(tp, value, sigma=1):    
    return value**tp                        


def LFS(eigvalues, eigvectors, timeSamples, kernel, args=None):
    '''
    Computing the Laplacian based signatures of the nodes of a graph for which the eigen-pairs of its Laplacian matrix
    are provided with the -eigvalues- and -eigvectors-. 
    kernel = the function which will be used to compute the signatures.
    args   = any extra arguments the kernel function needs.
    This function "Easy to read"  and thus good for debugging, but it is too slow for applications.
    '''   
    assert (len(eigvalues) == eigvectors.shape[0])          #eigvectors = eigvectors x nodes
    n = eigvectors.shape[1]
    print "Computing LFS with %d eigenvalues for %d nodes" % (len(eigvalues), n)    
    signatures = np.zeros((n, len(timeSamples)))                 

    def calculateNorm(eigvalues, tp, kernel, args=None):
        norm = 0
        if args != None:
            for v in eigvalues:
                norm += kernel(tp, v, args)
            norm = float(1)/norm
        else:
            for v in eigvalues:
                norm += kernel(tp, v)
            norm = float(1)/norm
        return norm
    
    if args != None:
        for t, tp in enumerate(timeSamples):
            norm = calculateNorm(eigvalues, tp, kernel, args)
            for node in xrange(n):
                for vector, value in zip(eigvectors, eigvalues):
                    signatures[node, t] += norm * kernel(tp, value, args) * ((vector[node])**2)
    else:
        for t, tp in enumerate(timeSamples):
            norm = calculateNorm(eigvalues, tp, kernel)
            for node in xrange(n):                            
                for vector, value in zip(eigvectors, eigvalues):
                    signatures[node, t] += norm * kernel(tp, value) * ((vector[node])**2)
    return signatures


def HKS(eigvalues, eigvectors, time):
    assert (len(eigvalues) == eigvectors.shape[0])           # eigvectors.shape = #eigvectors x #nodes.
    print "Computing Heat Kernel Signature."
    n = eigvectors.shape[1]  # Number of nodes.
    e = math.exp(1)
    signatures = np.empty((n, len(time)))

    squaredEigVectors = eigvectors**2
    squaredEigVectors = np.transpose(squaredEigVectors)
        
    for t, tp in enumerate(time):
        interm = e**(-tp*eigvalues) 
        normFactor = 1 / np.sum(interm)     
        for i in xrange (n):
            signatures[i, t] = np.dot(interm, squaredEigVectors[i]) * normFactor
      
    return signatures


def WKS(eigvalues, eigvectors, energies, sigma=1, verbose=True):    
    assert (len(eigvalues) == eigvectors.shape[0])  #dim(eigvectors) = #eigvectors x dimension of space=#nodes        
    if verbose: print "Computing Wave Kernel Signatures with %d eigenpairs." % (len(eigvalues, )),
    n = eigvectors.shape[1]  #number of nodes
    e = math.exp(1)
    log = np.log
    signatures        = np.empty((n, len(energies)))
    
    start = time()        
    squaredEigVectors = eigvectors**2
    squaredEigVectors = np.transpose(squaredEigVectors)    

    sigma = 2 * (sigma**2)
    for t, tp in enumerate(energies):
        interm = e ** (-1 * ( ( (tp - log(eigvalues))**2) / sigma))
        normFactor = 1 /np.sum(interm)
        for i in xrange (n):
            signatures[i, t] = np.dot(interm, squaredEigVectors[i]) * normFactor

    assert(np.alltrue(signatures >= 0))
    if verbose: print "Time took =  %0.3f secs." % (time() - start)
    return signatures


def cycle_signature(inGraph, maxCycleLength, mode="transition", normalize=True):    
    assert(mode == "transition" or mode == "adjacency")
    if mode == "transition":
        A    = spectral.transition(inGraph)
    if mode == "adjacency":
        A    = spectral.adjacency(inGraph)
    
    temp = A.copy()
    n    = A.shape[0]    # Number of nodes.
    sigs = np.empty((maxCycleLength, n))

    if mode == "transition":                
        for t in xrange(maxCycleLength):
            temp = temp.dot(A)            
            sigs[t,:] = temp.diagonal()

        
    if mode == "adjacency":
        for t in xrange(maxCycleLength):
            temp = temp.dot(A)
            if normalize:                
                sigs[t,:] = temp.diagonal() / np.sum(temp.diagonal())  # factorial(t+1)
            else:
                sigs[t,:] = temp.diagonal()

    return np.transpose(sigs)


def simrank(G, r=0.9, max_iter=10, eps=1e-10): 
    '''
    Easy to read but slow implementation.
    '''    
    nodes = [v for v in G.vertices()]
    nodes_i = {k: v for(k, v) in [(nodes[i], i) for i in range(0, len(nodes))]}
    sim_prev = np.zeros(len(nodes))
    sim = np.identity(len(nodes))
    
    for i in xrange(max_iter):        
        if np.allclose(sim, sim_prev, atol=eps): break
        sim_prev = np.copy(sim)
        for u, v in itertools.product(nodes, nodes):
            if u is v: continue
            s_uv = sum([sim_prev[nodes_i[u_n]][nodes_i[v_n]] for u_n, v_n in itertools.product(u.out_neighbours(), v.out_neighbours())])
            sim[nodes_i[u]][nodes_i[v]] = (r * s_uv) / (u.out_degree() * v.out_degree())
    return sim

    
def logTime(tmin, tmax, nstep):
    '''
    Given a range [tmin, tmax) and the number of samples -nstep- sample logarithmically.    
    '''    
    stepsize = (math.log(tmax) - math.log(tmin)) / nstep;
    logts = [math.log(tmin)]    
    for k in xrange(nstep):
        logts.append(logts[-1]+stepsize)     
    return np.array([math.exp(i) for i in logts])    
    
def wks_energy_generator(minEval, maxEval, timePoints, shrink=1):
    emin = math.log(minEval)   
    if shrink != 1:
        emax = math.log(maxEval) / float(shrink)
    else:
        emax = math.log(maxEval)
    
    emin = abs(emin)
    emax = abs(emax)
    
    if emax <= emin:
        print "Warning: too much shrink. - Will be set manually."
        emax = emin + 0.05*emin
    
    delta = ( emax - emin ) / timePoints    
    sigma = 10*delta
    
    res = [emin]
    for t in xrange(1, timePoints):
        res.append(res[-1] + delta)
    assert(is_increasing(res))
    return res, sigma        

def hks_time_sample_generator(minEval, maxEval, timePoints, shrink=1):    
    assert (minEval <= maxEval)    
    tmin = (1 * math.log(10)) / maxEval;              # Bounds based on largest-smallest eigenvalues used.
    tmax = (1 * math.log(10)) / (minEval*shrink);     
    
    if tmax <= tmin:
        print "Warning: too much shrinking. - Will be set manually."
        tmax = tmin + 0.30*tmin

    assert(tmax >= tmin and tmin > 0 )
    
    stepsize = (math.log(tmax) - math.log(tmin)) / timePoints;
    logts    = [math.log(tmin)]
    
    for k in xrange(timePoints):
        logts.append(logts[-1]+stepsize)    
     
    return [math.exp(i) for i in logts] 

# def hks_time_sample_generator(minEval, maxEval, timePoints, shrink=1):
#     # Idea. Doesn't seem fruitful.
#     return range(1, timePoints)


def laplacian_pseudo_inverse(inGraph=None, evals=None, evecs=None):
    
    if evals!=None and evecs !=None:
        nodes = evecs.shape[1]
        Lpseudo = np.zeros((nodes, nodes))
        for i in xrange(len(evals)):
            Lpseudo += (1.0/evals[i]) * evecs[i,:].dot(evecs[i,:])        
    else:
        L = spectral.laplacian(inGraph)
        Lpseudo = np.linalg.pinv(L.todense())
        
    return Lpseudo
        
        
def computeSpectrum(A, eigsNum, small=True, verbose=True):
    '''
    For input (sparse) matrix A, find its -eigsToUse- eigenvalues.
    If small=True => find the eigenvalues closest to zero, else find those with largest magnitude.
    '''
    start = time()
    if small==True:        
        evals, evecs = linalg.eigsh(A, eigsNum, sigma=0)
    else:
        evals, evecs = linalg.eigsh(A, eigsNum)        
    elapsed = time()                   
    if verbose:
        print "Time Spent for finding the top %d eigenvectors = %0.3f s" %(eigsNum, (elapsed - start))     
    evecs = np.transpose(evecs)
    return evals, evecs
  
  
def laplacian_spectrum(inGraph, graphName, eigsToCompute="all", edgeW=None, tryLoad=True, save=False, saveAt=None):
    '''
    Assumes the inGraph is connected.
    '''
    
    if eigsToCompute == "all":
        eigsToCompute = inGraph.num_vertices() - 1             

    eigsWhere    = "../Data/Spectrums/SM/"+graphName+"_"+str(eigsToCompute)+"_Eigs"
    
    if tryLoad:
        if os.path.isfile(eigsWhere):
            print "Loading Spectrum..."
            evals, evecs = IO.load_data(eigsWhere)            
            print "Spectrum Loaded."
            return evals, evecs             #If the data are loaded then the function exist without re-saving the loaded data
            
    if edgeW == None: #construct the Laplacian matrix
        print "Using simple Laplacian."
        Lg = spectral.laplacian(inGraph)
    else:
        print "Using weighted Laplacian."
        Lg = spectral.laplacian(inGraph, weight=edgeW)
    
    if eigsToCompute == inGraph.num_vertices() - 1: # We will compute the largest eigenvalues for efficiency.
        evals, evecs = computeSpectrum(Lg, eigsToCompute, small=False)
        assert(is_increasing(evals))
    else:
        evals, evecs = computeSpectrum(Lg, eigsToCompute+1, small=True)
        #Discard The 1st eigenpair (0-constant)
        evals = evals[1:]; evecs = evecs[1:, :]
    
    if save:
        if saveAt != None:
            IO.save_data(saveAt, evals, evecs)
        else:
            IO.save_data(eigsWhere, evals, evecs)
        
    return evals, evecs


def filterEigenPair(evals, evecs, keep, magnitude="small", balance=None):
    '''
    Keep -keep- many eigenvalues and corresponding eigenvectors from the input matrices evals and evecs.
    If magnitude == small => keep the smallest eigenvalues (starting from the one at evals[0])
    If magnitude == large => keep the largest eigenvalues (starting from the one at evals[0])
    '''
    assert ( (len(evals) == evecs.shape[0]) and keep <= len(evals) and is_increasing(evals) )
    
    if magnitude == 'small':
        return evals[:keep], evecs[:keep, :]
    elif magnitude == 'big':
        return evals[-keep:], evecs[-keep:, :]
    
    elif magnitude == 'both':
        if balance == None:         
            small = int(np.ceil(keep / 2.0))
            big   = keep - small
        else:
            small = balance[0]
            big   = balance[1]
            
        vals = np.hstack( (evals[:small], evals[-big:]) )
        vecs = np.vstack( (evecs[:small, :], evecs[-big:, :]) )
        return vals, vecs   
    else:
        assert(False)



def neighbor_signatures(inSigs, graph, aggregate):
    assert(aggregate == "average" or aggregate == "sum" or aggregate == "max" or aggregate == "min")
    n, t = inSigs.shape
    assert(graph.num_vertices() == n) 

    res = np.zeros((n,t))       
    if aggregate == "average":
        print "Aggregating average neighbor signatures."
        for v in graph.vertices():
            for nv in v.out_neighbours():
                res[int(v)] += inSigs[int(nv)]
            res[int(v)] = res[int(v)] / float(v.out_degree())            
    elif aggregate == "sum":
        for v in graph.vertices():
            for nv in v.out_neighbours():
                res[int(v)] += inSigs[int(nv)]
    elif aggregate == "min":
        for v in graph.vertices():
            for nv in v.out_neighbours():
                res[int(v)] = np.minimum(res[int(v)], inSigs[int(nv)])                
            res[int(v)] = res[int(v)] / float(v.out_degree())
            
    else: #max
        for v in graph.vertices():
            for nv in v.out_neighbours():
                res[int(v)] = np.maximum(res[int(v)], inSigs[int(nv)])
            res[int(v)] = res[int(v)] / float(v.out_degree())
    return res
    
    
    
        


#########################  CODE THAT NEEDS CLEANING #########################
def signLessLaplacianSpectrum(inGraph, graphName, eigsToCompute="all", save=True):
    '''
    Assumes the inGraph is connected.
    '''
    if eigsToCompute == "all":
        eigsToCompute = inGraph.num_vertices() - 1             
    eigsWhere    = "../Data/Spectrums/SignLess_L/"+graphName+"_"+str(eigsToCompute)+"_Eigs"
    if os.path.isfile(eigsWhere):
        print "Loading Spectrum"
        evals, evecs = IO.load_data(eigsWhere)
    else:
        Ld  = spectral.laplacian(inGraph).diagonal()       #TODO replace with degree vector to avoid building the laplacian 
        Ls  = spectral.adjacency(inGraph)
        Ls.setdiag(Ld)        
        evals, evecs = computeSpectrum(Ls, eigsToCompute, small=False)  #For sign less Laplacian big eigpairs matter
        #sort so that big eigen-pairs are first
        
        evals = evals[::-1]
        evecs = np.flipud(evecs)        
        if save:
            IO.save_data(eigsWhere, evals, evecs)
    return evals, evecs

def countCorrespondencesManyLevels(distanceMatrix, forwardMap, topKs, fast=True):
    '''
    topKs: How many mismatches we will search over before we give up.  
    '''
    total = distanceMatrix.shape[0] #number of points        
    topMatches = list()
    maxTopK = max(topKs)
    
    #for each data point find its K closest neighbors (quadratic)
    if fast:
        for i, pt in enumerate(distanceMatrix):        
            topMatches.append(nsmallest(maxTopK, range(total), key=distanceMatrix[i,:].__getitem__))
    else:
        for i, pt in enumerate(distanceMatrix):
            topMatches.append([j[0] for j in sorted(enumerate(pt), key=lambda x:x[1])][:maxTopK])            

    res = []
    for topK in topKs:
        matches = 0    
        for i in xrange(total):    
            if forwardMap[i] in topMatches[i][:topK]:
                matches += 1
        res.append((topK, matches, total, matches/float(total)))
    return res

def measureSignaturePreservation(from1To2Map, g1Sigs, g2Sigs, distanceF):
    n = g1Sigs.shape[0] 
    assert(n == g2Sigs.shape[0] == len(from1To2Map))
    distances = 0 
    for pt in xrange(len(from1To2Map)):
        distances += distanceF(g1Sigs[pt,:], g2Sigs[from1To2Map[pt],:])
    return distances / float(n) 
                
def measureMappingRobustness(from1To2Map, nodeMapping, sigFrom, pedantic=False):
        
    if pedantic:
        perservL = []
        distroL = []
        missedL = []
        inside = 0
        for pt in xrange(len(from1To2Map)):
            ptKNN = nodeMapping[pt][1]        
            if from1To2Map[pt] == ptKNN[0]:    #mapping was the best 
                perservL.append((nodeMapping[pt][0][1] - nodeMapping[pt][0][0]) / nodeMapping[pt][0][0])
                inside += 1
            elif from1To2Map[pt] not in ptKNN:  #was completely missed                            
                missedBy = np.linalg.norm(sigFrom[pt] - sigFrom[ptKNN[0]])  / nodeMapping[pt][0][0]                
                missedL.append(missedBy)
                inside += 1                                
            else:  
                posOfHit = np.where(ptKNN==from1To2Map[pt])[0][0]             
                distroL.append((nodeMapping[pt][0][posOfHit] - nodeMapping[pt][0][0]) / nodeMapping[pt][0][0])
                inside += 1        
        print inside, len(perservL) + len(missedL) + len(distroL)
        return perservL, distroL, missedL
    
    else:
    
        perserv = perservC = distor = distorC = missedBy = missedC = 0 

        for pt in xrange(len(from1To2Map)):
            ptKNN = nodeMapping[pt][1]        
            if from1To2Map[pt] == ptKNN[0]:    #mapping was the best 
                perserv += (nodeMapping[pt][0][1] - nodeMapping[pt][0][0]) / nodeMapping[pt][0][0]
    #             print "match " +str(from1To2Map[pt])
    #             print "best " + str(perserv)
                perservC += 1
            elif from1To2Map[pt] not in ptKNN:  #was completely missed                            
                missedBy += np.linalg.norm(sigFrom[pt] - sigFrom[ptKNN[0]])  / nodeMapping[pt][0][0]                
                missedC  +=1                    
            else:  
                posOfHit = np.where(ptKNN==from1To2Map[pt])[0][0]             
    #             print "hit at "+ str(posOfHit)
                distor += (nodeMapping[pt][0][posOfHit] - nodeMapping[pt][0][0]) / nodeMapping[pt][0][0]
    #             print "dist "+ str(distor)
                distorC += 1
        
        
        return (perservC and perserv/perservC), (distorC and distor/distorC), (missedC and missedBy/missedC), perservC, distorC, missedC 
                  
def computeISOSignatures(laplacian1, laplacian2, kernel, eigsToUse, timeSamples):
    '''Input:   -laplacian1-, -laplacian2-  are the Laplacians corresponding two 2 graphs
                -kernel- the type of kernel function to be used for the signature construction
                -eigsToUse- number of eigenvectos to be used
                -timeSamples-  time points where the kernel will be  evaluated
      Output:
               sig1, sig2: The signatures of all the nodes corresponding two the two graphs
          
    '''
    
    start = time()    
    evals, evecs = linalg.eigsh(laplacian1, k=eigsToUse, sigma=0)

    elapsed = time()       
    if "DEBUG" == True:
        print "Time Spent for finding top %d eigenvectors of 1st Laplacian= %0.3f s" %(eigsToUse, (elapsed - start))
        
    evecs = np.transpose(evecs)
    
    start = time()              
    g1Sigs = LFS(evals, evecs, timeSamples, kernel)
    elapsed = time()    
    if "DEBUG" == True:
        print "Time Spent for Making Signature with %d time samples for 1st graph = %0.3f s" %(len(timeSamples), elapsed - start)
    
    if laplacian2 ==None:
        return g1Sigs, None
    
    evals, evecs = linalg.eigsh(laplacian2, k=eigsToUse, sigma=0)       
    evecs = np.transpose(evecs)
    g2Sigs = LFS(evals, evecs, timeSamples, kernel)
    
    return g1Sigs, g2Sigs
        
def signatureBasedMapping(sig1, sig2, kernel, kneighbors):
    t0 = time()
    res = []
    if kernel == heatKernel:
        KDTree=spatial.cKDTree(sig2, leafsize=100) #TODO: Play with the leafsize to get the fastest result for your dataset
        n = sig1.shape[0]    
        for i in xrange(n):
            res.append(KDTree.query(sig1[i,:], k=kneighbors))                    
    else:
        assert(False)
    
    if "DEBUG":
        print "Time Spent for Calculating the Mapping between the two graphs = " + str(time()-t0)
            
    return res
          
def BF_SignatureBasedMapping(sig1, sig2, kernel, kneighbors):
    if kernel == heatKernel:
        from sklearn.metrics.pairwise import euclidean_distances
        distanceMatrix = euclidean_distances(sig1, sig2)
        topMatches = []        
        for i, pt in enumerate(distanceMatrix):        
            topMatches.append(nsmallest(kneighbors, range(distanceMatrix.shape[0]), key=distanceMatrix[i,:].__getitem__))
        return topMatches

def mappingAccuracy(mappings, groundTruth, topKs, equivsG=None, verbose=True):
    res = []
    n = len(mappings)

    for topK in topKs:
        robustNode  = 0
        mappedToEC  = 0
        matches     = 0    
        for i in xrange(n):    
            if groundTruth[i] in mappings[i][1][:topK]:
                matches += 1
                if equivsG != None and groundTruth[i] in equivsG.keys():
                    robustNode += 1 #This node was correctly matched despite the fact that it is in an non trivial Equivalence Class.
                continue
            if equivsG != None:
                try:
                    for eq in equivsG[i]:
                        if groundTruth[eq] in mappings[i][1][:topK]:
                            mappedToEC += 1
                            matches    += 1
                            break
                except: #node -i- does not have any equivalent nodes
                    continue
        if verbose and equivsG!=None:
            print "Out of the %d Matched Nodes, %d were mapped in a node in their Equivalence class and %d despite the fact the were inside a non-trivial \n"\
            "Equiv. Class, were mapped to the exact isomorphic node directly. (TopK=%d). (Nodes in some EC=%d)" %(matches, mappedToEC, robustNode, topK, len(equivsG.keys()))
            
        res.append((topK, matches, n, matches/float(n)))
            
    return res


def isLaplacian(L, verbose =True):        
    n = L.shape[0]
    if verbose:  
        e = sum(L.diagonal()) / 2
        print "Inside: \'isLaplacian\': The graph has %d node and  %d edges" % (n, e)

    for i in xrange(n):
        row = L.data[L.indptr[i] : L.indptr[i + 1]] #Corresponds to the ith-row of the matrix, only the non-zero entries        
        diagFound = False        
        row_sum   = 0
        for item in row:
            if (item != -1 and diagFound): return False
            if item != -1:
                diagFound = True                
            row_sum  += item            
        if row_sum != 0: return False
    return True

    

    
def waveTimeSample(minEval, maxEval, timePoints, sigma=None):
    assert (minEval < maxEval)
    
    if maxEval-minEval <= 0.01:
        print "Minimum eigenvalue very close to maximum one. We will manually adjust."
        maxEval = maxEval + ((5/float(100))*maxEval)

    if sigma == None:
        sigma = (7 / (timePoints + float(29)) ) * ( math.log(maxEval) - math.log(minEval))
    
    assert (sigma >= 0)
    
    emin = math.log(minEval) + (2 * sigma) #minimum energy
    emax = math.log(maxEval) - (2 * sigma)
    
    inc  = (emax - emin) / float(timePoints)
    assert (inc >= 0)

    res = [emin]
    for t in xrange(1, timePoints):
        res.append(res[-1] + inc)
        
    assert (is_increasing(res))
    
    #Mitigate negative values
    if res[0] < 0:
        res = [abs(t) for t in res]
    res.sort()
    assert(is_increasing(res))            
    return res, sigma
        

def waveTimeSample2(minEval, maxEval, timePoints):
    '''
    Implementing the strategy for energy sampling as was proposed in 
    
    '''
    assert (minEval < maxEval)
    #Strategy used in their 2nd paper: 
    emin = math.log(minEval)
    emax = math.log(maxEval) / 1.02
    
    delta = ( emax - emin ) / timePoints
    sigma = 7*delta
    
    res = [emin]
    for t in xrange(1, timePoints):
        res.append(res[-1] + delta)
    assert(is_increasing(res))
    return res, sigma


#########################  CODE THAT NEEDS CLEANING </end> #########################

if __name__ == "__main__":
    print "Running main() of -graphKernels- module."
    graphName = 'Zarathustra_English'
    eigsToCompute = "all"
    print "Will compute %s eigen-pairs." % eigsToCompute,
    inGraph = IO.load_data('../Data/Graphs/Novels/'+graphName+".GT.graph").next()
    weights = inGraph.ep['co-occurrence']    
    evals, evecs = laplacian_spectrum(inGraph, graphName, eigsToCompute=eigsToCompute, edgeW=weights, tryLoad=False, save=False)
    eigsWhere    = "../Data/Spectrums/SM/"+graphName+"_weighted_all_Eigs"
    IO.save_data(eigsWhere, evals, evecs)    
    sys.exit("Exiting Successfully.")
    