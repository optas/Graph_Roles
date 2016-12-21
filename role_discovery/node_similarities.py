'''
Created on Feb 16, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import numpy as np
from role_discovery import graph_kernels as gk
from graph_tool import spectral
from itertools import product
from graph_analysis import IO
from scipy.sparse import linalg, identity, diags
from basic_utils import is_increasing


def heat_sim(inGraph, graphName, eigs="all", timePoints=100, strategy=None, shrink=1):
    evals, evecs       = eigen_pairs(inGraph, graphName, eigs, strategy)
    timeSamples        = gk.hks_time_sample_generator(evals[0], evals[-1], timePoints, shrink=shrink)
    sig                = gk.HKS(evals, evecs, timeSamples)
    return sig


def wave_sim(inGraph, graphName, eigs="all", energies=100, strategy=None):
    evals, evecs        = eigen_pairs(inGraph, graphName, eigs, strategy)
    timeSamples, sigma  = gk.wks_energy_generator(evals[0], evals[-1], energies)    
    sig                 = gk.WKS(evals, evecs, timeSamples, sigma)
    return sig


def role_sim(graphName):    
    roleSimFile = "../Data/Sigs/RoleSim/" + graphName + ".greach-rs-I1-M1-dst.txt"
    distanceMatrix = IO.readRoleSimMatrix(roleSimFile)    
    return distanceMatrix


def sim_rank(graphName):    
    simRankFile = "../Data/Sigs/RoleSim/" + graphName + ".greach-sr-I1-M1-dst.txt"
    distanceMatrix = IO.readRoleSimMatrix(simRankFile)
    return distanceMatrix     


def refex_sim(graphName):
    refexSigFile       = "../Data/Sigs/refex_sim/" + graphName
    sig                = IO.load_data(refexSigFile).next()
    return sig
    
    
def vertex_sim(inGraph, alpha=.95, dist=False):
    '''
    Implementation of the node similarity metric described in
    "Vertex Similarity in networks, Leicht, Holme, Newman."
    '''    
    n = inGraph.num_vertices()
    A = spectral.adjacency(inGraph)
    evals, evecs = linalg.eigsh(A, k=1, which='LM', tol=1E-3)    
    scale = alpha / evals[0]
    S = linalg.inv(identity(n, format='csc') - scale * A.tocsc())
    d = inGraph.degree_property_map("total").a
    S = diags(1.0/d, 0) * S * diags(1.0/d, 0)
    S[np.diag_indices_from(S)] = 1              # Each node is self-similar.
    
    if dist:
        S = 1 - S.toarray()
        
    S = 0.5 * (S + S.T)                         # Force S to be symmetric and guards to numerical errors.
    
    return S

    
def cycle_sim(inGraph, cycleLength, matrixType, keepEven=False, clearConverged=False):

    if keepEven:
        sigs = gk.cycle_signature(inGraph, cycleLength, matrixType)
    else:
        sigs = gk.cycle_signature(inGraph, cycleLength*2, matrixType)
        sigs = sigs[:, [i for i in range(sigs.shape[1]) if i%2==1]]

    from numpy.linalg import norm
    print "Norm of signature matrix = " + str(norm(sigs, 'fro'))
    
    if  clearConverged:
        return remove_not_varying_parts(sigs, 0.001)
    else:        
        return sigs
     
     
def raw_laplacian_sim(inGraph, graphName, eigs="all", energies=None, strategy=None):
    '''
        Embed each graph node in its projections of laplacian eigenvectors.
    '''
    evals, evecs = eigen_pairs(inGraph, graphName, eigs, strategy)
    sig          = evecs.T
    return sig


def commute_time_dist(inGraph):
    '''
        Embed each graph node in its projections of laplacian eigenvectors.
    '''
    Lp           = gk.laplacian_pseudo_inverse(inGraph=inGraph)
    n            = inGraph.num_vertices()
    graphVolume  = np.sum(inGraph.degree_property_map("total").a)    
    dist         = np.zeros((n,n))
    for i in xrange(n):
        for j in xrange(i+1, n):
            dist[i,j] = graphVolume * (Lp[i,i] + Lp[j,j] - (2 * Lp[i,j]))
        
    dist[np.diag_indices_from(dist)] = 0
    dist = dist + dist.T    
    return dist

         
def compress_spectral_signatures(signatures, evals, sigma=.1):
    assert(signatures.shape[1] == len(evals))
    assert(sigma <= 1)
    
    allDifs = np.abs(np.array(evals[1:]) - np.array(evals[:-1]))
    averDff = np.average(allDifs)
    compressLocations = allDifs < (averDff * sigma)
    toBeCompressed    = np.sum(compressLocations)
    
    if toBeCompressed == 0:
        print "no compression at this sigma level."
        return signatures
    
    newSignatures = np.empty((signatures.shape[0], signatures.shape[1]-toBeCompressed))

    lastLocation       = 0                # Last compressed column.
    newSignatures[:, lastLocation] = signatures[:, 0];                      
    for i, merge in enumerate(compressLocations):   # i sequentially moves along the columns of the original signature (we use i+1 to start from 1).
        if merge:
            newSignatures[:, lastLocation] += signatures[:, i+1]
        else:
            lastLocation +=1
            newSignatures[:, lastLocation]  = signatures[:, i+1]
        
    return newSignatures
    

 
def panos_sim(inGraph, graphName, eigs="all", strategy=None, compressed=True, sigma=0.1):
    evals, evecs = eigen_pairs(inGraph, graphName, eigs, strategy)
    
    degrees      = inGraph.degree_property_map("total").a
    n            = evecs.shape[1]  # Number of nodes.
    signatures   = np.empty((n, len(evals)))
    
    squaredEigVectors = evecs**2        
    squaredEigVectors = np.transpose(squaredEigVectors)
                       
    for node in xrange(n):
        signatures[node, :] = squaredEigVectors[node, :] * (degrees[node] - evals)**2
  
    if compressed:
        signatures = compress_spectral_signatures(signatures, evals, sigma=sigma)
  
    print "Signature was created."
    
    return signatures


def panos_sim_2(inGraph, graphName, eigs="all", strategy=None, compressed=True, sigma=0.1):
    evals, evecs = eigen_pairs(inGraph, graphName, eigs, strategy)
    degrees      = inGraph.degree_property_map("total").a
    n            = evecs.shape[1]  # Number of nodes.
    signatures   = np.empty((n, len(evals)))
    
    squaredEigVectors = evecs**2        
    squaredEigVectors = np.transpose(squaredEigVectors)
                       
    for node in xrange(n):
        signatures[node, :] = squaredEigVectors[node, :] * ((degrees[node] - evals)**2) * degrees[node]
  
    if compressed:
        signatures = compress_spectral_signatures(signatures, evals, sigma=sigma)
  
    print "Signature was created."    
    return signatures



def panos_sim_3(inGraph, graphName, eigs="all", strategy=None, compressed=True, sigma=0.1):
    evals, evecs = eigen_pairs(inGraph, graphName, eigs, strategy)
    degrees      = inGraph.degree_property_map("total").a
    n            = evecs.shape[1]  # Number of nodes.
    signatures   = np.empty((n, len(evals)))
    
    squaredEigVectors = np.abs(evecs)        
    squaredEigVectors = np.transpose(squaredEigVectors)
                       
    for node in xrange(n):
        signatures[node, :] = squaredEigVectors[node, :] * (np.abs(degrees[node] - evals)) * degrees[node]
  
    if compressed:
        signatures = compress_spectral_signatures(signatures, evals, sigma=sigma)
  
    print "Signature was created."    
    return signatures



def panos_sim_4(inGraph, graphName, eigs="all", strategy=None, compressed=True, sigma=0.1):
    
    A        = spectral.adjacency(inGraph)
    degrees  = inGraph.degree_property_map("total").a
    M        = A + diags(degrees+1, 0)
    
    if eigs == "all":
        eigs = len(degrees) - 1
    
    evals, evecs = gk.computeSpectrum(M, eigs, small=False, verbose=True)
    print evals
    
    n            = evecs.shape[1]  # Number of nodes.
    signatures   = np.empty((n, len(evals)))
    
    squaredEigVectors = evecs ** 2
    squaredEigVectors = np.transpose(evecs)
                       
#     for node in xrange(n):
#         signatures[node, :] = squaredEigVectors[node, :] * ((degrees[node] - evals)**2)
#         
#     for i in xrange(10):
#         signaturesNew = gk.neighbor_signatures(signatures, inGraph, "average")    
#         if np.allclose(signatures, signaturesNew):
#             break
#     signatures = signaturesNew
  
    if compressed:
        signatures = compress_spectral_signatures(signatures, evals, sigma=sigma)
  
    print "Signature was created."    
    return signatures



def neighbors_matter(inGraph, dists):
    n           = dists.shape[0]
    sim         = np.exp(- (dists**2) / (2*(np.median(dists)**2)))
    res1        = np.zeros((n, n))
    
#     res2        = np.zeros((n, n))
    for i in xrange(n):
        iNeighbs = [int(ni) for ni in inGraph.vertex(i).out_neighbours()]
        for j in xrange(i+1, n):
            jNeighbs = [int(nj) for nj in inGraph.vertex(j).out_neighbours()]            
            for ni, nj in product(iNeighbs, jNeighbs):            
                res1[ni, nj] += sim[i,j]
#                 res2[ni, nj] *= sim[i,j]
    
    res1        = res1 + res1.T
    res1        = (res1.T / np.array(np.max(res1, 1), dtype=np.float)).T
    res1        = 0.5 * (res1 + res1.T)
    
    res1        = 1 - res1
    res1[np.diag_indices_from(res1)] = 0
    
#     res2        = res2 + res2.T
#     res2        = 1 - res2
#     res2[np.diag_indices_from(res2)] = 0
#     return res1, res2
    return res1
    
    
def pick_eigs_exhaust_derivative_(evals, stop=10):    
    for i in range(1, len(evals)-1):
        previousDer = abs(evals[i] - evals[i-1])
        nextDer     = abs(evals[i+1] - evals[i])
        if (nextDer / float(previousDer)) > stop:
            break
    return i+1
            

def load_web_snap_shot_eigs_():
    '''
    Since websnapshot eigenpairs are too big to be handled by cPickle, this specialized function was used to wrap np.savez
    '''

    # np.savez_compressed("WebSnapShot_all_Eigs", evals, evecs)    #how they were saved
    eigsWhere = "../Data/Spectrums/SM/WebSnapShot_all_Eigs.npz"
    data = np.load(eigsWhere)
    evals = data['arr_0']
    evecs = data['arr_1']    

    return evals, evecs


def eigen_pairs(inGraph, graphName, eigs, strategy=None, balance=None):
    '''
    Returns:
                evecs    -    (num_eigenvectors X num_nodes) matrix carrying the Laplacian eigenvectors.
    '''
    if graphName == "WebSnapShot" and eigs == "all":
        evals, evecs = load_web_snap_shot_eigs_()
    else: # Load "all" the eigenpairs and trim later.
        evals, evecs = gk.laplacian_spectrum(inGraph, graphName, "all", edgeW=None, tryLoad=True, save=True)
    
    assert(is_increasing(evals))
    
    if eigs == "all" and strategy != "derivative":
        return evals, evecs

    if type(eigs) != str and eigs <= 1:                 # Adjust eigenvalues is a percent [0, 1] is requested.
        totalEigs = len(evals)
        eigs      = int(np.ceil(eigs * totalEigs))   
        
    if strategy == "small" or strategy == None:        
        evals, evecs       = gk.filterEigenPair(evals, evecs, keep=eigs, magnitude="small")
    elif strategy == "big":
        evals, evecs       = gk.filterEigenPair(evals, evecs, keep=eigs, magnitude="big")    
    elif strategy == "derivative":
        eigsToKeep         = pick_eigs_exhaust_derivative_(evals)
        print "Eigenvalues kept " +str(eigsToKeep)
        evals, evecs       = gk.filterEigenPair(evals, evecs, keep=eigsToKeep, magnitude="small")
    elif strategy == "both_ends":
        evals, evecs       = gk.filterEigenPair(evals, evecs, keep=eigs, magnitude="both", balance=balance)
    else:
        assert(False)

    return evals, evecs






################################################################################
######################## METHODS WHICH MIGHT DIE FOLLOW ########################
################################################################################



def panos_sim_kavla(inGraph, graphName, eigs="all", energies=10, strategy=None):
        
    evals, evecs = eigen_pairs(inGraph, graphName, eigs, strategy)        
    degrees      = inGraph.degree_property_map("total").a    
    n            = evecs.shape[1]
    
    dists        = np.zeros((n,n))

    for i in xrange(energies):    
        norm_eigenV  = evecs[i,:] * (degrees - evals[i])
        dists += np.abs( np.subtract.outer(norm_eigenV, norm_eigenV))
                 
    return dists     



def spectralSim_PP_old(inGraph, graphName, eigs="all", energies=100, strategy=None):
    evals, evecs       = gk.laplacian_spectrum(inGraph, graphName, eigs, tryLoad=True)
    timeSamples, sigma = gk.wks_energy_generator(evals[0], evals[-1], energies)
    sig                = gk.WKS(evals, evecs, timeSamples, sigma)

    sig2               = gk.neighbor_signatures(sig, inGraph, "average") 
    sig                = np.hstack((sig, sig2)); 

    sig3               = gk.neighbor_signatures(sig, inGraph, "average")
    sig = np.hstack((sig, sig3)); 

    sig4               = gk.neighbor_signatures(sig, inGraph, "average")
    sig = np.hstack((sig, sig4)); 
     
    sig5               = gk.neighbor_signatures(sig, inGraph, "average")
    sig = np.hstack((sig, sig5)); 
                    
    sig2               = gk.neighbor_signatures(sig, inGraph, "average")
    sig3               = gk.neighbor_signatures(sig2, inGraph, "average")
    sig4               = gk.neighbor_signatures(sig3, inGraph, "average")
    sig5               = gk.neighbor_signatures(sig4, inGraph, "average")
 
    
    hs = np.hstack
    sig                = hs ((sig, hs((sig2, hs((sig3, hs((sig4, sig5))))))))
    print sig.shape

    return sig


def spectralSim_PP(inGraph, graphName, eigs="all", energies=100, strategy=None):
    evals, evecs       = eigen_pairs(inGraph, graphName, eigs, strategy)
    timeSamples, sigma = gk.wks_energy_generator(evals[0], evals[-1], energies)
    sig                = gk.WKS(evals, evecs, timeSamples, sigma)
    sig                = neighbor_signature_propagation(sig, inGraph, 1)
    return sig

def heatSim_PP(inGraph, graphName, eigs="all", energies=100, strategy=None, shrink=1):        
    evals, evecs       = eigen_pairs(inGraph, graphName, eigs, strategy)
    timeSamples        = gk.hks_time_sample_generator(evals[0], evals[-1], energies, shrink=50)
    sig                = gk.HKS(evals, evecs, timeSamples)
    sig                = neighbor_signature_propagation(sig, inGraph, 1) 
    return sig

def neighbor_signature_propagation(initialSig, inGraph, iterations):
    result = initialSig

    for i in xrange(iterations):
        sigNew    = gk.neighbor_signatures(initialSig, inGraph, "average")        
        result    = np.hstack((result, sigNew))
        initialSig = sigNew
        
    return result

def remove_not_varying_parts(A, thres):
    ''' Clear the ending portion of vectors that appears to not vary.
        Input:
                A      -     (m X n) array like. Lines encode the variables and columns the features.
                thres  -     (float)
        Output:  
                res    -     A list of numpy vectors. res[i] corresponds to the i-th line of A, with the trailing not varying features 
                             removed.
    '''
    m, n = A.shape
    assert(n > 2)
    res  = []
    for row in xrange(m):
        stopped = False
        diffOld = abs(A[row,1] - A[row,0])               
        for col in xrange(1, n-1):
            diffNew = abs(A[row, col] - A[row, col+1])
#             print diffOld, diffNew
#             raw_input()        
#             if abs(diffNew - diffOld) < thres * diffOld:
            if abs(diffNew - diffOld) < thres:
                res.append(np.array(A[row, 0:col]))
                stopped = True
#                 print col
                break;
            else:
                diffOld = diffNew        
        if stopped == False:
            res.append(np.array(A[row,:]))
    return res

