'''
Created on Feb 16, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

from .node_similarities import *
from clustering.distances import pdist_to_p_rank
from scipy import spatial as spt
from EMD import calc_emd
from clustering.distances import jaccard_vectorial
import warnings
import graph_tool as gt
import graphs_with_roles as roles 
from basic_utils import is_increasing

def is_spectral(methodName):
    '''
    returns True iff the method designated by the methodName is using the spectrum of the Adjacency/Laplacian matrix.
    '''
    if methodName in ["wave_sim", "heat_sim",  "spectralSim_PP", 'heatSim_PP',  'spectral_degree_sim', "raw_laplacian_sim"]: return True 
    if methodName.startswith("panos"): return True
    
    return False

def gives_embedding(methodName):
    '''
    Returns true iff the corresponding method can generate an embedding for each node.
    '''
    if methodName in ["sim_rank", "role_sim", "vertex_sim", "commute_time_dist"]: return False
    return True


def method_dispatcher(methodName):

    return globals()[methodName]
    
    if methodName == "heat_sim"         : return heat_sim    
    if methodName == "wave_sim"         : return wave_sim
    if methodName == "role_sim"         : return role_sim
    if methodName == "sim_rank"         : return sim_rank
    if methodName == "refex_sim"        : return refex_sim
    if methodName == "vertex_sim"       : return vertex_sim
    if methodName == "cycle_sim"        : return cycle_sim
    
    if methodName == "panos_sim"        : return panos_sim        
    if methodName == "panos_sim_2"      : return panos_sim_2
    if methodName == "panos_sim_3"      : return panos_sim_3
    if methodName == "panos_sim_4"      : return panos_sim_4

    if methodName == "commute_time_dist"      : return commute_time_dist            
    if methodName == "raw_laplacian_sim"      : return raw_laplacian_sim
    assert(False)    


def default_distance_function(methodName):
    if methodName == "heat_sim"             : return 'euclidean'
    if methodName == "wave_sim"             : return 'canberra'    
    if methodName == "refex_sim"            : return 'euclidean'    
    if methodName == "cycle_sim"            : return 'cityblock'      
    if methodName == "raw_laplacian_sim"    : return 'euclidean'    
                  
    if methodName == "panos_sim"            : return 'correlation'
    if methodName == "panos_sim_2"          : return 'correlation'
    if methodName == "panos_sim_3"          : return 'correlation'
    if methodName == "panos_sim_4"          : return 'correlation'
    assert(False)


def pairwise_distances_(signatures, methodName, distFunction="default"):
    '''
    Computes the pairwise distances between the given signatures.
    '''

    print "Computing distances with " + str(distFunction) + "."
        
    if distFunction == calc_emd or distFunction == jaccard_vectorial:
        #These functions do not work (yet) with pdist
        n = len(signatures)
        allDists = np.empty((n*(n-1)/2))
        m = 0
        for i in xrange(n):    
            for j in xrange(i+1, len(signatures)):
                allDists[m] = distFunction(signatures[i], signatures[j])
                m += 1         
        return spt.distance.squareform(allDists)
         
    if distFunction == "default" or distFunction == None:
        distanceMatrix = spt.distance.pdist(signatures, default_distance_function(methodName))
    else:        
        distanceMatrix = spt.distance.pdist(signatures, distFunction)
            
    distanceMatrix     = spt.distance.squareform(distanceMatrix)
    distanceMatrix[abs(distanceMatrix)<1e-9] = 0.0;      # Safeguarding against floating point numerical errors
    return distanceMatrix
    
    
def execute_method(methodName, inGraph, graphName, distances=True, distFunction="default", methodParams="default"):
    print "Executing " + methodName
    method = method_dispatcher(methodName)

    if methodName.startswith("panos_sim"):
        if methodParams == "default" or methodParams == None:
            signatures  = method(inGraph, graphName)            
        else:
            eigs, strategy, compressed, sigma = methodParams
            signatures  = method(inGraph, graphName, eigs, strategy, compressed, sigma)

    elif is_spectral(methodName):
        if methodParams == "default" or methodParams == None:
            signatures  = method(inGraph, graphName)
        else:
            eigs, energies, strategy = methodParams
            print "Using " +str(eigs) +" eigenvalues for a %d-dimensional embedding via the eigenvalue strategy \'%s\'" % (energies, strategy)
            signatures   = method(inGraph, graphName, eigs, energies, strategy)

    elif methodName == "refex_sim":
        signatures = method(graphName)
        
    elif methodName == "cycle_sim":                        
        cycleLength, matrixType, keepEven, clearConverged = methodParams
        if clearConverged and distances:
            assert(distFunction == calc_emd)
        signatures = method(inGraph, cycleLength, matrixType, keepEven, clearConverged)

    elif methodName in ["sim_rank", "role_sim", "vertex_sim", "commute_time_dist"]:
        if distances == False: # These methods only return a distance per pair of nodes.
            warnings('%s only provides a distance matrix.' %(methodName, ))
        if methodName == "vertex_sim":
            distanceMatrix = method(inGraph, dist=True)    
        elif methodName == "commute_time_dist":
            distanceMatrix = method(inGraph)
        else:
            distanceMatrix = method(graphName)
        distanceMatrix[abs(distanceMatrix)<1e-9] = 0.0;      # Safeguarding against floating point numerical errors.
        return distanceMatrix
    else:
        assert(False)
          
    if distances:                        
        distanceMatrix = pairwise_distances_(signatures, methodName, distFunction=distFunction)        
        return distanceMatrix
             
    else:
        return signatures


def prepare_input_graph(graphName, metric, verbose=True):
    storedFolder            = roles.graph_folder(graphName)
    try:
        inGraph             = IO.load_data("../Data/Graphs/"+storedFolder+"/"+graphName + ".GT.graph").next()
    except:
        inGraph             = gt.load_graph("../Data/Graphs/"+storedFolder+"/"+graphName + ".graph.xml")

    if verbose:    
        print "Loaded Input Graph.\nName = %s. \nMetric = %s. \n#Nodes = %d. #Edges = %d." % (graphName, metric, inGraph.num_vertices(), inGraph.num_edges())
    groupTaxa, blackList    = roles.graph_node_clusters(graphName, inGraph, metric)
    
    xTickMarks              = roles.taxa_names(graphName)
    
    if verbose:
        if blackList != None:
            print "True Number of Clusters = " +str(len(set(groupTaxa)) - len(blackList)) + "\n"
        else:
            print "True Number of Clusters = " +str(len(set(groupTaxa))) + "\n"
    
    return inGraph, groupTaxa, blackList, xTickMarks
    
# 
# def extensive_distance_based_comparison(graphName, allMethods = "all", ranks=False, parallel=False):
#     if allMethods == "all":
#         allMethods = ["roleSim", "refexSim", "simRank", "spectralSim", "spectralSim_PP", "heatSim", "spectralProjectionSim"]
#     
#     if parallel:
#         sys.stdout = open("../Output/"+graphName+"_distanceBased_comparison", "w")
#         
#     printFancy(graphName, "-")
#     inGraph, groupTaxa, blackList, xTickMarks = prepare_input_graph(graphName, metric="default")
#     gDiameter = topology.pseudo_diameter(inGraph)[0]
#     edgeDensity = inGraph.num_edges() / (2*float(inGraph.num_vertices()))
#     clusterCoef = clustering.global_clustering(inGraph)[0]
#     
#     print "clusterCoef = %0.3f" %(clusterCoef, )
#     print "Pseudo-diameter = ", gDiameter
#     print "Edge density = %0.3f" % (edgeDensity,)
#         
#     for methodName in allMethods:
#         print "---------------------------"
#         print "\t" + methodName
#         print "---------------------------"
#         
#         if is_spectral(methodName):
#             eigs = "all"
#             strategy = None;
#             if inGraph.num_vertices() > 5000:
#                 energies = 35
#             else:
#                 energies = int (gDiameter)
#                 print energies, int(1*gDiameter)                
# 
# #             energyMultiple = 3
# #             energies = int(energyMultiple * int(gDiameter))
#             methodParams = [energies, "all", strategy]            
#             print methodParams
# #             print "\nNumber of Energies = %d X %s" % (energyMultiple, " Diameter.")            
#             distMatrix = execute_method(methodName, inGraph, graphName, distances=True, ranks=ranks, distFunction="default", methodParams=methodParams)
#             
# #             myPlt.plt.matshow(distMatrix)
# #             myPlt.plt.show()
#                         
#             evaluate_distance_matrix(distMatrix, groupTaxa, "affinity", ranks)
#             inner_intra_distances(distMatrix, groupTaxa, blackList=None, ranks=ranks)
#             clustering_c_index(distMatrix, groupTaxa)
#             
# #         if methodName == "raw_laplacian_sim":
# #                     break
# 
#         else: #non Spectra methods            
#             distMatrix = execute_method(methodName, inGraph, graphName, distances=True, ranks=ranks, distFunction="default", methodParams="default")
#             evaluate_distance_matrix(distMatrix, groupTaxa, "all", ranks)
#             inner_intra_distances(distMatrix, groupTaxa, blackList=None, ranks=ranks)
#             clustering_c_index(distMatrix, groupTaxa)
#     
#     
# #     if parallel:
# #         myUtils.emailPanos(" ", graphName +"-distance based- done.")

