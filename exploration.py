'''
Created on Dec 31, 2014

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

from graph_tool import topology
from clustering.distances import exclude_black_listed
from clustering.evaluation import evaluate_distance_matrix, cluster_distance_matrix, evaluate_unsup_clustering
import numpy as np
import sys
from role_discovery.comparisons import pairwise_distances_, is_spectral
from graph_analysis.node_clustering import count_edges_between_clusters, charachteristics_of_node_clusters, average_signature_of_clusters
from role_discovery.graphs_with_roles import plot_roles_on_graph
from role_discovery.node_similarities import eigen_pairs, commute_time_dist, neighbors_matter, vertex_sim
from graph_analysis.IO import compute_or_load, save_data, load_data
from numpy.linalg import svd
import matplotlib.pyplot as matplt
from basic_utils import p
from role_discovery.comparisons import prepare_input_graph, execute_method
from nodeClusteringResults import nodeClusteringResults
from basic_utils import print_fancy

allGraphs  = ["serengeti-foodweb", "Carribean_FoodWeb", "ml_giant",  'kd_giant', 'db_giant', 
             "crcl_1", "manufacturing_1", "manufacturing_2", "StMarks", "polblogs", "celegansneural", 
             "netscience", "adjnoun", "polbooks", "polblogs", "Sawyer_English", "12_Years_Slave", "E_Coli", "Mapk",
             "Plato_Republic", "Zarathustra_English", "stratified_1", "crcl_1k_5k_10ccl_20_nodes_each_0", "crcl_1k_5k_10ccl_20_nodes_each_1", "WebSnapShot"]


someMethods = [ "wave_sim", "heat_sim", "role_sim", "sim_rank", "refex_sim", "cycle_sim", "panos_sim", "vertex_sim", "commute_time_dist"]


#             if methodName == "refex_sim":
#                 methodParams = None
# 
#             elif methodName.startswith("panos"):
#                 methodParams = [0.05, "Big", False, None]
#                 print methodParams 
# 
#             else:
#                 assert(False)


if __name__ == '__main__':    

    methodName = sys.argv[1]    
    eigs       = sys.argv[2]
    
    if eigs != "all":
        eigs = float(eigs)
        assert(eigs < 1 and eigs >0 )
    
    strategy  = sys.argv[3]
    
    print methodName
    
    for graphName in allGraphs:
        try:
            print "\n\n"        
                
            inGraph, trueClusters, blackList, xTickMarks = prepare_input_graph(graphName, metric="default")
            graphDiameter    = int(topology.pseudo_diameter(inGraph)[0])
            trueClustersSize = len(np.unique(trueClusters))

            methodParams = [eigs, graphDiameter, strategy]                        # SET MANUALLY
            print methodParams
                
            fileName         = "../Data/Sigs/" + methodName + "/" + graphName + "_" + '_'.join( str(e) for e in  methodParams)            
            signatures       =  compute_or_load(execute_method, fileName + "_sigs",  True, methodName, inGraph, graphName, distances=False, distFunction="default",  methodParams=methodParams)
            distanceMatrix   =  pairwise_distances_(signatures, methodName, distFunction="default")
            
            if blackList != None:
                trueClustersSize -= len(blackList)
                distanceMatrix, trueClusters  = exclude_black_listed(distanceMatrix, blackList, trueClusters)
            
            res = nodeClusteringResults(methodName, graphName)
            
            for clusteringMethod in ["spectral", "affinity", "dbscan", "hierarchical-complete", "hierarchical-average", "hierarchical-single"]:
          
                if clusteringMethod.startswith("hierarchical"):
                    paramKey =  clusteringMethod.split("-")[1]
                    clusteringMethod = "hierarchical"
                else:
                    paramKey = "default"
          
                predClusters                 = cluster_distance_matrix(distanceMatrix, clusteringMethod, linkage=paramKey, clusterNum=trueClustersSize)
                homo, comp, vmea, aran, amut = evaluate_unsup_clustering(trueClusters, predClusters)
             
                res.add_clustering_labels(clusteringMethod, paramKey, predClusters)
                res.add_clustering_evaluation(clusteringMethod, paramKey,  {"Homogeneity" : homo, "Completeness": comp, "V-Measure": vmea, "A-Rand":aran, "A-Mutual_Information": amut})
                       
            save_data("../Output/"+methodName + "_" + '_'.join(str(e) for e in  methodParams) + "_" + graphName, res)
    
        except Exception as e:
            print e
            continue
     
    sys.exit("Success.")
