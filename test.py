'''
Created on Jan 2, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

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
from basic_utils import p, print_fancy, relative_to_worst
from role_discovery.comparisons import prepare_input_graph, execute_method, gives_embedding
from nodeClusteringResults import nodeClusteringResults
from numpy import linalg
from matplotlib.mlab import PCA


allGraphs   = ["serengeti-foodweb", "Carribean_FoodWeb", "ml_giant",  "kd_giant", "db_giant",
             "crcl_1", "manufacturing_1", "manufacturing_2", "StMarks", "polblogs", "celegansneural", 
             "netscience", "adjnoun", "polbooks", "polblogs", "Sawyer_English", "12_Years_Slave", "E_Coli", "Mapk",
             "Plato_Republic", "Zarathustra_English", "stratified_1", "crcl_1k_5k_10ccl_20_nodes_each_0", "crcl_1k_5k_10ccl_20_nodes_each_1", "WebSnapShot"]


currentGraphs = allGraphs
 


someMethods = [ "wave_sim", "heat_sim", "role_sim", "sim_rank", "refex_sim", "cycle_sim", "panos_sim", "vertex_sim", "commute_time_dist"]

currentMethods = ["panos_sim", "sim_rank", "role_sim", "vertex_sim", "refex_sim", "commute_time_dist"]

if __name__ == '__main__':    
    
    for graphName in currentGraphs:
        print_fancy(graphName)
        v_measAllMethods = []
        for methodName in currentMethods:
            
            print "\n\n"                
            inGraph, trueClusters, blackList, xTickMarks = prepare_input_graph(graphName, metric="default", verbose=False)
            graphDiameter    = int(topology.pseudo_diameter(inGraph)[0])
            trueClustersSize = len(np.unique(trueClusters))
                         
            if methodName.startswith("panos"):
                methodParams = ["all", "small", False, None]
#                 print methodParams 
            else:     
                methodParams = None
                        
            
            try:
                if methodParams == None:
                    res = load_data("../Output/"+methodName +"_" + graphName).next()
                else:     
                    inFile = "../Output/" + methodName +"_"  + '_'.join( str(e) for e in  methodParams) + "_" + graphName                
                    res = load_data(inFile).next()
            except:
                print methodName + " None!"
                continue
            
            for clusteringMethod in ["spectral"]: #, "affinity", "dbscan", "hierarchical-complete", "hierarchical-average", "hierarchical-single"]:                    
                if clusteringMethod.startswith("hierarchical"):
                    paramKey =  clusteringMethod.split("-")[1]
                    clusteringMethod = "hierarchical"
                else:
                    paramKey = "default"    
    
#             print res.get_clustering_labels(clusteringMethod, paramKey)      
#                 print "Homogeneity",  res.get_clustering_evaluation(clusteringMethod, paramKey,  "Homogeneity")
#                 print "Completeness", res.get_clustering_evaluation(clusteringMethod, paramKey,  "Completeness")
#                 print "V-Measure",    res.get_clustering_evaluation(clusteringMethod, paramKey,  "V-Measure")

            v_measAllMethods.append(res.get_clustering_evaluation(clusteringMethod, paramKey,  "V-Measure"))
            
        print relative_to_worst (v_measAllMethods)


    sys.exit("Success.")
          



#       

 


# 
# 
# 
# 
# 
#         
# 
#   
#     count_edges_between_clusters(inGraph, trueClusters)
#     charachteristics_of_node_clusters(inGraph, trueClusters)
#     nd1 = neighbors_matter(inGraph, distanceMatrix)
#     plot_roles_on_graph(inGraph, graphName, metric="default")
#     aveSig = average_signature_of_clusters(signatures, trueNodeClusters)
#     
#     
#     
# 
# #     for sig in aveSig.keys():
# #         print aveSig[sig]
# 
# #     U, S, V = svd(signatures, full_matrices=False, compute_uv=1)
# #     print S
# #     S[10:] = 0
# #     S = np.diag(S)
# #     signaturesRed  = np.dot(U, np.dot(S, V))
# 




#     from sklearn.cluster import KMeans
#     features  = execute_method(methodName, inGraph, graphName, distances=False, ranks=False, distFunction=distFunction, methodParams=[eigs, energies, strategy])
#     estimator = KMeans(init='k-means++', n_clusters=4, n_init=10)
#     estimator.fit(features)                   
#     evaluate_unsup_clustering(trueClusters, estimator.labels_, n_clusters=4, verbose=True)
    
