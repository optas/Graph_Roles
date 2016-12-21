'''
Created on Apr 26, 2014

@author: optas
'''

import multiprocessing
import massiveComparisons as mc
import os


allGraphs = ["ml_giant",  'kd_giant', 'db_giant', "serengeti-foodweb", "crcl_1", "manufacturing_1", "manufacturing_2", "polbooks", "StMarks",
             "Carribean_FoodWeb", "polblogs", "netscience", "celegansneural", "Sawyer_English", "12_Years_Slave", "WebSnapShot"]



if __name__ == '__main__':
    parallel = True
    
#     Params for     mc.extensive_distance_based_comparison
#     ranks      = False    
#     allMethods = "all";
#     params = [[graph, allMethods, ranks, parallel]  for graph in allGraphs]
        
    allMethods = ["roleSim", "spectralSim", "heatSim", "simRank"]
    params = [[graph, allMethods, parallel]  for graph in allGraphs]
     
    for i in xrange(len(params)):
        p = multiprocessing.Process(target=mc.extensive_embedding_based_comparison, args=params[i]).start()
