'''
Created on Jul 15, 2014

@author: optas
'''

from graph_analysis import IO
from graph_analysis import basic_tools as myGT
from role_discovery import graph_kernels as gk
import multiprocessing
import time



#1. where to save the returned evals/evecs?
def egoSpectra(inGraph, node, k, tryLoad):
    egoNet        = myGT.egoNetwork(inGraph, node)    
    return        gk.laplacian_spectrum(egoNet, "", k,  edgeW=None, tryLoad=tryLoad, save=False)


    

def parallel_ego_signatures(inGraph):
#     multiprocessing.Process(target=egoSpectra, args=(inGraph, inGraph.vertices().next(), "all", False)).start()    
    for count, node in enumerate(inGraph.vertices()):
        multiprocessing.Process(target=egoSpectra, args=(inGraph, node, "all", False)).start()
#         if not count % 50: 
#             time.sleep(10)
        

import sys

if __name__ == "__main__":
    inGraph = graph_analysis.IO.load_GT_graph("karate", gcc=True, removeSL=True, removePE=True)
    megaB = 1024 *1024
    print sys.getsizeof(inGraph) / float(megaB)
    parallel_ego_signatures(inGraph)
 
    
#     
#     timeSamples, sigma = gk.waveTimeSample3(evals[0], evals[-1], energies, shrink)            
# #     timeSamples, sigma = gk.waveTimeSample(evals[0], evals[-1], energies)
# #     timeSamples, sigma = gk.waveTimeSample(evals[-1], evals[0], energies)
#     sig                = gk.WKS(evals, evecs, timeSamples, sigma)     
#     




# for node in inGraph.vertices():
#     multiprocessing.Process(target=egoSpectra, args=(inGraph, node, "all", False)).start()
    
    
    
    

