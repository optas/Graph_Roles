'''
Created on May 29, 2014

@author: optas
'''

import sys
from scipy import spatial as spt
import numpy as np
from role_discovery import graph_kernels as gk
import matplotlib.pyplot as plt    
from myUtils import pdist_to_pRank, inner_intra_distances
from graph_tool import generation
from graph_tool import Graph, draw
from graph_analysis.basic_tools import addSymmetricEdge
import myPloting as myPlt
from myUtils import p
from graph_tool.libgraph_tool_core import new_vertex_property
from _collections import defaultdict


def burtFig4(directed=False):
    '''
    Returns the graph presented at Burt's "Role Equivalence" paper, at Figure_x
    '''
    g = Graph(directed=directed)
    g.add_vertex(9)
    g.add_edge(0, 4); g.add_edge(1, 4); g.add_edge(2, 4); g.add_edge(3, 4)
    g.add_edge(4, 5); g.add_edge(4, 6)
    # 4-clique
    g.add_edge(5, 6); g.add_edge(5, 7); g.add_edge(5, 8)
    g.add_edge(6, 7); g.add_edge(6, 8);
    g.add_edge(7, 8)
    if directed:
        g.add_edge(6, 5); g.add_edge(7, 5); g.add_edge(8,5)
        g.add_edge(7, 6); g.add_edge(8, 6);
        g.add_edge(8, 7)    

    for ed in g.edges():
        print ed
    return g


def burtFig5(directed=False):
    '''
    Returns the graph presented at Burt's "Role Equivalence" paper, at Figure_x
    '''
    g = Graph(directed=directed)
    g.add_vertex(10)
    if directed:
        addSymmetricEdge(g, 0, 1)
        addSymmetricEdge(g, 1, 2)
        addSymmetricEdge(g, 2, 3)
        addSymmetricEdge(g, 4, 6)
        addSymmetricEdge(g, 6, 5)
        addSymmetricEdge(g, 7, 9)
        addSymmetricEdge(g, 9, 8)
        g.add_edge(4, 0); g.add_edge(5, 0);
        g.add_edge(7, 3); g.add_edge(8, 3);
    else:
        g.add_edge(0, 1); g.add_edge(1, 2); g.add_edge(2, 3); g.add_edge(4, 6);
        g.add_edge(6, 5); g.add_edge(7, 9); g.add_edge(9, 8); g.add_edge(4, 0); 
        g.add_edge(5, 0); g.add_edge(7, 3); g.add_edge(8, 3);
        
    return g
    

def batonGraph(drawGraph=False):    
    '''
    A graph with NO symmetries. 
    '''
    g = generation.lattice((9, 1))
    v = g.add_vertex()
    g.add_edge(v, g.vertex(6))    
    
    if drawGraph == True:
        draw.graph_draw(g, vertex_text=g.vertex_index, edge_color="black", output="batton.pdf")        

    return g
    
    
def balancedBinaryTree(h, drawGraph=False):
    '''
    h - the height of the tree
    '''
    g = Graph(directed=False)
    g.add_vertex(2**h - 1)    
    
    for i in xrange(1, 2**(h-1)):
        lc = 2*i - 1
        rc = lc + 1
        g.add_edge(i-1, lc)
        g.add_edge(i-1, rc)
    
    hIndex = g.new_vertex_property("int") #associate with each node the height at which it lies
    k = 2; m = 1;
    for i in xrange(1, len(hIndex.a)):
        hIndex.a [i] = m;
        k -= 1; 
        if k==0: m += 1; k = 2**m             
    g.vp['height'] = hIndex
    
    if drawGraph == True:
        draw.graph_draw(g, vertex_text=g.vertex_index, edge_color="black", output="binaryTree_h_"+str(h)+"_.pdf")        

    return g


def starGraph(numOfVertices):
    '''
    numOfVertices = Number of vertices including the center node of the star.
    dependencies: graph_tool
    '''
    g = Graph(directed=False)
    centerNode = g.add_vertex()
    for i in xrange(2, numOfVertices+1):
        node = g.add_vertex()
        g.add_edge(centerNode, node)            
    return g



import presentation as mc
if __name__ == "__main__":
    inGraph                = balancedBinaryTree(h=6, drawGraph=False)
    graphName              = "balanced_tree_h_6"
    allEvals, allEvecs     = gk.laplacian_spectrum(inGraph, graphName, "all", tryLoad=False, save=True)

    mc.optimize_class_frequencies(graphName, minEig=2, maxEig=10, energies=29, taxon=3)    
#     
# 
#     innerAv = []
#     intraAv = []
#     print allEvals
#     
#     for throw in xrange(len(allEvals)-1):        
#         evals, evecs                         = gk.filterEigenPair(allEvals, allEvecs, keep=len(allEvals)-throw, magnitude="small")
#         
#         tols = 1e-7
#         evecs[np.abs(evecs)<tols] = 0        

            
#         timeSamples, sigma            = gk.waveTimeSample2(evals[0], evals[-1], timePoints=10)
#         sig                           = gk.WKS(evals, evecs, timeSamples, sigma)
#         dists                         = spt.distance.pdist(sig, 'canberra')       
#         distanceMatrix                = spt.distance.squareform(dists)        
#         innerDists, intraDists        = inner_intra_distances(distanceMatrix, groupTaxa=inGraph.vp['height'].a, blackList=None)
#  
#         innerDists[innerDists < 1e-6] = 0
#         intraDists                     = np.sum(spt.distance.squareform(intraDists),1) / (len(innerDists)-1)
#                  
#         innerAv.append(np.average(innerDists))
#         intraAv.append(np.average(intraDists))
#          
#         if len(allEvals) - throw == 11:
#             myPlt.heatMap(distanceMatrix, save=False, saveAt="9e.pdf")
#             p("yo")
# #     distanceMatrix     = pdist_to_pRank(distanceMatrix)
#     print distanceMatrix
         
 
#         print innerAv
#         print intraAv
#         print [k/float(m) for k, m in zip(innerAv, intraAv)]    
#         print range(2, len(allEvals)+1)[::-1]
#         
# 
# #     YtoPlot = [innerAv, intraAv, [k/float(m) for k, m in zip(innerAv, intraAv)]]    
#     YtoPlot = innerAv
#     myPlt.x_Vs_y("innerAverage", range(2, len(allEvals)+1)[::-1], YtoPlot, "NumOfEigs", 'Toy tree', legendLabel=['Inner_Average', 'Intra_Average', 'Inner/Intra'], stamp="-",
#                  save=True, saveAt="tree_toyinnerAverage")
# 
#     


#     myPlt.heatMap(distanceMatrix, save=False, saveAt="")
#     distanceMatrix     = pdist_to_pRank(distanceMatrix)
#     print distanceMatrix

#     waveNorm = inGraph.new_vertex_property("double")
#     for i, sig in enumerate(sig):
#         waveNorm[inGraph.vertex(i)]  = np.linalg.norm(sig)
#         


    sys.exit("Exiting Successfully.")