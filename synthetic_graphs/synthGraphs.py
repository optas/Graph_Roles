'''
Created on May 9, 2014

@author: optas
'''

import warnings, sys
from graph_analysis import IO
import numpy as np
from role_discovery import graph_kernels as gk
import myPloting as myPLT
from matplotlib.mlab import PCA
from graph_analysis.basic_tools import make_simple_graph

try:
    from graph_tool import Graph, topology, spectral, generation
    from graph_tool.all import random_graph, random_rewire
except:
    print "Graph_Tool is not installed in this environment."
     
def erdos_renyi_graph(n, e, directed=False, gcc=True):        
    g = Graph(directed=directed)
    g.add_vertex(n)
    
    rint = np.random.randint
    edge_list = [[x,y] for x,y in zip(rint(0, n, size=e), rint(0, n, size=e))]
    g.add_edge_list(edge_list)        

    random_rewire(g, model="erdos")      
    g = make_simple_graph(g, undirected=1-directed, gcc=gcc)    
    return g
                   

def erdos_circles_cliques(bbSize, nCircles, circleSize, nCliques, cliqueSize, save=True):
    '''
    Makes a random (Erdos-Renyi) graph connected with circles and cliques.
    Input:
        -bbSize[0]-   how many nodes the random Graph will have
        -bbSize[1]-   how many edges the random Graph will have
        -nCliques-    how many cliques to add of size -cliqueSize-
        -nCircles-    how many circles to add of size -circleSize-    
    '''

    backbone = erdos_renyi_graph(bbSize[0], bbSize[1], directed=False, gcc=True)
    print "The random backbone graph has %d vertices and %d edges." % (backbone.num_vertices(), backbone.num_edges())    
    
    assert( np.sum(topology.label_largest_component(backbone).a) ==  backbone.num_vertices())
    
    if backbone.num_vertices() < nCircles + nCliques:        
        warnings.warn("The erdos part of the graph is too small to add the requested circles/cliques.")
        
    
    backbone.vp["bb"]    = backbone.new_vertex_property("short")         # Mark all nodes belonging in the random backbone.
    backbone.vp["bb"].a  = np.ones(backbone.num_vertices())
            
    gCircle              = add_in_graph(None, nCircles, generation.circular_graph, circleSize)
    gCircle.vp["circ"]   = gCircle.new_vertex_property("short")          # Mark all nodes belonging in the Circles.
    gCircle.vp["circ"].a = np.ones(gCircle.num_vertices())
    
    gCliq                = add_in_graph(None, nCliques, generation.complete_graph, cliqueSize)
    gCliq.vp["cliq"]     = gCliq.new_vertex_property("short")            # Mark all nodes belonging in the Cliques.
    gCliq.vp["cliq"].a   = np.ones(gCliq.num_vertices())
    
    concat1 = generation.graph_union(backbone, gCliq, internal_props=True)
    gFinal  = generation.graph_union(concat1, gCircle, internal_props=True)

    assert (sum(gFinal.vp['cliq'].a == 1) == gCliq.num_vertices() and \
            sum(gFinal.vp['circ'].a == 1) == gCircle.num_vertices())
            
    comp, hist  = topology.label_components(gFinal)                     # Mark to which CC every node is.    
    numOfCC     = max(comp.a)
    assert(numOfCC == nCircles + nCliques)
  
    bbNodes = np.where(gFinal.vp["bb"].a == 1)[0]
    np.random.shuffle(bbNodes); k = 0;
    
    gFinal.vp["attachments"] = gFinal.new_vertex_property("short")      # Bookkeeping which nodes of the backbone where used to connect with the circles/cliques.

    for cc in range(1, numOfCC+1):
        atNode = np.where(comp.a == cc)[0][0]  # Since all nodes of the added graphs are equivalent we can pick the 1st to make the attachment.
        gFinal.add_edge(atNode, bbNodes[k]); k+=1
        gFinal.vp["attachments"].a[atNode] = 1
        gFinal.vp["attachments"].a[bbNodes[k]] = 1      
    
    assert (topology.label_components(gFinal)[1][0] == gFinal.num_vertices()) # gFinal must be Fully Connected

    print "The graph with the cliques and circles has in total %d vertices and %d edges." % (gFinal.num_vertices(), gFinal.num_edges())
    
    return gFinal
    
    
def add_in_graph(gInit, n, graphGen,  *args):
    '''
    Add into initial input graph -gInit-, -n- disconnected graphs that are produced by the graphGen(*args) function.
    Input:
        -gInit- initial graph. If equals to None, then the output Graph will be comprised by -n- disconnected components each produced by graphGen(*args).
        -n- the number of graphs that we will add to gInit.
        -graphGen- a function that is producing a specific type of graph, e.g. create_star_graph
        -*args- the arguments with witch graphGen will be called.
        
    Example:
        Calling add_in_graph(None, 10, generation.complete_graph, 5)
        produces a graph that is comprised by 10 complete graphs of size 5 each.        
    '''

    if   n <= 1 and gInit != None : return generation.graph_union(graphGen(*args), gInit)
    elif n <= 1 and gInit == None : return graphGen(*args)    
    else: 
        return generation.graph_union(graphGen(*args), add_in_graph(gInit, n-1, graphGen, *args))
    
    
def append_lattices(g, nLatte = 20, latticeSize=[5,5]):
    def isCorner(gLatte, node):
        if gLatte.vertex(node).out_degree() == 2: return True
        else: return False
            
    gLatte = add_in_graph(None, nLatte, generation.lattice, latticeSize)
    gLatte.vp["late"]   = gLatte.new_vertex_property("short")
    gLatte.vp["late"].a = np.ones(gLatte.num_vertices())
    
    gFinal = generation.graph_union(g, gLatte, internal_props=True)
    comp, hist  = topology.label_components(gFinal)    
    numOfCC     = max(comp.a)
    assert(numOfCC == nLatte)
    
    bbNodes  = set(np.where(gFinal.vp["bb"].a == 1)[0])             
    attached = set(np.where(gFinal.vp["attachments"].a == 1)[0])
    freeBBnodes = bbNodes - attached
    if freeBBnodes < nLatte:
        warnings.warn("Not enough free nodes in the erdos graph to attach all the lattices.")
        return None
    
    freeBBnodes = list(freeBBnodes)
    np.random.shuffle(freeBBnodes);     
    k = 0;        
    for cc in range(1, numOfCC+1):        
        atNode = np.where(comp.a == cc)[0][0] 
        assert(isCorner(gFinal, atNode))      # The first node of a generation.lattice graph is a corner.        
        gFinal.add_edge(atNode, freeBBnodes[k]); k+=1
        gFinal.vp["attachments"].a[atNode] = 1
        gFinal.vp["attachments"].a[freeBBnodes[k]] = 1        
    assert (topology.label_components(gFinal)[1][0] == gFinal.num_vertices()) # gFinal must be Fully Connected.
    return gFinal


    
def palets(mapType):
    if mapType == "intToColor":
        return {0 :"b", 1: "r", 2:"g", 3:"y", 4:"m", 5:"c", 6:"k", 7:"burlywood", 8:"chartreuse", 9:"#FF00CC", 10:"#FF6600"}
    elif mapType == "typeToColor":
        return {"Gnp_Nodes" :"blue",  "Cliques": "red", "Circles":"green", "Lattices":"magenta"}
    elif mapType == "posToColor":
        return {"N" :"blue",  "V": "red", "J":"green"} #N=noun, V=verb


def colorNodes(g, nodeTypes):            
    colorMap = np.zeros(g.num_vertices())   #all non-specified nodes are coloro blue
    
    if "cliq" in nodeTypes:
        colorMap[np.where(g.vp["cliq"].a == 1)] = 1
    if "circ" in nodeTypes:
        colorMap[np.where(g.vp["circ"].a == 1)] = 2
    if "star" in nodeTypes:
        colorMap[np.where(g.vp["star"].a == 1)] = 3
    if "late" in nodeTypes:
        colorMap[np.where(g.vp["late"].a == 1)] = 4

    colorMap = colorMap.astype(np.int8)    
    colorIndex = ["blue", "red", "green", "yellow", "magenta"]
    colors = [colorIndex[x] for x in colorMap]
    return colors
    


def example_signaturesAndPCA():
#     Calculate signatures of random+clique+lattices+circles mixture graph and plot their PCA                                        
    crc     = erdos_circles_cliques([1000, 5000], nCliques=20, cliqueSize=15, nCircles=20, circleSize=30, save=False)
    gFinal  = append_lattices(crc, nLatte = 20, latticeSize=[5,5])        
    colors  = colorNodes(gFinal, ["cliq", "circ", "late"])
    Lg      = spectral.laplacian(gFinal)        

    eigsToCompute = 100
    eigsWhere ="../Data/Spectrums/LM/mix2"+str(eigsToCompute+1)+"_Eigs"    
    evals, evecs = graph_analysis.IO.compute_or_load(gk.computeSpectrum, eigsWhere, False, Lg, eigsToCompute+1, small=True)
    timeSample = gk.heatTimeSample3(evals[0], evals[-1], 10)    
    sig = gk.HKS(evals, evecs, eigsToCompute, timeSample)    

    sigsPCA = PCA(sig) 
    print "Fraction of Variance in the first three PCA dimensions: " +str(sum(sigsPCA.fracs[0:3]))
    saveAt = "gnp_cliques_circles__lattices_1K_5K_20_15_20_30_20_25.pdf"
    myPLT.plot_PCA_Res(sigsPCA, dim=3, legendLabel=["Gnp_Nodes", "Cliques", "Circles", "Lattices"], colorMap=[colors, palets("typeToColor")], save=False, saveAt=saveAt)
    

def make_crcl_graphs(nExamples, erdosSize, clic_circ_lat_N, clic_circ_lat_Dim, saveAtFolder=None, graphID='crcl'):
    nCliques, nCircles, nLatte          = clic_circ_lat_N
    cliqueSize, circleSize, latticeSize = clic_circ_lat_Dim
    allGraphs = []
    for i in xrange(nExamples):
        inGraph     = erdos_circles_cliques(erdosSize, nCliques, cliqueSize, nCircles, circleSize)
        inGraph     = append_lattices(inGraph, nLatte, latticeSize)    
        allGraphs.append(inGraph)        

        if saveAtFolder != None:
            graphName   = graphID + "_" + str(i)
            fileName    = saveAtFolder + "/" + graphName
            inGraph.save(fileName+".graph.xml")
            IO.from_GT_To_Greach(inGraph, fileName + ".greach.graph")
            IO.from_GT_To_Snap  (inGraph, fileName + ".snap.graph")


if __name__ == "__main__":                 
    saveAtFolder = "../../Data/Graphs/Synthetic/"
    make_crcl_graphs(2, erdosSize=(1000, 5000), clic_circ_lat_N=(10,10,10), clic_circ_lat_Dim=(20, 20, (4,5)), saveAtFolder=saveAtFolder, graphID="crcl_1k_5k_10ccl_20_nodes_each")
    sys.exit("Exiting Successfully.")