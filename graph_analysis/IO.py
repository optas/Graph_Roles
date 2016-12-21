'''
Created on Mar 24, 2014

Contact: pachlioptas@gmail.com

Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import random, math, itertools, cPickle, os
import numpy as np
from string import rstrip

try:
    from graph_tool import collection, topology, spectral
    from graph_tool import stats as gtStats
    from graph_tool import *
except:
    "Graph_Tool is not installed in this environment."

    
def save_data(fileName, *args):
    myFile = open(fileName, "w")
    cPickle.dump(len(args), myFile)
    for item in args:
        cPickle.dump(item, myFile)
    myFile.close()    

def load_data(fileName):    
    inFile = open(fileName, "r")    
    size = cPickle.load(inFile)            
    for i in xrange(size):
        yield cPickle.load(inFile)        
    inFile.close()

def load_GT_graph(graphExample, gcc=False, removeSL=False, removePE=False):
    '''
    Input:  - graphExample,  graph example from Graph-tool collections  (e.g.,  'cond-mat-2003', 'adjnoun' 'karate' 'netscience') or a graphfile in .gml format
            - gcc = True if only the giant connected component should be returned
            - removeSL = True if any self-loops must be removed 
            - removePE = True if any parallel-edge must be removed
    Output: the corresponding graph_tool graph object
    '''

    if graphExample[-4:] == ".gml":
        g = load_graph(graphExample)
    else:
        g = collection.data[graphExample]
        
    if g.is_directed:
        g.set_directed(False)
#         g = Graph(g, directed=False)
    
    if removePE:
        gtStats.remove_parallel_edges(g)    
    if removeSL:        
        gtStats.remove_self_loops(g) 
    if gcc:    
        l = topology.label_largest_component(g)  #Keep Largest Connected Component
        g.set_vertex_filter(l)
        #g = GraphView(g, vfilt=l)
    g.purge_vertices()
    
    return g

def compute_or_load(computeFunction, fileName, save=True, *args, **kargs):
    '''
    It returns the objects saved at fileName or the fileName does not correspond to an actual file, it executes the computeFunction and returns the 
    nodeClusteringResults.
    Precondition: if fileName exists, then should corresponds to a pickle where the first thing that is saved is the number of the pickled objects.
    Safe to use for a fileName that was created with saceData
          
    '''
    if os.path.isfile(fileName):
        print "Loading " + fileName
            
        inFile = open(fileName, "r")
        size   = cPickle.load(inFile)
        
        if size == 1:
            load = cPickle.load(inFile)
            inFile.close()
            return load
        else:
            inFile.close()
            return load_data(fileName)
    else:               
        print "Computing " + computeFunction.func_name
        if kargs == None:
            res = computeFunction(*args)
        else:
            res = computeFunction(*args, **kargs)
        
        if save:
            myFile = open(fileName, "w")
            if isinstance(res, tuple):    # I.e., the function returned an actual tuple, or many objects in a tuple.
                cPickle.dump(len(res), myFile)
                for item in res:
                    cPickle.dump(item, myFile)
            else:                         # The function returned a single object.
                cPickle.dump(1, myFile)
                cPickle.dump(res, myFile)
            myFile.close()               
        return res
        
def makeFolder(folderPath):
    if not os.path.exists(folderPath):
        os.makedirs(folderPath)
        
def rearrangeRowCols(A, perm):
    '''
    Input:
        A: Symmetric matrix
        perm: vector of permuted indices
    Output:
        Ap: A permuted matrix of -A- according to the order of -perm- 
    '''
    
    n = A.shape[0]
    assert(len(perm) == n)
    Ap = np.zeros_like(A)
 
    for x in xrange(n):
        for y in xrange(x,n): # Only need to traverse half the matrix since it is symmetric
            Ap[x, y] = Ap[y, x] = A[perm[x], perm[y]]    
    return Ap
                   

def mcSherryFromGraph(inGraph):    
    lng = spectral.laplacian(inGraph, normalized=True)
    lng = -lng
    lng.setdiag(np.zeros(lng.diagonal().shape))
    return lng


def mcSherryFromLaplacian(Lg):
    assert(Lg.shape[0]==Lg.shape[1])
    Mg = Lg.copy()    
    
    for i in xrange(Mg.shape[0]):
        for j in xrange(i+1, Mg.shape[0]):
            if Mg[i,j] == -1: 
                newV = 1 / math.sqrt(Mg[i,i]*Mg[j,j])
                Mg[i,j] = Mg[j,i]= newV
                     
    Mg.setdiag(np.ones(Mg.diagonal().shape))
    return Mg
    
      
def adjacencyToGraph(A):
    '''
    Input: A, a symmetric adjacency matrix of an undirected graph 
    Output: Graph corresponding to -A-
    Dependencies: graph-tool
    '''
    n = A.shape[0]
    assert n == A.shape[1]     
    g = Graph(directed=False)
    g.add_vertex(n)
    
    for i in xrange(n):
        for j in xrange(i,n):
            if A[i,j] == 1: g.add_edge(g.vertex(i), g.vertex(j))
            else:
                assert (A[i,j] == 0)         
    return g
        
         
          
def isBinary(A):
    for i in xrange(A.shape[0]):
        for j in xrange(A.shape[1]):
            if A[i,j] not in [0,1]: 
                return False
    return True 



    

def createIsoGraph(adjacencyM, assertIso=False, inGraph = None):
    '''
    Input:   - adjacencyM = Adjacency matrix of graph for which we will create an isomorphic one
             - assertIso = if True assert the isomorphy between the produced graph and the -inGraph- (very_slow, use only for small graphs)
             - inGraph = graph corresponding to -adjacencyM- used only for asserting isomorphy
                          
    Output: 1) An Isomorphic graph to inGraph, 
            2) The mapping from vertices of the 'example' graph to its isomorphic one.  
    '''
        
    perm       = range(adjacencyM.shape[0]); random.shuffle(perm)       #a permutation vector over the rows (integers) of the Laplacian            
    forwardMap = [i[0] for i in sorted(enumerate(perm), key=lambda x:x[1])]  #which node corresponds to which 
    Ap         = rearrangeRowCols(adjacencyM, perm)            
    gp         = adjacencyToGraph(Ap)
    if assertIso:
        assert(inGraph != None)
        assert(topology.isomorphism(inGraph, gp))
    return gp, forwardMap


def laplacianToGraph(L):
    '''
    Input: L, a symmetric Laplacian matrix of an undirected graph 
    Output: Graph corresponding to -L-
    Dependencies: graph-tool
    '''
    n = L.shape[0]
    assert n == L.shape[1]     
    g = Graph(directed=False)
    g.add_vertex(n)
    
    for i in xrange(n):
        for j in xrange(i+1,n):
            if L[i,j] == -1: 
                g.add_edge(g.vertex(i), g.vertex(j))     
    return g
        

                                  
def makeISOLaplacians(graphExample, gcc, save=False, saveAt=None):
    
    g  = load_GT_graph(graphExample, gcc, removeSL=True, removePE=True) #For forming the Laplacian, we have to remove the Self-Loops and Parallel Edges if any.
    Ag = spectral.adjacency(g).todense()    #Despite this todense() the Laplacian is going to be sparse. It is though taking more time to create it.  
    gp, forwardMap = createIsoGraph(Ag)
     
    Lg  = spectral.laplacian(g)  #create the Laplacians  
    Lgp = spectral.laplacian(gp)
    
    if save:
        assert(saveAt != None)
        save_data(saveAt, Lg, Lgp, forwardMap)
    
    return Lg, Lgp, forwardMap
 
 
def from_GT_To_Snap(inGraph, outFile):
    '''
    Takes a graph in graphTool's format and produces a file that encodes the graph in Snap readable format.
    '''
    totalNodes = inGraph.num_vertices()
    totalEdges = inGraph.num_edges()
    directed   = inGraph.is_directed()
    
    with open(outFile, "w") as fOut:
        
        fOut.write("#Directed = " + str(directed) + "\n")
        fOut.write("#Nodes = " + str(totalNodes) + "\n")
        fOut.write("#Edges = " + str(totalEdges) + "\n")
        
        for edge in inGraph.edges():        #TODO check if Snap relies on the order the edges are listed?
            fOut.write(str(edge.source()) +" "+ str(edge.target()) + "\n")


def from_GT_To_Greach(inGraph, outFile):
    '''
    Takes a graph in graphTool's format and produces a file that can be used by greach (presentation).
    '''
    with open(outFile, "w") as fOut:
        fOut.write("graph_for_greach\n")
        fOut.write(str(inGraph.num_vertices())+"\n")
        i = 0
        for v in inGraph.vertices():
            assert(int(v)==i)
            fOut.write(str(i)+": ")
            neighbors = [int(n_v) for n_v in v.out_neighbours()]
            neighbors.sort()
            for n in neighbors:
                fOut.write(str(n)+" ")
            fOut.write("#\n")
            i += 1
       


def readRoleSimMatrix(inFile):
    '''
    Created to read the output of presentation (i.e., the square matrix with presentation or simRank scores)
    '''
    inMatrix = list()
    with open(inFile, "r") as inF:
        for line in inF:
            row = list()
            for value in line.split():
                row.append(float(value))
            inMatrix.append(row)
    return np.array(inMatrix)



def loadGraphWithAnnotations(graphFile):
    '''
    Used to read the graphs provides by RoleSim people, regarding scientific collaborations
    and the g/h index.
    '''
    g = Graph(directed=False)
    with open(graphFile, "r") as inF:
        num_nodes = int(inF.readline().split()[1])
        g.add_vertex(num_nodes)
        g_names  = g.new_vertex_property("string")
        g_H_ind  = g.new_vertex_property("int")
        g_G_ind  = g.new_vertex_property("int")
        
        for i, line in enumerate(inF):                        # Read Meta-Data of Nodes
            if rstrip(line) == "*Edges": break
            contents = rstrip(line).split("\t")
            gID, name, gIndex, hIndex = contents[0], contents[1], int(contents[2]), int(contents[3])        
            assert(gID == str(i))            
#             print gID, name, gIndex, hIndex
            g_names[g.vertex(i)] = name
            g_H_ind.a[i] = gIndex
            g_G_ind.a[i] = hIndex
        
            
        for i, line in enumerate(inF):                        # Read Edges
            tokens = line.split()
            fromE, toE = int(tokens[0]), int(tokens[1])            
            g.add_edge(fromE, toE)
        
        g.vp["names"] = g_names
        g.vp["h-Index"] = g_H_ind
        g.vp["g-Index"] = g_G_ind

    gtStats.remove_parallel_edges(g)
    gtStats.remove_self_loops(g)
    l = topology.label_largest_component(g)  #Keep Largest Connected Component    
    g.set_vertex_filter(l)
    g.purge_vertices()
    
    return g


def saveLaplacianToMatlab(graph, normal=False, saveAt=None):
    '''
    Save the laplacian of the -graph- in the sparse format that is used by Matlab.
    '''
    nodes = graph.num_vertices()
    outFormat = "%d %d %d\n"
    with open(saveAt, "w") as outF:
        outF.write(outFormat % (nodes, nodes, 0)) #specify size of full matrix
        cx = spectral.laplacian(graph, normalized=normal).tocoo()        
        for i,j,v in itertools.izip(cx.row, cx.col, cx.data):
            outF.write(outFormat % (i+1, j+1, v))               #Matlab indexes arrays from 1. 


    


#     print "IO module: Running main()" 
#     import myBlockModels as myBM
#     N = 1000                    #Keep Note
#     imageMatrix, edgeDisMatrix       = myBM.manufacturingImage()
#     roleDistribution                 = [0.10, 0.15, 0.15, 0.25, 0.10, 0.25]    
#     graphName                        = "manufacturing_2" 
#     inGraph                          = myBM.buildModel(N, imageMatrix, roleDistribution, edgeDisMatrix)
#     inGraph                          = myBM.applyManufacturingConstraints(inGraph)
#     
#     save_data("../Data/Graphs/synthetic/" + graphName + ".GT.graph", inGraph)
# 
#     from_GT_To_Greach(inGraph, "../Data/Graphs/synthetic/" + graphName + ".greach.graph")
#     from_GT_To_Snap(  inGraph, "../Data/Graphs/synthetic/" + graphName + ".snap.graph")
                                                                

# if __name__ == "__main__":        
# 
#     graph_name = "serengeti-foodweb"    
#     save_at    =  graph_name + "_adjacency_matlab"
# 
#     in_graph, group_taxa, blackList, x_tick_marks = mc.prepare_input_graph(graph_name, metric="default")
#     print blackList
# 
#     nodes = in_graph.num_vertices()    
#     edges = in_graph.num_edges()
#     output_format = "%d %d %d\n"
#     
#     with open(save_at, "w") as outF:
#         outF.write(output_format % (nodes, 0, 0)) # Specify number of nodes.
#         outF.write(output_format % (edges, 0, 0)) # Specify number of edges.
#         cx    = spectral.adjacency(in_graph).tocoo()
#         for i,j,v in itertools.izip(cx.row, cx.col, cx.data):
#             outF.write(output_format % (i+1, j+1, v))         # Matlab indexes arrays from 1.
#     
# 
#     sys.exit("Normal Exiting")