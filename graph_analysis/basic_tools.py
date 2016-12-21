'''
Created on May 1, 2014

@author: optas
'''

try:
    from graph_tool import topology, GraphView
    from graph_tool import stats as gtStats
except:
    print "graph_tool library is not installed in this environment."

import numpy as np
from scipy import sparse
import copy
import myUtils as myUt
from collections import defaultdict


def make_simple_graph(g, undirected=True, gcc=True):
    '''
    Returns input graph -g- in a version without parallel edges or self-loops.
    If undirected = True, returned graph is also undirected.
    If gcc        = True, returned graph is giant connected component of g.
    '''

    if undirected and g.is_directed:
        g.set_directed(False)
        
    gtStats.remove_self_loops(g)        
    gtStats.remove_parallel_edges(g)
    
    if gcc:    
        l = topology.label_largest_component(g)  # Keep Largest Connected Component.
        print "Nodes in largest connected component: " + str(np.sum(l.a))
        g.set_vertex_filter(l)
    g.purge_vertices()
    return g
    
    

def addSymmetricEdge(g, fromN, toN):
    g.add_edge(fromN, toN)
    g.add_edge(toN, fromN)


def areIsoLaplacians(Lap1, Lap2, from1to2Map, verbose=True):    
    degrees1 = Lap1.diagonal()        
    degrees2 = Lap2.diagonal()
    
    isodegrees = [degrees2[iso] for iso in from1to2Map]    

    if np.any(degrees1 != isodegrees):
        if verbose:
            print "Laplacians do not match in the diagonals."
        return False
    
    if verbose:
        print "Diagonals are correct."
    
    L1  = Lap1.todense()
    L2  = Lap2.todense()
    
    rows = L1.shape[0]
    cols = L1.shape[1]    
    assert (rows == cols == L2.shape[0] == L2.shape[1] == len(from1to2Map))
    
    for i in xrange(rows):
        neihgb = []
        for j in xrange(cols):
            if L1[i,j] == -1:
                neihgb.append(j)
        
        isoNeighb = [from1to2Map[j] for j in neihgb]
        isoNeighb.sort()                 

        whatIhave = []        
        for k in xrange(cols):
            if L2[from1to2Map[i], k] == -1:
                whatIhave.append(k)

        whatIhave.sort()
        
        if whatIhave != isoNeighb:
            return False
        
    return True
        

def verifyIsomorphy(graph1, graph2, map1to2, verbose=True):
    '''TODO: implement for directed graph'''

    if graph1.num_edges() != graph2.num_edges():
        return False
    
    if graph1.num_vertices() != graph2.num_vertices():
        return False

    for e in graph1.edges():
        if graph2.edge(map1to2[int(e.source())], map1to2[int(e.target())]) == None :
            return False
    return True

def haveSameNeighbors(v1, v2, g):
    '''
    True if v1 has exactly same neighbors with v2 in g (excluding a potential connection with each other)
    Precondition: g is not directed
    '''
    if v1.out_degree() != v2.out_degree():
        return False    
    for v1n in v1.all_neighbours():
        if v1n == v2: continue
        if g.edge(v1n, v2) == None:
            return False
    return True


def isValidStructEquivClass(eqC, g):
    '''
    Checks if the proposed structurally equivalent class, is a valid one.
    Provides a fast way for ensuring that only false negatives might have occurred when creating the class.
    '''
    
    for node1, equivNode in eqC.iteritems():
        for node2 in equivNode:            
            if node1 == node2: continue
            if not haveSameNeighbors(g.vertex(node1), g.vertex(node2), g): return False    
    return True

def checkEquivNodesStandBetweenIso(g1classes, g2classes, forwardMap):
              
    if len(g1classes.keys()) != len(g2classes.keys()):
        return False
        
    for node, equivs in g1classes.iteritems():            
            nodeIso   = forwardMap[node]
            equivsIso = g2classes[nodeIso]
            if len(equivs) != len(equivsIso):
                return False
            for n in equivs:
                if forwardMap[n] not in equivsIso:
                    return False
    return True


def structEquivNodesSLOW(g, asList=False):
    
    equiV1 = defaultdict(list)
    done = set()
    for vertex1 in g.vertices():        
        if vertex1 in done: continue
        for vertex2 in vertex1.all_neighbours():   #neighbors of vertex1
            if int(vertex2) not in done and haveSameNeighbors(vertex1, vertex2, g):                
                equiV1[int(vertex1)].append(int(vertex2)); done.add(int(vertex2))
                 
            for vertex3 in vertex2.all_neighbours():
                if vertex3 == vertex1 or int(vertex3) in done: continue                
                if haveSameNeighbors(vertex1, vertex3, g):
                    equiV1[int(vertex1)].append(int(vertex3)); done.add(int(vertex3))
        done.add(int(vertex1))
    
    return myUt.propagateLabels(equiV1)
   
def structEquivNodes(g):
    equivClasses = list()
    ecChecked    = set()

    for vertex1 in g.vertices():    #test_myRolx equivalence BETWEEN all the nodes that are neighbors of vertex1                 
        v1NeighbSq = dict()                    
        for vertex2 in vertex1.all_neighbours():            
            if int(vertex2) in ecChecked: continue
            neighV2 = [int(n) for n in vertex2.all_neighbours()]
            neighV2.sort()
            try:            
                v1NeighbSq[tuple(neighV2)].append(int(vertex2))
            except:
                v1NeighbSq[tuple(neighV2)] = [int(vertex2)]
        
        for eq in v1NeighbSq.itervalues():
            if len(eq) > 1:
                equivClasses.append(eq)
                for node in eq:
                    ecChecked.add(node)
            ecChecked.add(eq[0])
                    
#     equivClasses = myUt.listGroupingToDict(equivClasses)    
    
    ecChecked = set()
    adjacentEquiv = dict()
    for vertex1 in g.vertices():    #test_myRolx equivalence between Vertex1 and its neighbors
        if int(vertex1) in ecChecked: continue
        adjacentEquiv[int(vertex1)] = [int(vertex1)]
        for vertex2 in vertex1.all_neighbours():
            if int(vertex2) in ecChecked: continue
            if haveSameNeighbors(vertex1, vertex2, g):
                ecChecked.add(int(vertex2))                
                adjacentEquiv[int(vertex1)].append(int(vertex2))
        if len(adjacentEquiv[int(vertex1)]) == 1: #doesn't have any "adjacent" equiv. node
            del adjacentEquiv[int(vertex1)]
            

    updatedEC = copy.deepcopy(equivClasses)
    for adjEquivs in adjacentEquiv.itervalues():
        merged = False;
        for node in adjEquivs:
            if merged: break;
            for i, oldEc in enumerate(equivClasses):
                if node in oldEc:
                    updatedEC[i] = updatedEC[i] + adjEquivs
                    merged = True; break;
        #didn't got merged
        updatedEC.append(adjEquivs)
         
    updatedEC = myUt.listGroupingToDict(updatedEC)
            
    return updatedEC
                             
def equivellantRows(A):  
    '''
    An O(nodes^2) way for finding nodes belonging at the same equivalence classes in the graph. - has not being tested - 
    '''
    
    B = A.todense()
    
    equiv = dict()
    seen = set()
    for i in xrange(B.shape[0]):
        if i in seen: continue
        seen.add(i)
        for j in xrange(i+1, B.shape[0]):
            if np.any(B[i,:] != B[j,:]) :
                continue 
            else:
                seen.add(j)
                try:
                    equiv[i].append(j)
                except:
                    equiv[i] = [j]

    return equiv




def swapEdges(g, noiseLevel):
    assert(noiseLevel>0 and noiseLevel<1)
        
    gnoisy = g.copy()
        
    e = gnoisy.num_edges()
    myBuffer = 2 #
    toSwap = np.ceil(e * noiseLevel)
    
    randomIndex = set(np.random.choice(e, size=toSwap*myBuffer, replace=False))
        
    edgesToSwap = []
    for i, edge in enumerate(gnoisy.edges()):
        if i in randomIndex:
            edgesToSwap.append(edge)
            
    swapped = 0    
    for i in xrange(0, len(edgesToSwap)-1, 2):        
        edge1 = edgesToSwap[i]
        edge2 = edgesToSwap[i+1]
        a = edge1.source(); b = edge1.target(); c = edge2.source(); d = edge2.target(); 
        if gnoisy.edge(a,c) == None and gnoisy.edge(b,d) == None:            
            gnoisy.add_edge(a,c)
            gnoisy.add_edge(b,d)
            gnoisy.remove_edge(edge1)
            gnoisy.remove_edge(edge2)
            swapped = swapped + 2
        if swapped >= toSwap:
            break
        
    if swapped < toSwap:
        print ("Warning, wanted to swap %d, but only $d were actually swapped."(toSwap, swapped)) 
        
    return gnoisy
    

def addUniformRandomEdges(laplacianG, noiseLevel, verbose=True):
    '''
    Modifies the Laplacian matrix of a symmetric graph in place, by adding uniformly random edges between the nodes.
    -noiseLevel- percentage over the total original edges to be added    
    '''
    assert (isinstance(laplacianG, sparse.lil.lil_matrix))
    
    n = laplacianG.shape[0]
    e = sum(laplacianG.diagonal())/2
    toAdd = np.ceil(e * noiseLevel)
    
    if verbose:
        print "Inside AddingUniformEdges: Input graph has %d nodes %d edges. Adding %d more edges" % (n, e, toAdd) 
            
    while True:        
        fromN = np.random.randint(low=0, high=n, size=toAdd)
        toN   = np.random.randint(low=0, high=n, size=toAdd)
        
        for i in xrange(len(fromN)):
            if laplacianG[fromN[i], toN[i]] == 0 and fromN[i] != toN[i]:
                
                laplacianG[fromN[i], toN[i]]  = -1                                
                laplacianG[fromN[i],fromN[i]] += 1
                laplacianG[toN[i], fromN[i]]  = -1  #symmetric graphs 
                laplacianG[toN[i],toN[i]] += 1
                toAdd -= 1
                
        if toAdd == 0: break
            
    return laplacianG.tocsr(), None          #You could leave in the discretion of the caller
    



def egoNetwork(inGraph, node):
    '''
    Compute the ego-network subgraph of the -inGraph- where the ego is the -node-.    
    Precondition: inGraph is undirected
    '''
    neighbors    = [int(n) for n in node.out_neighbours()]
    neighborhood = inGraph.new_vertex_property("bool")
    neighborhood.a[neighbors] = 1
    neighborhood.a[int(node)] = 1
    return GraphView(inGraph, vfilt = neighborhood) 
    
