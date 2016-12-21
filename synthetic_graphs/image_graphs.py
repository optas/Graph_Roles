'''
Created on Jun 25, 2014

@author: optas
'''

import sys
from graph_tool.all import *
from myUtils import *
from scipy import spatial as spt
from role_discovery import graph_kernels as gk
from synthetic_graphs import smallGraphs as sg
from graph_analysis import IO
from scipy.stats import pearsonr


def manufacturing_image_graph(probabilistic=None):
    '''
    Return the image graph encoded as a matrix for the manufacturing model.
    It is also generating the edge densities.
    '''    
    r = 6 # Total number of roles in a manufacturing graph.
    imageMatrix = np.zeros((r,r))    
    imageMatrix[0,1] = 1; 
    imageMatrix[1,2] = 1; 
    imageMatrix[2,3] = imageMatrix[2,4] = 1; 
    imageMatrix[3,4] = 1;
    imageMatrix[3,5] = 1;    
    imageMatrix = imageMatrix + np.transpose(imageMatrix) # The model is undirected that the image graph is undirected.
    
    edgeDisMatrix = dict()
    edgeDisMatrix[(0,1)] = (4,4)
    edgeDisMatrix[(1,2)] = (3,5)
    edgeDisMatrix[(2,3)] = (6,6)    
    edgeDisMatrix[(2,4)] = (2,2)
    edgeDisMatrix[(3,4)] = (3,3)
    edgeDisMatrix[(3,5)] = (2,10) 
     
    return imageMatrix, edgeDisMatrix

    
def assign_node_roles(graph, roleDistribution):
    assert (abs(sum(roleDistribution) - 1) < 1e-3)    # roleDistribution is a probability measure.

    N       = graph.num_vertices()
    roles   = graph.new_vertex_property("int")    
    num_roles   = len(roleDistribution)
    roleGroupID = dict()
    lastLabeled = 0
    for i in range(num_roles-1):        
        num_i = int(np.floor(roleDistribution[i] * N))                       # Number of nodes with role -i-.
        
        roles.a[lastLabeled:lastLabeled+num_i] = i
        roleGroupID[i] = (lastLabeled, lastLabeled+num_i)        
        lastLabeled = lastLabeled+num_i

    roles.a[lastLabeled:] = i+1                                              # The remaining nodes get the last role.
    roleGroupID[i+1]  = (lastLabeled, N)
    graph.vp['roles'] = roles
    return roleGroupID
    
    
    
    
def orderAdmissibleRolePairs(numRoles, imageGraph):
    res = []
    for i in range(numRoles): #iterate over all roles
        admissible = np.where(imageGraph[i,:] == 1)[0]
        admissible = admissible[admissible >= i]
        for j in admissible:
            res.append((i,j))
    return res
    



def buildModel(N, imageGraph, roleDistribution, edgeDistribution):
    m, n = imageGraph.shape
    assert (m == n)     # Square_matrix
    assert (len(roleDistribution) == m)
    
    G            = Graph(directed=False)
    G.add_vertex(N)        
    roleGroupID  = assign_node_roles(G, roleDistribution)
    rolePairs    = orderAdmissibleRolePairs(m, imageGraph)
    
    print rolePairs
    
    
    for i, j in rolePairs:      #Role-i CAN connect to role-j             

        nodesTypeI = range(roleGroupID[i][0], roleGroupID[i][1])
        nodesTypeJ = range(roleGroupID[j][0], roleGroupID[j][1])
        minDeg, maxDeg = edgeDistribution[(i,j)]

        for nodeI in nodesTypeI:
            degree = np.random.randint(minDeg, maxDeg+1)
            np.random.shuffle(nodesTypeJ)
            connectTo = nodesTypeJ[:degree]
#             print nodeI, connectTo
            for nodeJ in connectTo:
                G.add_edge(nodeI, nodeJ)
        
    return G



def applyManufacturingConstraints(graph):
    blocks = graph.vp["blockIDs"]
    blackList = list()
#     blackListEdges = graph.new_edge_property("short")
#     blackListEdges.a = np.ones(len(blackListEdges.a))
    
    for role in [1,2,4,3,5]:                              #we must visit 4 before 3 (a shop can get products either from delivery or from warehouse)
        for v in np.where(blocks.a == role)[0]:           #nodes of particular role (in increasing order)
            neighborTypes = []
            for neighb in graph.vertex(v).out_neighbours():
                if neighb not in blackList:
                    neighborTypes.append(blocks[neighb])
            
            if (role == 1 and 0 not in neighborTypes) \
            or (role == 2 and 1 not in neighborTypes) \
            or (role == 4 and 2 not in neighborTypes) \
            or (role == 3 and 2 not in neighborTypes and 4 not in neighborTypes) \
            or (role == 5 and 3 not in neighborTypes):
                blackList.append(int(v))
#                 for e in graph.vertex(v).out_edges():
#                     print type(e)
#                     blackListEdges[e] = 0
    

    bl = graph.new_vertex_property("bool")
    bl.a = np.ones(len(bl.a))        
    bl.a[blackList] = 0    
    print bl.a
    print graph
    graph.set_vertex_filter(bl)
    print graph.num_vertices()
#     graph = GraphView(graph, vfilt=bl.a)
    graph.purge_vertices(in_place=False)
    print graph
    return graph
    



def example_blockModel():
    '''
    Traditional assortative block model, quick and dirty example.
    '''
    N            = 100   #number of Nodes
    blocksNum    = 5
    averageEdges = 6

    def corr(a, b):
        if a == b:          #inner-block Probability of connection
            return 0.9
        else:
            return 0.001    #cross-block Probability of connection

    g, bm = random_graph(N, lambda: np.random.poisson(averageEdges), directed=False, model="blockmodel-traditional",\
                            block_membership = lambda: np.random.randint(blocksNum), vertex_corr=corr)
    print bm
    
    for i in xrange(blocksNum):
        print "Block %d has %d nodes." % (i, len(np.where(bm.a == i)[0]))

    print bm
    print type(bm)

    graph_draw(g, vertex_fill_color=bm, edge_color="black", output="blockmodel.pdf")        
    return g        

    
            

if __name__ == "__main__":        


#     N = 1000        
#     imageMatrix, edgeDisMatrix       = manufacturing_image_graph()
#     roleDistribution                 = [0.15, 0.10, 0.15, 0.25, 0.10, 0.25]    
#     graphName                        = "manufacturing_1" 
#     G                                = buildModel(N, imageMatrix, roleDistribution, edgeDisMatrix)
#     G                                = applyManufacturingConstraints(G)


#     N = 1000        
#     imageMatrix, edgeDisMatrix       = manufacturing_image_graph()
#     roleDistribution                 = [0.10, 0.15, 0.15, 0.25, 0.10, 0.25]    
#     graphName                        = "manufacturing_2" 
#     inGraph                          = buildModel(N, imageMatrix, roleDistribution, edgeDisMatrix)
#     inGraph                          = applyManufacturingConstraints(inGraph)


    N = 1000        
    imageMatrix, edgeDisMatrix       = manufacturing_image_graph()
    roleDistribution                 = [0.10, 0.15, 0.15, 0.25, 0.10, 0.25]    
    graphName                        = "manufacturing_3" 
    inGraph                          = buildModel(N, imageMatrix, roleDistribution, edgeDisMatrix)
    print inGraph
#     inGraph                          = applyManufacturingConstraints(inGraph)

    sys.exit("Exiting Successfully.")
    
