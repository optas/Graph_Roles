import graph_tool as gt
import random
import numpy as np
from graph_analysis import IO
    
def sample_with_replacement(population, k):
    '''
    Chooses k random elements (with replacement) from a population.
    Input:
            population  -   a non empty list carrying objects corresponding to elements of the
                            population to be sampled.
            k           -   (positive int) the number of samples to be taken.
    Output:
                        -   a list with the sampled elements.

    '''
    n = len(population)
    _random, _int = random.random, int  # Speed hack.
    return [population[_int(_random() * n)] for i in xrange(k)]

def stratified_network_exponential_preference(nodes, classes, p0=0.12, alpha=2.0):
    '''
    Creates an undirected graph with stratified nodes. The probability of an edge between two nodes, drops exponentially
    in the distance between their corresponding classes.

    Input:
            nodes   -   (positive int), the number of nodes in the final graph.
            classes -   (list of positive int), classes[i] is the variable of stratification corresponding to
                        class i. The probability of an edge between a node belonging in classes[i] and a node
                        belonging in classes[j] is:
                                    p0 * exp(- alpha * (classes[i]-classes[j]) ).
            p0      -   (float, default = 0.2), parameter of the aforementioned probability.
            alpha   -   (float, default = 2.0), parameter of the aforementioned probability.

    Output:
                    -   (graph_tool graph) with specified number of nodes. Each node is associated with
                        a single class (stratification variable), encoded in a vertex_property_map named
                        "classes".

    ---------------------------------------------------------------------------------------------------------
    Example of Usage:
        g = stratified_network_exponential_preference(1000, range(10), p0=0.12, alpha=2.0)

        g now has 1000 nodes each belonging in one of 10 classes. Each of the ten classes correspond to a
        variable that takes one of 10 integer values from 1 to 10 (e.g., can be age of a child).

        In this example the probability of an edge drops by a factor of exp(1)^2 for every additional
        unit separating the classes of the nodes.

    '''

    def corr(a, b):
        return p0 * np.exp(-alpha * abs(a-b))
    
    classLabels = sample_with_replacement(classes, nodes)
    
    g = gt.Graph(directed=False)
    g.add_vertex(nodes)
        
    for nodeA in xrange(nodes):
        for nodeB in xrange(nodeA+1, nodes): 
            pAB = corr(classLabels[nodeA], classLabels[nodeB])
            if np.random.uniform() <= pAB:
                g.add_edge(nodeA, nodeB)    
    
    g.vp["classes"]    = g.new_vertex_property('int')
    g.vp["classes"].a  = np.array(classLabels)
            
    return g

from role_discovery import graphs_with_roles as roles

def make_stratified_graph(nodes, classes, graphName):    
    '''
    Paper oriented, make a stratified graph and save/export it in all relevant formats.
    '''
    inGraph         = stratified_network_exponential_preference(nodes, classes)    
    saveAtFolder    = roles.graph_folder(graphName)
    
    inGraph.save("../../Data/Graphs/"+saveAtFolder+"/"+graphName + ".graph.xml")
    graph_analysis.IO.from_GT_To_Greach(inGraph, "../../Data/Graphs/" + saveAtFolder + "/" + graphName + ".greach.graph")
    graph_analysis.IO.from_GT_To_Snap(  inGraph, "../../Data/Graphs/" + saveAtFolder + "/" + graphName + ".snap.graph")
    return inGraph


import sys

if __name__ == "__main__":                         
    make_stratified_graph(1000, range(1,11), 'stratified_1')
    sys.exit("Exiting Successfully.")
