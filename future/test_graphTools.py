'''
Created on May 6, 2014

@author: optas
'''
import unittest, sys
from graph_analysis import IO
from graph_analysis.basic_tools import structEquivNodes, haveSameNeighbors, isValidStructEquivClass, egoNetwork
from synthetic_graphs.smallGraphs import starGraph
from unittest.suite import TestSuite

try:
    from graph_tool import generation, collection
except:
    print "Graph_Tool is not installed in this environment."
    sys.exit(1)



class Test(unittest.TestCase):
    
    def testSameNeighbors(self, n=10):
        '''Testing  graphTools.haveSameNeighbors() with a clique example.'''        
        g = generation.complete_graph(n, self_loops=False, directed=False)  
        for i in xrange(g.num_vertices()):
            for j in xrange(i+1, g.num_vertices()):
                self.assertTrue(haveSameNeighbors(g.vertex(i), g.vertex(j), g))
        g.remove_edge(g.edge(0,n-1))
        self.assertTrue(haveSameNeighbors(g.vertex(0), g.vertex(n-1), g))        
        for i in xrange(1, g.num_vertices()-1):            
            self.assertFalse(haveSameNeighbors(g.vertex(i), g.vertex(0), g))
            self.assertFalse(haveSameNeighbors(g.vertex(i), g.vertex(n-1), g))
        
    def test_CompleteGraph_equivClasses(self, n=10):
        ''' Testing structurally equivalent classes of a complete graph.'''        
        g = generation.complete_graph(n, self_loops=False, directed=False)        
        gClasses  = structEquivNodes(g)        
        for i in xrange(n):
            self.assertTrue(sorted(gClasses[i]) == range(n), "In a complete graph all nodes form 1 Equivalence Class")
        
        self.assertTrue(isValidStructEquivClass(gClasses, g))

    def test_StarGraph_equivClasses(self, n=10):
        ''' Testing structurally equivalent classes of a star graph'''        
        g         = starGraph(n)
        gClasses  = structEquivNodes(g)
        self.assertTrue(len(gClasses.keys()) == n-1)
        self.assertTrue(0 not in gClasses.keys())
        for i in xrange(1, n):
            self.assertTrue(sorted(gClasses[i]) == range(1,n))
        self.assertTrue(isValidStructEquivClass(gClasses, g))
    
    def test_examples_with_no_equiv_Classes(self, n = 100):
        g = generation.circular_graph(n)
        gClasses  = structEquivNodes(g)
        self.assertTrue(len(gClasses) == 0)
        
        g = generation.lattice([n,n])
        gClasses  = structEquivNodes(g)
        self.assertTrue(len(gClasses) == 0)
        
    def test_equivClassesForIsoGraphs(self, exampleG='cond-mat-2003'):
        ''' Do two isomorphic graphs have the same equivalence classes of nodes? '''        
        gcc = True; save = False;
        Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Matrices/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", save, exampleG, gcc)
        g                   = graph_analysis.IO.load_data("../Data/Graphs/"+exampleG+"_GCC_purged").next()
        gp                  = graph_analysis.IO.load_data("../Data/Graphs/"+exampleG+"_GCC_Iso").next()
                
        gClasses  = structEquivNodes(g)              
        self.assertTrue(isValidStructEquivClass(gClasses, g))
        
        gpClasses = structEquivNodes(gp)        
        self.assertTrue(isValidStructEquivClass(gpClasses, gp))
                
        self.assertEqual(len(gClasses.keys()), len(gpClasses.keys()))
        
        for node, equivs in gClasses.iteritems():            
            nodeIso   = forwardMap[node]
            equivsIso = gpClasses[nodeIso]
            self.assertEqual(len(equivs), len(equivsIso))
            for n in equivs:
                self.assertTrue(forwardMap[n] in equivsIso)

    def test_egoNetwork(self, graphExample = 'karate'):
        '''The egoNetwork of node -n- must at least have the same number of nodes as n's neighbors + 1 '''
        inGraph = collection.data[graphExample]
        if inGraph.is_directed:
            inGraph.set_directed(False)

        for node in inGraph.vertices():
            ego = egoNetwork(inGraph, node)            
            self.assertTrue(ego.num_vertices() == node.out_degree() + 1)
            

fast = TestSuite()
fast.addTest(Test.test_egoNetwork)

if __name__ == "__main__":
#     unittest.TestSuite(fast)    
    unittest.main()
    