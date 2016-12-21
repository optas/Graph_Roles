'''
Created on Jul 22, 2014

Contact: pachlioptas@gmail.com

Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import sys
from graph_analysis import IO
from graph_tool import topology, Graph
from myUtils import p
from presentation import *




def import_St_Mark_Data(save=True, export=True):
    
    saveLoadFolder = "StMarks"    
    graphName = "StMarks"
    graphFile = "../Data/Graphs/"+saveLoadFolder+"/StMarks_adjusted_raw.txt"

    g = Graph(directed=False)
    with open(graphFile, "r") as inF:
        num_nodes = int(inF.readline().split()[0])                
        g.add_vertex(num_nodes)        
        for line in inF:            
            if line[0] == "*": break
            node1, node2, unk = line.split()
            fromE, toE = int(node1)-1, int(node2)-1            
            g.add_edge(fromE, toE)

        gMass  = g.new_vertex_property("double") #Read-record Bio-mass
        for i, line in enumerate(inF):                   
            gMass[g.vertex(i)] = float(line.split()[0])
    
    trophicF = graphFile = "../Data/Graphs/"+saveLoadFolder+"/10_Isotrophic_Groups.txt"

    troClass  = g.new_vertex_property("int") #Read-record Trophic Classes

    with open(trophicF, "r") as inF:
        for i, line in enumerate(inF):
            classMembers = line.split(",")
            for m in classMembers:                
                troClass.a[int(m)-1] = i
#                 print m, troClass.a[int(m)-1]    
    
    g.vp["mass"] = gMass
    g.vp["trophicClass"] = troClass
    
    g = graph_analysis.IO.make_simple_graph(g, undirected=True, gcc=True)

    if save:        
        graph_analysis.IO.save_data("../Data/Graphs/"+ saveLoadFolder + "/" + graphName +".GT.graph", g)
    if export:
        exportToSnapAndGreach(graphName, saveLoadFolder)
    
    return g
        

    
def import_e_coli_ppi(save=False, export=False):
    '''
    Imports the dataset E_Coli and saves it as a graph (Snap, GTools and Greach format).
    '''
    saveLoadFolder = "E_Coli"
    graphName      = "E_Coli"
    graphFile      = "../Data/Graphs/"+saveLoadFolder+"/E_Coli_Edge_List.txt"
#     labelsFile     = "../Data/Graphs/"+saveLoadFolder+"/E_Coli_Essential_Genes.txt"
    labelsFile     = "../Data/Graphs/"+saveLoadFolder+"/E_Coli_Labels.csv"
    
    g = Graph(directed=False)
#     interactionWeight = g.new_edge_property("double")
    
    proteinNameToNode = dict()
    with open(graphFile, "r") as inF:
        for line in inF:
#             splitted = line.rstrip().split('|')
#             fromNode = splitted[1].strip()
#             toNode   = splitted[3].strip()
#             weight   = float(splitted[4])

            fromNode, toNode = line.strip().split()
            fromNode = fromNode.strip()
            toNode = toNode.strip()
            
#             print fromNode, toNode


#             print  fromNode, toNode, weight
            
            if fromNode not in proteinNameToNode:
                newNode = g.add_vertex()
                proteinNameToNode[fromNode] = int(newNode)
            if toNode not in proteinNameToNode:
                newNode = g.add_vertex()
                proteinNameToNode[toNode] = int(newNode)
            source = proteinNameToNode[fromNode]
            target = proteinNameToNode[toNode]
            edge = g.add_edge(g.vertex(source), g.vertex(target))
#             interactionWeight[edge] = weight

    essentiality = g.new_vertex_property("short")
    essentiality.a = 0
    symbolToInt = {'N':0, 'E':1, '?':'2', 'X':3}
    
    print g
    
    import csv 
    with open(labelsFile, "r") as inFile:
        count = 0 
        data = [row for row in csv.reader(inFile.read().splitlines())]
        for pair in data:
            proteinName, attribute = pair
            try:
                essentiality.a[proteinNameToNode[proteinName.lower()]] = symbolToInt[attribute]
            except:
                count +=1 
        print count 

    g.vp["essentiality"] = essentiality
        
    lethalOrNot = essentiality.a == 0
    lethalOrNot += essentiality.a == 1
    
    
    lethality    = g.new_vertex_property("boolean")
    lethality.a  = lethalOrNot

    
    g.set_vertex_filter(lethality)
    g.purge_vertices()
    print g
    p()

    lethality.a = 0
    lethality.a[essentiality.a==1] = 1
    
    g.vp["lethality"] = lethality

    
    
#     lethality    = g.new_vertex_property("boolean")
#     lethality.a  = 0
    
    
#     with open(labelsFile, "r") as inF:
#         for line in inF:
#             try:
#                 nodeID = proteinNameToNode[line.rstrip()]
#                 lethality.a[nodeID] = 1
#             except: #we don't have this node in the PPI net
#                 pass
#     
#     
#     g.vp["lethality"] = lethality
# #     g.ep["weights"]   = interactionWeight


        
    g = graph_analysis.IO.make_simple_graph(g, undirected=True, gcc=True)

    if save:        
        graph_analysis.IO.save_data("../Data/Graphs/"+ saveLoadFolder + "/" + graphName +".GT.graph", g)
    if export:
        exportToSnapAndGreach(graphName, saveLoadFolder)

    return g
            
            
def mapk_protein_network(save=False, export=False):
    saveLoadFolder = "Mapk_new"
    graphName      = "Mapk"
    graphFile      = "../Data_new/Graphs/"+saveLoadFolder+"/Protein_Edge_List.csv"
    labelsFile     = "../Data_new/Graphs/"+saveLoadFolder+"/Protein_Classes.txt"

    #read-load graph
    g                 = Graph(directed=False)
    proteinName       = g.new_vertex_property("string")
    proteinNameToNode = dict()
    nodesIn           = 0
    with open(graphFile, "r") as inF:
        for line in inF:
            if line[0] == "#": continue            
            node1ID, node2ID, node1Name, node2Name, garbage  = line.strip().split(",")
            if node1Name not in proteinNameToNode:
                newNode = g.add_vertex()
                proteinName[newNode] = node1Name
                proteinNameToNode[node1Name] = int(newNode)
                assert(proteinNameToNode[node1Name] == nodesIn)
                nodesIn +=1
            if node2Name not in proteinNameToNode:
                newNode = g.add_vertex()
                proteinName[newNode] = node2Name
                proteinNameToNode[node2Name] = int(newNode)
                assert(proteinNameToNode[node2Name] == nodesIn)
                nodesIn +=1
                
            source = proteinNameToNode[node1Name]
            target = proteinNameToNode[node2Name]
            g.add_edge(g.vertex(source), g.vertex(target))

    print g
    
    nodeClasses = g.new_vertex_property("short")
    nodeClasses.a = np.zeros_like(nodeClasses.a)    
    classI = 1
    #read-load labels    
    with open(labelsFile, "r") as inF:        
        for line in inF:
            if line[0] == "#": continue            
            labels  = line.strip().split("\t")
            for i in range(1, len(labels)):
                try:
                    nodeClasses[g.vertex(proteinNameToNode[labels[i]])] = classI                    
                except:
                    print "not found ", labels[i]
            classI += 1
                
    g.vp["proteinName"] = proteinName
    g.vp["Classes"]     = nodeClasses

    graph_analysis.IO.make_simple_graph(g, undirected=True, gcc=True)
        
    if save:        
        graph_analysis.IO.save_data("../Data_new/Graphs/"+ saveLoadFolder + "/" + graphName +".GT.graph", g)
    if export:
        exportToSnapAndGreach(graphName, saveLoadFolder)

    return g

    
# import graph_kernels as gk
import paper
if __name__ == "__main__":

    graphName = "Plato_Republic"
    inGraph, groupTaxa, blackList, xTickMarks = prepare_input_graph(graphName, metric="default")
    roles.plot_roles(inGraph, groupTaxa, graphName)
    
#     for graphName in paper.graphs:
#         inGraph, groupTaxa, blackList, xTickMarks = prepare_input_graph(graphName, metric="default")
# #         
#         storedFolder                              = roles.graphFolder(graphName)
#         inGraph.save("../Data/Graphs/"+storedFolder+"/"+graphName + ".graph.xml")
#  
    sys.exit("Exiting Successfully.")