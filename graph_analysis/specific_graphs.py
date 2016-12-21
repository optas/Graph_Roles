'''
Created on Feb 20, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''
from graph_tool import Graph, topology
from IO import save_data, from_GT_To_Greach, from_GT_To_Snap, load_GT_graph
from graph_analysis.basic_tools import make_simple_graph

saveLoadPath   = "../../Data/Graphs/"
    
def embed_k_core(graphName):
    inGraph              = load_GT_graph(graphName, gcc=True, removeSL=True, removePE=True)
    kcore                = topology.kcore_decomposition(inGraph, deg="total")
    inGraph.vp["kcore"]  = kcore
    return inGraph


def import_adj_noun(save=False):
    graphName       = "adjnoun"
    saveLoadFolder  = "GT_Coll"
    g               = load_GT_graph(graphName, gcc=True, removeSL=True, removePE=True)    
    if save:
        save_data(saveLoadPath + saveLoadFolder + "/" + graphName + ".GT.graph", g)
        g.save   (saveLoadPath + saveLoadFolder + "/" + graphName + ".graph.xml", fmt="xml")
        from_GT_To_Greach(g, saveLoadPath + saveLoadFolder + "/" + graphName + ".greach.graph")
        from_GT_To_Snap  (g, saveLoadPath + saveLoadFolder + "/" + graphName + ".snap.graph")
    else:
        return g
    
        
def import_carribean_food_web_graph(save=True, export=True):    
    saveLoadFolder = "Carribean_FoodWeb"    
    graphFile      = saveLoadPath + saveLoadFolder + "/Carribean_Adjacency_Matrix_raw.txt"
    g              = Graph(directed=False)
    edgeWeights    = g.new_edge_property("double")
    counter        = -1
        
    with open(graphFile, "r") as inF:
        for line in inF:
            if line[0] == "#":  # This line is a header. Skip it.
                continue        

            if counter == -1:   # First non header line revels all the species/categories.
                categories = line.split()
                num_nodes = int(len(categories))
                print num_nodes
                counter += 1                
                g.add_vertex(num_nodes)        
                continue
            
            splitted = line.split()
            category = splitted[0]             
            assert(category == categories[counter])
                        
            for neighbor, weight in enumerate(splitted[1:]):
                if weight != "0":
                    e = g.add_edge(g.vertex(neighbor), g.vertex(counter))  # Neighbor eats him.
                    edgeWeights[e] = float(weight)
            counter += 1
    
    taxaToInt = {"D":0, "A":1, "I":2, "F":3, "R":4, "B":5}
    troClass  = g.new_vertex_property("int")    
    for i, categ in enumerate(categories):
        troClass.a[i] = taxaToInt[categ[0]]

    g.vp["trophic_class"] = troClass
    g.ep["edge_weight"]   = edgeWeights

    g = make_simple_graph(g, undirected=True, gcc=True)
    
    graphName = saveLoadFolder
    if save:
        save_data(saveLoadPath + saveLoadFolder + "/" + graphName +".GT.graph", g)
        g.save   (saveLoadPath + saveLoadFolder + "/" + graphName +".graph.xml", fmt="xml")
    if export:    
        from_GT_To_Greach(g, saveLoadPath + saveLoadFolder + "/" + graphName + ".greach.graph")
        from_GT_To_Snap  (g, saveLoadPath + saveLoadFolder + "/" + graphName + ".snap.graph")

    return g

if __name__ == '__main__':
#     import_carribean_food_web_graph(save=True, export=False)
#     import_adj_noun(True)
    pass