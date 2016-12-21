'''
Created on Jul 10, 2014

Contact: pachlioptas@gmail.com

Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import numpy as np
try:
    from graph_tool import draw
    from myUtils import cluster_based_on_percentile
    from future.languageGraphs import part_of_speech_int_map
except:
    pass
      
#Keep the names of all the graphs represented in this modules
allGraphs = ["ml_giant",  'kd_giant', 'db_giant', "serengeti-foodweb", "Carribean_FoodWeb"
             "crcl_1", "manufacturing_1", "manufacturing_2", "StMarks", "polblogs", "celegansneural", "E_Coli",
             "netscience", "adjnoun", "polbooks", "polblogs", "Sawyer_English", "12_Years_Slave", "WebSnapShot",
             "Mapk", "Plato_Republic", "Zarathustra_English", "stratified_1"]


def stratified_node_clusters(inGraph):        
    blackList = None
    return inGraph.vp["classes"].a, blackList


def crcl_node_clusters(inGraph, attach=True):
    '''
    Roles = Erdos - Cliques - Circles - Lattices.
    attach == False => nodes used to attach nodes from different classes will be blackListed.    
    '''
    allTaxa   = inGraph.new_vertex_property("short")    
    allTaxa.a = inGraph.vp['bb'].a                          #bb = backbone nodes pf the erdos-renyi random graph
    allTaxa.a[np.where(inGraph.vp['cliq'].a == 1)[0]] = 2
    allTaxa.a[np.where(inGraph.vp['circ'].a == 1)[0]] = 3
    allTaxa.a[np.where(inGraph.vp['late'].a == 1)[0]] = 4

    if attach:
        blackList = None
    else:        
        blackList = np.where(inGraph.vp['attachments'].a == 1)[0]
    return allTaxa.a, blackList


def manufacturing_node_clusters(inGraph):
    '''
    inGraph: A graph representing a manufacturing process.
    It is assumed that was created via a random blockModel [buildModel(manufacturingImage())].
    '''
    blackList = None
    return inGraph.vp["blockIDs"].a, blackList


def kCore_node_clusters(inGraph):
    '''
    Taxa are the k-core shells and no nodes are blackListed
    '''
    blackList = None            #Clean
    return inGraph.vp["kcore"].a, blackList


def adj_noun_clusters(inGraph):
    '''
    Adjective or Noun: it refers to the small example imported from GT (which is imported by Newman).
    '''
    blackList = None
    vals = inGraph.vp['value'].a
    vals = np.array(vals, dtype=np.int16)    
    return vals, blackList


def cluster_based_on_percentile(a, perMin=0, perMax=100):    
    '''
    Categorize each entry of a in a percentile.
    '''
    percentiles = list()
    
    for percent in range(perMin, perMax, 10):
        pd = np.percentile(a, percent)        
        percentiles.append(pd)        
#         print percent, pd, sum(a <= pd)

    percentiles = np.unique(percentiles)
    percentClass = np.digitize(a, bins=percentiles)
#     labels = np.unique(percentClass)

#     print labels
#     for i in labels:
#         print i 
#         print len(np.where(percentClass == i)[0])
#         print np.unique(a[np.where(percentClass == i)[0]])


    return percentClass


def co_author_node_clusters(inGraph, metric = "g-Index"):
    blackList = None
    if  metric == "g-Index":
        return cluster_based_on_percentile(inGraph.vp["g-Index"].a), blackList
    else:
        return cluster_based_on_percentile(inGraph.vp["h-Index"].a), blackList


def serengeti_node_clusters(inGraph):
    groupTaxa    = inGraph.vp["group"].a
    blackList    = None
    return groupTaxa, blackList


def stMarks_node_clusters(inGraph, metric="Trophic Taxa"):
    if metric == "mass":
        groupTaxa    = cluster_based_on_percentile(inGraph.vp["mass"].a)
    else:
        groupTaxa    = inGraph.vp["trophicClass"].a
    blackList    = None
    return groupTaxa, blackList


def carribean_food_web_node_clusters(inGraph, metric="Trophic Taxa"):
    groupTaxa    = inGraph.vp["trophic_class"].a
    blackList    = None
    #TODO - implement blacklist and top predators vs. least (see Yu email)
    return groupTaxa, blackList


def novel_node_clusters(inGraph):
    groupTaxa    = inGraph.vp["partOfSpeach_encoded"].a
    blackList    = None
    return groupTaxa, blackList


def e_coli_node_clusters(inGraph, metric):
    return inGraph.vp["lethality"].a, None


def mapk_node_clusters(inGraph, metric):
    blackList = [0]                     # 0 are the nodes for which we do NOT have a class. 7 corresponds to the "non clear" class of Enzymes.
    return inGraph.vp["Classes"].a, blackList
    
def taxa_names(graphName):
    if graphName in ['ml_giant', 'kd_giant', 'db_giant']:
        taxaNames = ["10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%"]
    elif graphName[:4] == "crcl":
        taxaNames = ["Erdos", "Clique", "Circle", "Lattice"]
    elif graphName[:8] == "manufact":
        taxaNames = ["Supply", "Assembly", "Delivery", "Shop", "Warehouse", "Consumer"]    
    elif graphName == "serengeti-foodweb" or "StMarks":
        taxaNames =  ["Taxa_"+str(i) for i in range(25)]
    elif graphName == "Carribean_FoodWeb":
        taxaNames =  ["Detritus", "Autotroph", "Invertebrate", "Fish", "Reptile", "Bird"]
    elif graphName in ["WebSnapShot" , "polbooks", "polblogs", "netscience", "celegansneural"]:
        taxaNames = ["Shell_"+str(i) for i in range(1, 30)]
    elif graphName == "adjnoun":
        taxaNames = ["Noun", "Adjective"]
    elif graphName in ["12_Years_Slave", "Sawyer_English", "Zarathustra_English", "Plato_Republic"]:
        taxaNames = []
        posMap = part_of_speech_int_map(posToInt=False)
        allKeys = posMap.keys()
        for posType in sorted(allKeys):
            taxaNames.append(posMap[posType])
    elif graphName == "E_Coli":
        return ["Non Lethal", "Lethal"]
    elif graphName == "Mapk":
        return ["TF", "MAPK-L1", "MAPK-L2", "MAPK-L3", "MAPK-L4", "Membrane", "Enzyme1", "Enzyme2"]
    else: 
        assert(False)
        
    return taxaNames


def graph_node_clusters(graphName, inGraph, metric):
    if graphName in ['ml_giant', 'kd_giant', 'db_giant']:
        return co_author_node_clusters(inGraph, metric)
        
    elif graphName == "serengeti-foodweb":
        return serengeti_node_clusters(inGraph)

    elif graphName[:4] == "crcl":
        return crcl_node_clusters(inGraph, attach=True)
    
    elif graphName[:8] == "manufact":
        return manufacturing_node_clusters(inGraph)
    
    elif graphName[:10] == 'stratified':
        return stratified_node_clusters(inGraph)

    elif graphName in ["WebSnapShot", "polblogs" , "polbooks",  "netscience", "celegansneural"]:        
        return kCore_node_clusters(inGraph)
            
    elif graphName == "adjnoun":
        return adj_noun_clusters(inGraph)
    
    elif graphName in ["Zarathustra_English", "Sawyer_English", "12_Years_Slave", 'Plato_Republic']: 
        return novel_node_clusters(inGraph)
        
    elif graphName == "StMarks":
        return stMarks_node_clusters(inGraph, metric)

    elif graphName == "Carribean_FoodWeb":
        return carribean_food_web_node_clusters(inGraph, metric)
    
    elif graphName == "E_Coli":
        return e_coli_node_clusters(inGraph, metric)
    
    elif graphName == "Mapk":
        return mapk_node_clusters(inGraph, metric)
        
    else:
        assert(False)

def graph_folder(graphName):
    '''
    Returns the name of the folder the graph resides. I.e., ../Data/Graphs/<returned_value>/
    '''    

    if graphName in ['ml_giant', 'kd_giant', 'db_giant']:
        return  'Co-authors'
    elif graphName in ['serengeti-foodweb', 'adjnoun', 'polblogs', 'polbooks', 'netscience', 'celegansneural']:
        return 'GT_Coll'
    elif graphName[:len('crcl')] == 'crcl' or graphName[:len('manufact')] == 'manufact' or graphName[:len('stratified')] == 'stratified':
        return 'Synthetic'
    elif graphName == "WebSnapShot" or graphName == "StMarks" or graphName =="Carribean_FoodWeb"\
        or graphName =="Mapk" :  #Folder and Graph have same name
        return graphName
    elif graphName == "12_Years_Slave" or graphName == "Sawyer_English" or graphName == "Zarathustra_English" or graphName == "Plato_Republic":
        return 'Novels'
    elif graphName[:6] == "E_Coli":        
        return "E_Coli"    
    else:
        assert(False)

def plot_roles_on_graph(inGraph, graphName, metric="default"):
    nodeClusters, blackList = graph_node_clusters(graphName, inGraph, metric=metric)
#     print type(nodeClusters)
#     print nodeClusters
#     raw_input()
    #     pos = draw.sfdp_layout(inGraph)    
    nodeClustersProp  = inGraph.new_vertex_property('float')
    nodeClustersProp.a = np.array(nodeClusters)    

    draw.graph_draw(inGraph, vertex_fill_color = nodeClustersProp, edge_color="black", output = graphName + "_roles.pdf")