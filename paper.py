'''
Created on Aug 11, 2014

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

from myPloting import print_as_latex_table
from graph_tool import topology, clustering
from graph_analysis import IO
from role_discovery import graphs_with_roles as roles
import presentation as comp
import myUtils
from myUtils import p


graphs =  [ "crcl_1", "manufacturing_1", "manufacturing_2", "ml_giant",  "kd_giant", "db_giant",
            "StMarks", "serengeti-foodweb", "Carribean_FoodWeb", "E_Coli", "Mapk", "celegansneural", 
            "polbooks", "polblogs", "netscience", "WebSnapShot", "Sawyer_English",
            "12_Years_Slave", "Plato_Republic", "Zarathustra_English"
          ]

graphRoles = { "crcl_1" : "Synthetic: erdos", "manufacturing_1": "Synthetic: image", "manufacturing_2": "Synthetic: image",
               "ml_giant" : "H-index decile",  "kd_giant" : "H-index decile", "db_giant" : "H-index decile",
               "E_Coli": "Protein lethality", "Mapk":"Cancer association",
               "StMarks" : "Trophic taxa", "serengeti-foodweb":"Trophic taxa", "Carribean_FoodWeb":"Animal class",
               "celegansneural" : "K-core shell", "polbooks" : "K-core shell", "polblogs" : "K-core shell", "netscience" :"K-core shell", "WebSnapShot":"K-core shell",
               "Sawyer_English" : "Part of speech", "12_Years_Slave": "Part of speech", "Plato_Republic": "Part of speech", "Zarathustra_English": "Part of speech"
              }


def paper_graph_name(graphName):
    '''
    Gives a name to each graph appropriate for the presentation of the paper.
    '''
    if graphName[0:4] == "crcl"             : return "CRCL"
    if graphName[0:7] == "Manufac"          : return "Manufact-" + graphName[-1]
    if graphName == "Sawyer_English"        : return "Tom-Sawyer"
    if graphName == "Carribean_FoodWeb"     : return "Caribbean"
    if graphName == "serengeti-foodweb"     : return "Serengeti"
    if graphName == "kd_giant"              : return "KD"
    if graphName == "ml_giant"              : return "ML"
    if graphName == "db_giant"              : return "DB"
    if graphName == "celegansneural"        : return "C-Elegans"
    if graphName == "polbooks"              : return "Pol-Books"
    if graphName == "polblogs"              : return "Pol-Blogs"
    if graphName == "StMarks"               : return "St-Marks"
    if graphName == "netscience"            : return "Net-Science"
    if graphName == "12_Years_Slave"        : return "12-Years-Slave"
    if graphName == "WebSnapShot"           : return "Web-Snap-Shot"
    if graphName == "Zarathustra_English"   : return "Zarathustra"
    if graphName == "Plato_Republic"        : return "Republic"
    if graphName == "Mapk"                  : return "MAPK"

    return graphName


def graph_characteristics(graphName):
    '''
    Return all the characteristics of a graph that we present in the paper.
    '''
    storedFolder            = roles.graph_folder(graphName)
    
    inGraph                 = graph_analysis.IO.load_data("../Data/Graphs/"+storedFolder+"/"+graphName + ".GT.graph").next()
    groupTaxa, blackList    = roles.graph_node_clusters(graphName, inGraph, metric="default")
    res                     = []
    
    headers                 = ["Graph Name", "\#Nodes", "\#Edges", "Edge density", "Clustering Coef.", "Diameter", "Role Taxonomy",  "\#Clusters" ]    
    res.append(paper_graph_name(graphName))
    res.append(inGraph.num_vertices())
    res.append(inGraph.num_edges())
    res.append(round(inGraph.num_edges() / (2*float(inGraph.num_vertices())), 2)  )
    res.append(round(clustering.global_clustering(inGraph)[0], 3) )
    res.append(topology.pseudo_diameter(inGraph)[0])
    res.append(graphRoles[graphName])
    res.append(len(set(groupTaxa)))
    return res, headers
    

def clustering_tables(graphName, strategy=None):
    '''
    Create the data for the two tables of Beta_CV and C-Index.
    '''
    inGraph, groupTaxa, blackList, xTickMarks = cmp.prepare_input_graph(graphName, metric="default")
    gDiameter = topology.pseudo_diameter(inGraph)[0]
    ranks = False

    if blackList != None: print "BlackListed = " + str(blackList)
    allBetaCV = []
    allCIndex = []
#     for methodName in ["roleSim", "simRank", "heatSim", "spectralSim"]:
    for methodName in ["spectralSim"]:
        print methodName
        if cmp.is_spectral(methodName):
            energies = int(gDiameter)            
            methodParams = [energies, 10, strategy]
        else: #non Spectra methods         
            methodParams="default"
        
        if graphName in ["E_Coli", "Carribean_FoodWeb", "Mapk"]  and methodName in ["heatSim", "heatSim_PP" ]:
            distMatrix = cmp.execute_method(methodName, inGraph, graphName, distances=True, ranks=ranks, distFunction="canberra", methodParams=methodParams)
        else:
            distMatrix = cmp.execute_method(methodName, inGraph, graphName, distances=True, ranks=ranks, distFunction="default", methodParams=methodParams)
            
        allBetaCV.append(myUtils.inner_intra_distances(distMatrix, groupTaxa, blackList, ranks=ranks))
        allCIndex.append(myUtils.clustering_c_index(distMatrix, groupTaxa, blackList))
    
    #transform values relative to the second worst
    secondWorse   = sorted(allBetaCV)[-2]
    relativeBetas = [ ((secondWorse-i)/float(i)) * 100 for i in allBetaCV]

    secondWorse   = sorted(allCIndex)[-2]
    relativeCs    = [ ((secondWorse-i)/float(i)) * 100 for i in allCIndex]
    
    return allBetaCV, relativeBetas, allCIndex, relativeCs


#### Print statistics of all graphs
def all_graph_stats():
    allGraphsChars = []
    for graph in graphs:
        graphCharacter, headers  = graph_characteristics(graph)
        allGraphsChars.append(graphCharacter)        
    print_as_latex_table(allGraphsChars, headers)


#### Calculate c-index and betaCV for all graphs
def all_beta_and_cIndex(strategy=None):
    allBetaCV       = []
    relativeBetas   = []
    allCIndex       = []
    relativeCs      = []
    for graph in ["Mapk"]: #graphs:
#         if graph == "WebSnapShot":
#             continue
        myUtils.print_fancy(graph, "-")
        graphName = paper_graph_name(graph)
        res = clustering_tables(graph, strategy)
        allBetaCV.append([graphName] + [round(i,3) for i in res[0]])      #round the raw score to its third decimal
        relativeBetas.append([graphName] + [round(i,2) for i in res[1]])
        allCIndex.append([graphName] + [round(i,3) for i in res[2]])
        relativeCs.append([graphName] + [round(i,2) for i in res[3]])
     
    methodsUsed = ["roleSim", "simRank", "heatSim", "spectralSim"]

    try:
        print_as_latex_table(allBetaCV,     headers = methodsUsed)
        print_as_latex_table(relativeBetas, headers = methodsUsed)
        print_as_latex_table(allCIndex,     headers = methodsUsed)
        print_as_latex_table(relativeCs,    headers = methodsUsed)

    except: #running the script on mad-max
        print methodsUsed
        print allBetaCV
        print relativeBetas
        print allCIndex
        print relativeCs
         

if __name__ == '__main__':
    all_graph_stats()
#     all_beta_and_cIndex()