'''
Created on May 23, 2014

@author: optas
'''

import sys, nltk, itertools, codecs
from graph_analysis import IO
import numpy as np 
import myPloting as myPlt
from matplotlib.mlab import PCA
from role_discovery import graph_kernels as gk
from collections import defaultdict
from nltk.tokenize.punkt import PunktWordTokenizer
from nltk.tokenize.punkt import PunktSentenceTokenizer
from synthetic_graphs.synthGraphs import *
import scipy.cluster

try:
    from graph_tool import collection, topology, spectral, stats
    from graph_tool import stats as gpStats
    from graph_tool import *
except:
    "Graph_Tool is not installed in this environment."




def assignToKMeans(sig, centroids):
    res = []
    from sklearn.metrics.pairwise import euclidean_distances
    distanceMatrix = euclidean_distances(sig, centroids)
    for i, pt in enumerate(distanceMatrix):        
        print distanceMatrix[i]
        raw_input()
#         a.index(max(a))
#         res.append(nsmallest(kneighbors, range(distanceMatrix.shape[0]), key=distanceMatrix[i,:].__getitem__))
#         return topMatches
    


def readContent(inFile):
    '''
    Read inFile (e.g., a novel) and return all its lines as a single string from which the newlines
    have been replaced by a space.
    '''
    content = []
    with open(inFile, "r") as fIn:
        for line in fIn:
            content.append(line.strip()+" ")

    content = ''.join(content)  #keep content as a big string.

    return content
    
    

def partOfSpeech(inName, inFile, lang):
    content = readContent(inFile)
    print "Total characters of "+ inName + ": ", len(content)
    
    if lang == "English":
        sentences = PunktSentenceTokenizer().tokenize(content.decode("utf-8"))
#         with open("../Data/"+lang+"_"+inName+"_sentence_split.txt", "w") as fout:
#             for s in sentences:            
#                 fout.write(s); fout.write("\n")
        tagger = nltk.pos_tag        
    else:
        raise NotImplementedError
                    
    res        = dict()
    posMap     = dict()     #keep for each word its part of speech
    doubleMean = set()      #doubleMenaing: 2 or more part of speech roles
    coOc = []               #list with co-occurrences in sentences

    for s in sentences:
        tokens = nltk.word_tokenize(s)
        tagged = tagger(tokens)
        posInSent = set()

        for word, tag in tagged:            
            if tag[0] == 'N' or tag[0] == 'V': #Keep Nouns and Verbs
                word = word.lower()                                
                if word in posMap.keys() and posMap[word][0] != tag[0]: 
                    doubleMean.add(word)                
                    continue                
                posInSent.add(word)                
                posMap[word] = tag
        if len(posInSent) > 1:
            coOc.append(posInSent)
    
    print "Words that were both Verbs and Nouns and will be discarded: ", len(doubleMean)        

    fout = open("../Data/" + lang + "_" + inName +"_NounVsVerb.txt", 'w')
    for s in coOc:
        clearCoOc = []
        for word in s:
            if word not in doubleMean:
                clearCoOc.append(word)
        
        if len(clearCoOc) > 1:
            print >> fout, clearCoOc            
            for n1, n2 in itertools.combinations(clearCoOc, 2):
                key = tuple(sorted((n1,n2)))                
                if key in res:
                    res[key] += 1
                else: res[key] = 1    
    
    return res, posMap



def coOcDicToGraph(coOcDic, posMap):
    '''    
    Converts a dictionary keeping word-word co-occurrences to a graph object.

    Args:
        coOcDic: A dictionary where keys are tuples with 2 elements, (string1, string2) corresponding to the co-occurrence of two words. 
        The corresponding value is an integer capturing the times the co-occurrence happened (e.g., through out a book).
        posMap: A dictionary from word to Part Of Speech.
    Returns:
        A graph object.
    
    '''
    nounToNodeID = dict()                           #maps a noun (string) to an ID (int)
    g            = Graph(directed=False)
    eWeight      = g.new_edge_property("int")       #edges have weights capturing number of co-occurrences
    vNoun        = g.new_vertex_property("object")  #node is a (potentially unicode) string     
    pos          = g.new_vertex_property("string")
    nodeID = 0
    for edge in coOcDic.keys():                   #Each key is a (noun, noun) string. It will become an edge 
        if edge[0] not in nounToNodeID:
            nounToNodeID[edge[0]] = nodeID
            v = g.add_vertex()
            assert(str(v) == str(nodeID))
            vNoun[v] = edge[0]
            pos[v]   = posMap[edge[0]]
            nodeID += 1
        if edge[1] not in  nounToNodeID:
            nounToNodeID[edge[1]] = nodeID
            v = g.add_vertex()
            assert(str(v) == str(nodeID))
            vNoun[v] = edge[1]
            pos[v]   = posMap[edge[1]]
            nodeID += 1        
        source = nounToNodeID[edge[0]]
        target = nounToNodeID[edge[1]]
        e = g.add_edge(source, target)
        eWeight[e] = coOcDic[edge]

    g.edge_properties["weights"]   = eWeight
    g.vertex_properties["noun"]    = vNoun    
    g.vertex_properties["partOs"]  = pos
    return g



def printNounVerbGraphStats(g):
    print "%d Nodes, %d Edges, IsDirected? %s." % (g.num_vertices(), g.num_edges(), g.is_directed())
    nouns = verbs = 0
    
    for v in g.vertices(): 
        if g.vp['partOs'][v][0] == "N":
            nouns +=1
        elif g.vp['partOs'][v][0] == "V":
            verbs +=1
        else:
            assert(False)
    print "From the %d vertices %d, are Nouns and %d are Verbs." % (g.num_vertices(), nouns, verbs)

    verbVerb = nounNoun = nounVerb = totalW = 0
    for e in g.edges():
        totalW += g.ep['weights'][e]
        v1 = e.source(); v2 = e.target()
        if g.vp['partOs'][v1][0]    == "N" and  g.vp['partOs'][v2][0] == "N":
            nounNoun +=1 
        elif (g.vp['partOs'][v1][0] == "N" and g.vp['partOs'][v2][0] == "V" ) or (g.vp['partOs'][v1][0] == "V" and   g.vp['partOs'][v2][0] == "N"):
            nounVerb +=1
        else: verbVerb +=1
    
    print "Interactions:\n%d Noun-Noun Edges \n %d Noun-Verb Edges \n %d Verb-Verb Edges" % (nounNoun , nounVerb, verbVerb)
    print "Total Edge Weight = "+ str(totalW)
    



novel = "Zarathustra"
# novel = "Sawyer"
lang   = "English"
inFile = "../Data/" + novel +"_" + lang + ".txt"
# coOcDic, posMap = partOfSpeech(novel, inFile, lang)
# IO.saveData(novel+"_"+lang+"_nounVSverb.data" , coOcDic, posMap)
coOcDic, posMap = graph_analysis.IO.load_data(novel+"_"+lang+"_nounVSverb.data")

g = coOcDicToGraph(coOcDic, posMap)
l = topology.label_largest_component(g)  #Keep Largest Connected Component
g.set_vertex_filter(l); g.purge_vertices()
printNounVerbGraphStats(g)

nodesColored = []
nodesLabels  = []
palet = palets("posToColor")

allTypes = set()
for v in g.vertices():
    nodesColored.append(palet[g.vp['partOs'][v][0]])
    allTypes.add(g.vp['partOs'][v])

#Make binary label for SVM
for v in g.vertices():
    if g.vp['partOs'][v][0] == 'N':
        nodesLabels.append(1)
    else:
        nodesLabels.append(0)

eigsToCompute = 500
eigsToUse     = 50
timePoints    = 10
save  = True
small = True

eigsWhere       = "../Data/Spectrums/SM/" + novel + "_" + lang +"_NounVsVerb_"+"_Eigs_" + str(eigsToCompute + 1)
Lg              = spectral.laplacian(g, weight=g.ep['weights'])
evals, evecs    = graph_analysis.IO.compute_or_load(gk.computeSpectrum, eigsWhere, save, Lg, eigsToCompute+1, small=small)

timeSamples     = gk.logTimeScale(evals, eigsToUse, timePoints)    
sig             = gk.HKS(evals, evecs, eigsToUse, timeSamples)



sig = PCA(sig).Y[:, 0:3]
# print sigsPCA.shape
# p()

# 
# x = []; y = []; z = [];
# labels = kcore.a.tolist()
# labelsPruned = list()
# for i, item in enumerate(sigsPCA.Y):
#     if labels[i] in [1,2,3]:
#         x.append(item[0])
#         y.append(item[1])
#         z.append(item[2])
#         labelsPruned.append(labels[i])



import sklearn.cluster
km = sklearn.cluster.KMeans(n_clusters=2, init='k-means++')
predictions = km.fit_predict(sig)
 
print len(np.where(predictions == 0)[0])
print len(np.where(predictions == 1)[0])
# print len(np.where(predictions == 2)[0])
# print len(np.where(predictions == 3)[0])



from sklearn.neighbors import KNeighborsClassifier
knn = KNeighborsClassifier()


knn.fit(sig[0:4000,:], nodesLabels[0:4000])
KNeighborsClassifier(algorithm='auto', leaf_size=30, metric='minkowski', n_neighbors=1, p=4, weights='uniform')
yP = knn.predict(sig[4001:])

print sum (yP)
print len(yP)
print sum(yP == nodesLabels[4001:]) / float(len(yP))


# 
# X = sig
# y = np.array(nodesLabels)
# from sklearn.svm import SVC
# 
# clf = SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
# gamma=0.0, kernel='rbf', max_iter=-1, probability=False,
# random_state=None, shrinking=True, tol=0.001, verbose=False)
# clf.fit(sig, y)






# sigsPCA         = PCA(sig) 
# print "Fraction of Variance in the first three PCA dimensions: " +str(sum(sigsPCA.fracs[0:3]))
# print "Variance on 1dt Dim = ", sigsPCA.fracs[0]
# myPlt.plot_PCA_Res(sigsPCA, dim=1, legendLabel=None, colorMap=[nodesColored, palet], save=False, saveAt=None)    


 
