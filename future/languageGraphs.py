'''
Created on Jul 16, 2014

@author: optas
'''

try:
    import nltk
    from nltk.tokenize.punkt import PunktSentenceTokenizer
except:
    pass
    
from graph_tool import Graph
from collections import defaultdict
import urllib, itertools


def read_content_url(url, utf8=False, gutenbergClean=False, title=None, verbose=True):
    '''
    '''    
    response = urllib.urlopen(url).read()    
    raw = response.decode('utf8') if utf8 else response.decode()
        
    if verbose: print "Number of characters read = %d " % (len(raw), )
    
    if gutenbergClean and title != None:
        start = raw.find(title)
        print start
        print raw[: start + 100]
        raw_input()
        end   = raw.rfind("End of Project Gutenberg") #search from the end of string        
        print raw[end-100, :]
        
    return raw
    
    

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
    print "Content has %d total characters" % (len(content), )
    return content


def sentencesInFile(inFile, utf8=False):
    content    = readContent(inFile)
    if utf8:
        sentences  = PunktSentenceTokenizer().tokenize(content.decode("utf-8"))
    else:
        sentences  = PunktSentenceTokenizer().tokenize(content)
            
    return sentences


def removeAmbigousWords(posTags, ambiguousWords):
    '''
    Remove any word in the set -ambiguousWords- from the dictionary keeping the partOfSpeach for
    each word (key) -posTags-. Removal happens in place.
    '''

    for word in ambiguousWords:
        posTags.pop(word)

def POS_tagging(inFile):
    '''
    Tag the Nouns, Verbs, Adjectives and Adverbs of a corpus.
    Precondition: English Text File
    '''            
    posTags    = dict()     #Keep for each word its part of speech
    ambiguous  = set()      #Words with 2 or more parts of speech

    for s in sentencesInFile(inFile):
        tokens    = nltk.word_tokenize(s)
        tagged    = nltk.pos_tag(tokens)        
        for word, tag in tagged:            
            if tag[0] in ['N', 'V', 'J']: #or tag[0:2] == 'RB':  #Nouns, Verbs, Adjectives, Adverbs
                word = word.lower()                                
                if word in posTags.keys() and posTags[word][0] != tag[0]: 
                    ambiguous.add(word)                
                else:
                    posTags[word] = tag

    removeAmbigousWords(posTags, ambiguous)
    return posTags, ambiguous

                
def keepTaggedTokens(posTags, tokens):
    '''
    Filters the list of input -tokens- so that the returned list contains only tokens that are in 
    the posTags dict, i.e., have been tagged.
    '''
    res = []
    for tok in tokens:
        if tok in posTags:
            res.append(tok)
    return res


def build_co_occurrences(inFile, posTags):
    '''
    For every pair of POS tagged elements (i.e., keys in the posTags) return the non-zero
    counters that contain the number the pair was found together in a sentence of the inFile.
    '''        
    res       = dict()    
    for s in sentencesInFile(inFile):
        tokens    = nltk.word_tokenize(s)
        tokens    = keepTaggedTokens(posTags, tokens)
        if len(tokens) > 1:            
            for tok1, tok2 in itertools.combinations(tokens, 2):
                key = tuple(sorted((tok1, tok2)))                
                if key in res:
                    res[key] += 1
                else: res[key] = 1    
    
    return res

def part_of_speech_int_map(posToInt=True):
    '''
    Codes every part of speech (string) to a unique int identifier and vice versa.
    Nouns = 0
    Verbs = 1
    Adjectives = 2
    Adverbs = 3
    '''
    if posToInt:
        return {'N':int(0), 'V':int(1), 'J':int(2), 'R':int(3)}
    else:
        return {0:"Noun", 1:"Verb", 2:"Adjective", 3:"Adverb"}
    
def build_cooc_graph(coDic, posMap):
    '''    
    Converts a dictionary keeping word-word co-occurrences to a graph object where the edge weight
    between nodes (words) corresponds to their co-occurrence.

    Args:
        coOcDic: A dictionary where keys are tuples with 2 elements, (string1, string2) corresponding to the co-occurrence of two words. 
        The corresponding value is an integer capturing the times the co-occurrence happened (e.g., through out a book).
        posMap: A dictionary from word to Part Of Speech.
    Returns:
        A graph object.
    
    '''
    g               = Graph(directed=False)
    wordToNodeID    = dict()                           #maps a word to the ID of the node that it will be stored    
    eWeight         = g.new_edge_property("int")       #edges have weights capturing number of co-occurrences
    words           = g.new_vertex_property("object")  #keep for each node the (potentially unicode) corresponding word as an attribute
    POS             = g.new_vertex_property("string")  #keep the Part Of Speech
    nodeID = 0
    for word1, word2 in coDic.keys():                   #Each key is a (noun, noun) string. It will become an edge 
        if word1 not in wordToNodeID:
            wordToNodeID[word1] = nodeID
            v = g.add_vertex()
            assert(str(v) == str(nodeID))
            words[v] = word1
            POS[v]   = posMap[word1]
            nodeID  += 1
        if word2 not in wordToNodeID:
            wordToNodeID[word2] = nodeID
            v = g.add_vertex()
            assert(str(v) == str(nodeID))
            words[v] = word2
            POS[v]   = posMap[word2]
            nodeID += 1        
        source = wordToNodeID[word1]
        target = wordToNodeID[word2]
        e = g.add_edge(source, target)
        eWeight[e] = coDic[(word1, word2)]

    g.edge_properties["co-occurrence"]   = eWeight
    g.vertex_properties["word"]          = words
    g.vertex_properties["partOfSpeach"]  = POS
    
    #Encode the POS as a short number
    POS_encoded = g.new_vertex_property("short") 
    posEncoder  = part_of_speech_int_map(posToInt=True)

    for v in g.vertices():
        POS_encoded[v] = posEncoder[POS[v][0]]

    g.vertex_properties["partOfSpeach_encoded"]  = POS_encoded
    
    return g

def print_co_pos_stats(co_oc_Dic, posTags):       
    wordsOfCategory   = defaultdict(set)
    edgesOfCategory   = defaultdict(lambda : 0)
    weightOfCategory  = defaultdict(lambda : 0)
    
    for word1, word2 in co_oc_Dic.keys():
        pos1 = posTags[word1][0]
        pos2 = posTags[word2][0]
        wordsOfCategory[pos1].add(word1)
        wordsOfCategory[pos2].add(word2)
        edgesOfCategory[tuple(sorted([pos1,pos2]))] += 1
        weightOfCategory[tuple(sorted([pos1,pos2]))] +=  co_oc_Dic[(word1,word2)]
    
    print "#######~~~~~~~~~~~~########"
    for categ in wordsOfCategory.keys():
        print categ + " has %d members" %(len(wordsOfCategory[categ]),)    
    print "#######~~~~~~~~~~~~########"
    for categ1, categ2 in edgesOfCategory.keys():
        if categ1 != categ2:
            edgeDensity = edgesOfCategory[(categ1, categ2)] / float(len(wordsOfCategory[categ1]) * len(wordsOfCategory[categ2]))
        else:
            edgeDensity = 2*edgesOfCategory[(categ1, categ2)] / float(len(wordsOfCategory[categ1]) * (len(wordsOfCategory[categ2])-1))
            
        print "Between %s words  and %s words, the Number of co-occurrences = %d, the number of edges =%d and the edge density = %f." %(categ1, categ2, weightOfCategory[(categ1, categ2)], 
                                                                                                                                        edgesOfCategory[(categ1, categ2)], edgeDensity)        
    print "#######~~~~~~~~~~~~########"
        
def extract_POS_Graph(novelFile, novelName, tryLoad=False, save=False, verbose=True):
    if tryLoad:        
        try:                    
            posTags, ambiguous = IO.load_data("../Data/Novels/"+novelName+".posTags")
            print "posTags + ambiguous words were Loaded."
        except:
            print "posTags + ambiguous are being Calculated."
            posTags, ambiguous = POS_tagging(novelFile)
    else:
            print "posTags + ambiguous are being Calculated."
            posTags, ambiguous = POS_tagging(novelFile)
            
    if save:
        IO.save_data("../Data/Novels/"+novelName+".posTags", posTags, ambiguous)
            
    co_oc_Dic  = build_co_occurrences(novelFile, posTags)    
    if verbose:
        print_co_pos_stats(co_oc_Dic, posTags)
    g = build_cooc_graph(co_oc_Dic, posTags)    
    IO.make_simple_graph(g)

    return g
    
def test_extract_POS_Graph(novelFile, novelName):
    g = extract_POS_Graph(novelFile, novelName, tryLoad=True, verbose=True)

    for v in g.vertices():
        nv = [int(e.target()) for e in v.out_edges()]
        assert(len(nv) == len(set(nv)))      #no parallel edges 
  

    Lg  = spectral.laplacian(g, weight=g.ep['co-occurrence'])        
    diag = Lg.diagonal()                
    for i, v in enumerate(g.vertices()):      
        w = 0
        for e in v.out_edges(): w  +=  g.ep['co-occurrence'][e]
        assert(w == diag[i])        #The laplacian diagonal reflects the weigh as captures by total co-occurrences of word i
    #TODO
    #All the statistics that you print with  print_co_pos_stats
    #cross check them with those extracted from the graph object
    #Like
#     print len(np.where(g.vertex_properties["partOfSpeach_encoded"].a == 0)[0])
#     print len(np.where(g.vertex_properties["partOfSpeach_encoded"].a == 1)[0])
#     print len(np.where(g.vertex_properties["partOfSpeach_encoded"].a == 2)[0])
#     
#     NJ_filter = g.new_vertex_property("bool")
#     NJ_filter.a[np.where(g.vertex_properties["partOfSpeach_encoded"].a == 0)[0]] = 1
#     NJ_filter.a[np.where(g.vertex_properties["partOfSpeach_encoded"].a == 1)[0]] = 1
#     
#     NJ = GraphView(g, vfilt=NJ_filter)


# edgeW=None g.ep['co-occurrence'[]

from graph_analysis import IO
from time import time
from graph_tool import spectral, GraphView
import numpy as np
from role_discovery import graph_kernels as gk
from sklearn.cluster import KMeans


if __name__ == "__main__":        
    novelName = "Zarathustra_English"
    novelFile = "../Data/Novels/"+novelName+".txt"
    saveGraph = True; tryLoadPOS=True; savePOS=True; verbose=True;    
    g = IO.compute_or_load(extract_POS_Graph, "../../Data/Graphs/novels/"+novelName+".GT.graph", saveGraph, novelFile, novelName, tryLoadPOS, savePOS, verbose)
    print g
#     evals, evecs = gk.laplacian_spectrum(g, novelName, eigsToCompute = "all", edgeW=None, tryLoad=True, save=True)
#     print evecs.shape    
    