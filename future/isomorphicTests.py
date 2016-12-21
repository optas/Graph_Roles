'''
Created on Mar 24, 2014

@author: optas
'''


from myUtils import *
from graph_analysis import IO
import sys
import myPloting as plt
import pylab as plb
from scipy import spatial, sparse
from scipy.sparse import linalg
import random, scipy, math
from role_discovery.graph_kernels import *
from graph_tool.all import *



# Mental exercise for the usage of the forward/backward matching
#    Lg = spectral.laplacian(g).todense()
#    Lp = spectral.laplacian(gp).todense()
#
#     for i in xrange(Lg.shape[0]):
#         assert (Lg[i,i] == Lp[forwardMap[i], forwardMap[i]])
#     for i in xrange(Lp.shape[0]):
#         assert (Lp[i,i] == Lg[perm[i],perm[i]])    



   
 
    

    
   



def spectralStability(exampleG=None, nstep=None):

    eigsToUse = 11
    save = False
    kernel = heatKernel
    maxTopK = 16
    topKs = [1, 2, 4, 8, maxTopK]
    gcc = True
    smallest = True
#     
    if exampleG == None:
#         exampleG  = 'cond-mat-2003'
        exampleG  = 'email-Enron'
                
# #         = raw_input("Test module\n Give Graph To Load = ")
    if nstep == None:    
        nstep  = 10
#         #= int(raw_input("nstep"))
#     
#   

#     Lg, Lgp, forwardMap = IO.compute_or_load(IO.makeISOLaplacians, "../Data/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", False, exampleG, gcc)        
#     Mg                  = IO.compute_or_load(IO.mcSherryFromLaplacian, "../Data/Laplacians/"+exampleG+"_McSher", save, None)
#     Mgp                 = IO.compute_or_load(IO.mcSherryFromLaplacian, "../Data/Laplacians/"+exampleG+"_McSher_Iso", save, Lgp)
#          

#     Simple examples#
#     g = generation.lattice((10,10,10))
#     g = IO.loadGraph(exampleG, True, True, True)
#     dist, ends = topology.pseudo_diameter(g)
#     print dist
#     print ends
#     p('geia')


#     Ag = spectral.adjacency(g).todense()  
#     gp, forwardMap = IO.createIsoGraph(Ag, assertIso=True, inGraph = g)
# 
#     Lg  = spectral.laplacian(g)  #create the Laplacians  
#     Lgp = spectral.laplacian(gp)
    
#     Mg  = IO.mcSherryFromGraph(g)
#     Mgp = IO.mcSherryFromGraph(gp)

    
#     evals1, evecs1 = computeSpectrum (Lg, eigsToUse+1, smallest)
    
    
#     IO.save_data("../Data/Spectrums/"+exampleG+"McSherry_Eigs_"+str(eigsToUse), evals1, evecs1)
    evals1, evecs1 = graph_analysis.IO.load_data("../Data/Spectrums/Smallest_sigma/"+exampleG+"_GCC_"+str(eigsToUse)+"_Eigs")
    print "Eigenvalues ", str(evals1)
    p("geias")
#     evals2, evecs2 = computeSpectrum (Lgp, eigsToUse+1, smallest)
#     
#     np.flipud(evecs1)
#     np.flipud(evecs2)
#     evals1 = evals1[::-1]
#     evals2 = evals2[::-1]
# 
#     timeSamples1 = logTimeScale(evals1, eigsToUse, nstep)
#     sig1 = IO.compute_or_load(HKS, "../Data/", False, evals1, evecs1, eigsToUse, timeSamples1)        
#     timeSamples2 = logTimeScale(evals2, eigsToUse, nstep)                
#     sig2 = IO.compute_or_load(HKS, "../Data/_ISO", False, evals2, evecs2, eigsToUse, timeSamples2)        
# 
#     print sig1.shape
#     for i,pt in enumerate(sig1):
#         print pt
#         print sig2[forwardMap[i]]
        
# 
#     nodeMapping = IO.compute_or_load(signatureBasedMapping, "skata", False, sig1, sig2, kernel, 16)                            
#     resN = mappingAccuracy(nodeMapping, forwardMap, topKs)
#     print resN
    


#     print evecs1.shape#     
#     for eig in xrange(eigsToUse):
#         a = set(evecs1[:, eig])
#         b  = set([str(i) for i in evecs1[:,eig]])
#     print eig, evals1[eig], len(a), len(b)
#         p()
#     p()
#         print "eig " + str(eig+1) +" lambda=" +str(evals1[eig]) +" "+ str(evals2[eig])
#         for i in xrange(len(forwardMap)):
#             print evecs1[i,eig], evecs2[forwardMap[i],eig]
 
#     evecs1 = abs(evecs1)
#     evecs2 = abs(evecs2)
# 
    print evecs1.shape
     
    eps = 1e-18
    
    print "epsilon", str(eps)

    for k in xrange(1, eigsToUse):                                             
 
        sig1 = evecs1[:, eigsToUse-1-k:eigsToUse-1]
 
#         print sig1[0]
#         print sig1[1]
#         p()
          
# #         print k, sig1.shape
# #         print evals1[eigsToUse-1-k:eigsToUse-1]
# #         print evals1
#         sig2 = evecs2[:,eigsToUse-1-k:eigsToUse-1]
         
#         print "points live in ", sig1.shape
        KDTree=spatial.cKDTree(sig1, leafsize=100) #TODO: Play with the leafsize to get the fastest result for your dataset
 
        n = sig1.shape[0]    
        uniq = 0
        for i in xrange(n):
            pt = sig1[i,:]
#             print pt
            knn = KDTree.query(pt, k=2, distance_upper_bound=eps)
#             print knn
#             p()
            if np.inf in knn[0]: uniq +=1
#             print i, knn[1]
#             p()
 
         
        print "#eigs = "+str(k)+ " unique nodes "+ str(uniq)
                     
#         

#                                        
#         print mappingAccuracy(nodeMapping, forwardMap, topKs)

    

#     print evecs1.shape
#     print evecs2.shape
#     p()
#     
#     print evals1
#     print evals2
#     

    
    
    
   
   
def computeISOSpectrums(save = True):     
    exampleG  = raw_input("Test module\n Give Graph To Load = ")
    eigsToUse = int(raw_input("How many Eigenvalues?"))
    print "Eigenvalues that will be computed :" + str(eigsToUse)
                                 
    gcc = True
    smallest = True
    Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", False, exampleG, gcc)
     

    eigs1Where ="../Data/Spectrums/"+exampleG+"_GCC_"+str(eigsToUse+1)+"_Eigs"                
    evals1, evecs1 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs1Where, save, Lg, eigsToUse+1, smallest)
         
    eigs2Where ="../Data/Spectrums/"+exampleG+"_ISO_GCC_"+str(eigsToUse+1)+"_Eigs"
    evals2, evecs2 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs2Where, save, Lgp, eigsToUse+1, smallest) 
     




def plotEigDistribution(exampleG=None, allEigs=None, smallest=True, save=True):
    if exampleG == None:
        exampleG = raw_input("Test module\n Give Graph To Load = ")
    if allEigs == None:    
        allEigs = int(raw_input("#Eigs? "))    
    print "Compute *Smallest* Eigenvalues:" + str(smallest)
    
    gcc = True
    
    Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", False, exampleG, gcc)

    if smallest:  #smallest Eigenvalues    
        eigsSmalWhere ="../Data/Spectrums/Smallest_sigma/"+exampleG+"_GCC_"+str(allEigs+1)+"_Eigs"
        print "computing " +str(allEigs+1)                
        evalsSmal, evecsSmal = graph_analysis.IO.compute_or_load(computeSpectrum, eigsSmalWhere, save, Lg, allEigs+1, smallest)
        assert(len(evalsSmal) == allEigs+1) 
        evalToPlot = evalsSmal[1:]
        label = "Smallest"
        plb.figure()
        plb.ylim(ymin = evalToPlot[0]-0.5, ymax = evalToPlot[-1] + 0.5)
    else:        
        eigsLargWhere ="../Data/Spectrums/LargestMagn/"+exampleG+"_GCC_"+str(allEigs)+"_Eigs"                
        evalsLarg, evecsLarg = graph_analysis.IO.compute_or_load(computeSpectrum, eigsLargWhere, save, Lg, allEigs, smallest)
        assert(len(evalsLarg) == allEigs) 
        evalToPlot = evalsLarg[::-1]  #make it be decreasing
#         evalToPlot = evalsLarg
        label = "Largest"
        plb.figure()
        plb.ylim(ymin = evalToPlot[-1]-0.5, ymax = evalToPlot[0] + 0.5)
                
    
    plb.xlim(xmin = 0.5, xmax = allEigs + 0.5)

    if allEigs <= 20:
        plb.xticks(range(1, allEigs+1), range(1, allEigs+1))    
        plb.plot(range(1, allEigs+1), evalToPlot)
    else:
        
        plb.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off') #turn of the ticks
        plb.plot(range(1, allEigs+1), evalToPlot, "*")
    
     
    plb.xlabel("Rank of Eigenvalue")
    plb.ylabel("Magnitude")
    plb.title("<Graph = "+exampleG+">   "+label+" Eigenvalues Ranking")
               
    plb.savefig(exampleG+"_"+label+"_"+str(allEigs)+"_Eigen_Distribution.png")
    

    
     
     
def accurasyWithTimeSamples(exampleG=None, allEigs=None):
    if exampleG == None:
        exampleG = raw_input("Test module\n Give Graph To Load = ")
    if allEigs == None:    
        allEigs = int(raw_input("#Eigs? "))
        
    print exampleG, allEigs
    
    allNsteps = range(1,11)
    allNsteps.extend([20,30,40,50,60,70,80,100,200])

    kernel = heatKernel    
    
    maxTopK = 2 ; topKs = [1, maxTopK]

    save = False
    gcc = True
    smallest = True
     
    Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", False, exampleG, gcc)
     

    #---- Load MAX Number of Eigenvalues -----
    eigs1Where ="../Data/Spectrums/Smallest_sigma/"+exampleG+"_GCC_"+str(allEigs+1)+"_Eigs"                
    evals1, evecs1 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs1Where, False, Lg, allEigs+1, smallest)
     
    eigs2Where ="../Data/Spectrums/Smallest_sigma/"+exampleG+"_ISO_GCC_"+str(allEigs+1)+"_Eigs"
    evals2, evecs2 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs2Where, False, Lgp, allEigs+1, smallest) 
    
    allAccuracies = []
                   
    for nstep in allNsteps:        
        saveStamp = "HKS/"+exampleG+"_GCC_eigs_"+str(allEigs)+"_guibasTime_"+str(nstep)
              
        timeSamples1 = logTimeScale(evals1, allEigs, nstep)
        sig1 = graph_analysis.IO.compute_or_load(HKS, "../Data/sigs/"+saveStamp, save, evals1, evecs1, allEigs, timeSamples1)
        #--- ISOmorphic case ----          
        timeSamples2 = logTimeScale(evals2, allEigs, nstep)                
        sig2 = graph_analysis.IO.compute_or_load(HKS, "../Data/sigs/"+saveStamp+"_ISO", save, evals2, evecs2, allEigs, timeSamples2)

        nodeMapping = graph_analysis.IO.compute_or_load(signatureBasedMapping, "../Data/maps/"+saveStamp, save, sig1, sig2, kernel, maxTopK)
                            
        resN = mappingAccuracy(nodeMapping, forwardMap, topKs)
        print resN
        allAccuracies.append(resN)
 
    accForTop1 =[exp[0][3] for exp in allAccuracies]
    label = ["TopK=1"]
    inText = "Heat Kernel\n#Eigs="+str(allEigs)+"\nTopk=1"
    plt.xQuantityVsAccuracy(exampleG, allNsteps, accForTop1, "TimeSamples" ,label, inText, "*", save=True)
    
    accForTop1 =[exp[1][3] for exp in allAccuracies]
    label = ["TopK=2"]
    inText = "Heat Kernel\n#Eigs="+str(allEigs)+"\nTopk=2"
    plt.xQuantityVsAccuracy(exampleG, allNsteps, accForTop1, "TimeSamples" ,label, inText,"*", save=True)

     
     
     

def accurasyWithSmallestEigens(exampleG=None, nstep=None):
    if exampleG == None:
        exampleG = raw_input("Test module\n Give Graph To Load = ")
    if nstep == None:    
        nstep = int(raw_input("nstep"))
    
    kernel = heatKernel    
    maxTopK = 16
    topKs = [1, 2, 4, 8, maxTopK]
    save = True
    gcc = True
    smallest = True
     
    Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", False, exampleG, gcc)
     
    allEigs = 1000 #maximum Precomputed Eigenvalues
    
    #---- Load MAX Number of Eigenvalues -----
    eigs1Where ="../Data/Spectrums/Smallest_sigma/"+exampleG+"_GCC_"+str(allEigs+1)+"_Eigs"                
    evals1, evecs1 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs1Where, False, Lg, allEigs+1, smallest)
     
    eigs2Where ="../Data/Spectrums/Smallest_sigma/"+exampleG+"_ISO_GCC_"+str(allEigs+1)+"_Eigs"
    evals2, evecs2 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs2Where, False, Lgp, allEigs+1, smallest) 
    
    allAccuracies = []

    eigsToCheck = range(1, 10)  
    eigsToCheck.extend(range(10, allEigs, 10))
    
    #Use Subset of them to do the matching     
    for k in eigsToCheck:
        saveStamp = "HKS/"+exampleG+"_GCC_eigs_"+str(k)+"_guibasTime_"+str(nstep)
              
        timeSamples1 = logTimeScale(evals1, k, nstep)
        sig1 = graph_analysis.IO.compute_or_load(HKS, "../Data/sigs/"+saveStamp, save, evals1, evecs1, k, timeSamples1)
        #--- ISOmorphic case ----          
        timeSamples2 = logTimeScale(evals2, k, nstep)                
        sig2 = graph_analysis.IO.compute_or_load(HKS, "../Data/sigs/"+saveStamp+"_ISO", save, evals2, evecs2, k, timeSamples2)

        nodeMapping = graph_analysis.IO.compute_or_load(signatureBasedMapping, "../Data/maps/"+saveStamp, save, sig1, sig2, kernel, maxTopK)
                            
        resK = mappingAccuracy(nodeMapping, forwardMap, topKs)
        print resK
        allAccuracies.append(resK)
 
    accForAllTopK = []
    for i, topK in enumerate(topKs):
        accForAllTopK.append([exp[i][3] for exp in allAccuracies])
  
    label = ["TopK="+str(topK) for topK in topKs]
    inText = "Heat Kernel\nLogTimeSample |t|="+str(nstep) 
    plt.xQuantityVsAccuracy(exampleG, eigsToCheck, accForAllTopK, "EigenValues", label, inText, "*", save=True)

                                                   
                                                                       
def testKernelWithNoise():    
    DEBUG = True
    PEDANTIC = False
    save = True
        
    kernel = heatKernel
    timeSamples = [0.001, 0.002, 0.003, 0.005, 0.009, 0.010, 0.011, 0.012, 0.22]    
#     allEigs = [1, 10, 20, 100, 500, 1000]
    allEigs = [1, 10]    
    maxTopK = 16
    topKs = [1, 2, 4, 8, maxTopK]
            
    exampleG = raw_input("Heat_kernel module\n Give Graph To Load = ")
#     exampleG = 'karate'    
    gcc = True
            
    Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", save, exampleG, gcc)    
        
    noiseLevels = [0.01, 0.05, 0.1, 0.2]

    allAccuracies = {}
    for noiseLevel in noiseLevels:
        allAccuracies[noiseLevel] = []
        
        LgpNoisy, junk   = graph_analysis.IO.compute_or_load(addUniformRandomEdges, "../Data/Laplacians/Add_Random_Edges/"+exampleG+"_GCC_Noisy_Laplacian_1_AddedEdges_"+str(noiseLevel),                                                
                                                save, Lgp.tolil(), noiseLevel)         
        
        if PEDANTIC:
            assert(isLaplacian(LgpNoisy))
        
        for eigsToUse  in allEigs:

            sig1, sig2          = graph_analysis.IO.compute_or_load(computeISOSignatures, "../Data/sigs/"+exampleG+"_GCC_HKs_eigs_"+str(eigsToUse), 
                                                save, Lg, Lgp, kernel, eigsToUse, timeSamples)        
                        
            noisySig, junk      = graph_analysis.IO.compute_or_load(computeISOSignatures, "../Data/sigs/Add_Random_Edges/"+exampleG+"_GCC_HKs_eigs_"+str(eigsToUse)+"_AddedEdges_"+str(noiseLevel),
                                           save, LgpNoisy, None, kernel, eigsToUse, timeSamples)
            
                    
            nodeMapping, junk   = graph_analysis.IO.compute_or_load(signatureBasedMapping, "../Data/maps/Add_Random_Edges/"+exampleG+"_GCC_HKs_eigs_"+str(eigsToUse)+"_AddedEdges_"+str(noiseLevel)+
            "_Topk_"+str(maxTopK), save, sig1, noisySig, kernel, maxTopK)
        
            allAccuracies[noiseLevel].append(mappingAccuracy(nodeMapping, forwardMap, topKs))
     
    graph_analysis.IO.save_data(exampleG+"Noise_Accuracies_AddEdges", allAccuracies)

        
    for noiseLevel in noiseLevels:
        print "Noise Level = " + str(noiseLevel)    
        for i, exp in  enumerate(allAccuracies[noiseLevel]):
            print "Eigenvalues = " + str(allEigs[i]) 
            print exp
            

    for noiseL in noiseLevels:
        for j, topK in enumerate(topKs):                
            accuraciesAtTopKAtNoise = list()
            for i, exp in  enumerate(allAccuracies[noiseL]):  #exp - experiment based on eigevalue number
                accuraciesAtTopKAtNoise.append( exp[j][3] )
            
            plt.xQuantityVsAccuracy(exampleG, allEigs, accuraciesAtTopKAtNoise, ["TopK="+str(topK)],  
                                       "Heat Kernel\nAdd Edge Noise="+str(noiseL), save=True)
                


#     for i, eig in enumerate(allEigs):
#         for j, noiseL in enumerate(topK):
#             accuraciesForEigAtNoise = [eigAc[3] for eigAc in allAccuracies[i]]
#             plt.topKVsAccuracy(exampleG, topKs, accuraciesForEig, legendlabel=["Heat Kernel\n#Eigs="+str(eig)], save=True)    
             
    sys.exit("Exiting Successfully")
    




     

def accurasyWithMcSherry(exampleG=None, nstep=None):
    if exampleG == None:
        exampleG = raw_input("Test module\n Give Graph To Load = ")
    if nstep == None:    
        nstep = int(raw_input("nstep"))
    
    kernel = heatKernel    
    maxTopK = 16
    topKs = [1, 2, 4, 8, maxTopK]
    save = False
    gcc = True
    smallest = False
     

    Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", save, exampleG, gcc)
    Mg                  = graph_analysis.IO.compute_or_load(graph_analysis.IO.mcSherryFromLaplacian, "../Data/Laplacians/"+exampleG+"_McSher", save, Lg)
    Mgp                 = graph_analysis.IO.compute_or_load(graph_analysis.IO.mcSherryFromLaplacian, "../Data/Laplacians/"+exampleG+"_McSher_Iso", save, Lgp)
        
                            
    allEigs = 500 #maximum Precomputed Eigenvalues
    
    #---- Load MAX Number of Eigenvalues -----
    eigs1Where = "../Data/Spectrums/"+exampleG+"_McSher_"+str(allEigs)+"_Eigs"                
    evals1, evecs1 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs1Where, save, Mg, allEigs, smallest)
     
    eigs2Where = "../Data/Spectrums/"+exampleG+"_McSher_ISO_"+str(allEigs)+"_Eigs"
    evals2, evecs2 = graph_analysis.IO.compute_or_load(computeSpectrum, eigs2Where, save, Mgp, allEigs, smallest) 
    
    ####################################
    # Put inside other functions: 
    #Need to happend For McSherry (where the largest eigenpairs matter the most)
    np.flipud(evecs1)
    np.flipud(evecs2)
    evals1 = evals1[::-1]
    evals2 = evals2[::-1]

    print evals1
    ####################################
        
    allAccuracies = []

    eigsToCheck = range(1, 10)  
    eigsToCheck.extend(range(10, allEigs, 10))
    
    #Use Subset of them to do the matching     
    for k in eigsToCheck:
        saveStamp = "MHKS/"+exampleG+"_GCC_eigs_"+str(k)+"_McSherTime_"+str(nstep)              
        timeSamples1 = McSherryTimeScale(evals1, k, nstep)
        sig1 = graph_analysis.IO.compute_or_load(HKS, "../Data/sigs/"+saveStamp, save, evals1, evecs1, k, timeSamples1)
        print k
        print sig1[0]
        #--- ISOmorphic case ----          
        timeSamples2 = McSherryTimeScale(evals2, k, nstep)                
        sig2 = graph_analysis.IO.compute_or_load(HKS, "../Data/sigs/"+saveStamp+"_ISO", save, evals2, evecs2, k, timeSamples2)

        nodeMapping = graph_analysis.IO.compute_or_load(signatureBasedMapping, "../Data/maps/"+saveStamp, save, sig1, sig2, kernel, maxTopK)
                            
        resK = mappingAccuracy(nodeMapping, forwardMap, topKs)
        print resK
        allAccuracies.append(resK)
 
    accForAllTopK = []
    for i, topK in enumerate(topKs):
        accForAllTopK.append([exp[i][3] for exp in allAccuracies])
  
    label = ["TopK="+str(topK) for topK in topKs]
    inText = "McSherry Heat Kernel\nLogTimeSample |t|="+str(nstep) 
    plt.xQuantityVsAccuracy(exampleG, eigsToCheck, accForAllTopK, "EigenValues", label, inText, "*", save=True)

#     
 

'''
Created on Apr 30, 2014

@author: optas
'''

import sys, itertools
from graph_analysis import IO
import numpy as np 
from role_discovery import graph_kernels as gk
from graph_analysis import basic_tools as myGT
from myUtils import *
import myPloting as myPLT
from scipy.sparse.linalg.eigen.arpack.arpack import eigs


try:
    from graph_tool import collection, topology, spectral, centrality, generation
    from graph_tool import *
    import graph_tool as gt
except:
    "Graph tool not installed."

#     swap poison graphs    
#     save=gcc=True
#     smallest=True
#     nodes = 10000; averDeg = 10; 
# 
#     exampleG = "poisson_"+str(nodes)+"_nodes_"+str(averDeg)+"_avDeg"
#     g, gp, forwardMap, eqClass =  IO.loadData("../Data/Graphs/synthetic/"+exampleG)    
#     
#     for noiseLevel in [0.01, 0.02, 0.03, 0.04, 0.05, 0.10]:
#         gswapped = myGT.swapEdges(g, noiseLevel)
#         gswapped = IO.saveData("../Data/Graphs/synthetic/"+exampleG+"_swappedEdges_"+str(noiseLevel), gswapped)
#         
#     

#make poisson graphs
#     gcc=True
#     nodes = 10000
#     averDeg = 20
#     g, gp, forwardMap = makeISOPoisson(nodes, averDeg, gcc)    
#     gClasses = myGT.structEquivNodes(g)    
#     saveAt="../Data/Graphs/synthetic/poisson_"+str(nodes)+"_nodes_"+str(averDeg)+"_avDeg"        
# 
#     gClasses = myGT.structEquivNodes(g)    
#     assert(myGT.checkEquivNodesStand(gClasses, g))
#       
#     gpClasses = myGT.structEquivNodes(gp)    
#     assert(myGT.checkEquivNodesStand(gpClasses, gp))
#       
#     assert(myGT.checkEquivNodesStandBetweenIso(gClasses, gpClasses, forwardMap))
#     
#     IO.saveData(saveAt, g, gp, forwardMap, gClasses)


#double verify equiv. classes between isomorfphic are ok
#     gClasses = myGT.structEquivNodes(g)    
#     assert(myGT.checkEquivNodesStand(gClasses, g))
#      
#     gpClasses = myGT.structEquivNodes(gp)    
#     assert(myGT.checkEquivNodesStand(gpClasses, gp))
#      
#     assert(myGT.checkEquivNodesStandBetweenIso(gClasses, gpClasses, forwardMap))


##compute equiv.classes and save
#     exampleG = raw_input("example?")
#     gcc = True; save = False;
#     Lg, Lgp, forwardMap = IO.computeOrLoad(IO.makeISOLaplacians, "../Data/Matrices/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", save, exampleG, gcc)
#     g                   = IO.loadData("../Data/Graphs/"+exampleG+"_GCC_purged").next()                
#     gClasses  = myGT.structEquivNodes(g)              
#     IO.saveData(exampleG + "_GCC_equiv_Classes", gClasses)


def testing_Standard_HKS_New_input(exampleG, eigsToLoad, eigsToUse, nstep, smallest=True, plot=False):
    ###PARAMETERS###
    gcc = True
    kernel = gk.heatKernel
    topKs = [1,2]     
    #####-#####-#####    

    graphWhere ="../Data/Graphs/synthetic/"+exampleG           
    g, gp, forwardMap, equivsG = graph_analysis.IO.load_data(graphWhere)
       
    evals1, evecs1, evals2, evecs2 = computeSpectrums(spectral.laplacian(g), spectral.laplacian(gp), eigsToLoad, smallest, exampleG)
    
    res = []    
    if type(eigsToUse) == list: #many eigenvalues to test_myRolx
        for k in eigsToUse:
            timeSample1 = gk.logTimeScale(evals1, k, nstep)      
            sig1 = gk.HKS(evals1, evecs1, k, timeSample1)
                    
            timeSample2 = gk.logTimeScale(evals2, k, nstep)    
            sig2 = gk.HKS(evals2, evecs2, k, timeSample2)
                
            nodeMapping = gk.signatureBasedMapping(sig1, sig2, kernel, topKs[-1])    
            resK = gk.mappingAccuracy(nodeMapping, forwardMap, topKs, equivsG)
            print "Eigenvalues used: "+ str(k)
            print resK
            res.append(resK)
            
        
        accForAllTopK = []
        for i, topK in enumerate(topKs):
            accForAllTopK.append([exp[i][3] for exp in res])
        label = ["TopK="+str(topK) for topK in topKs]
        inText = "Heat Kernel\nLogTimeSample |t|="+str(nstep)
        myPLT.xQuantityVsAccuracy(exampleG, eigsToUse, accForAllTopK, "Eigenvalues", label, inText, save=True)            
    else:
        timeSample1 = gk.logTimeScale(evals1, eigsToUse, nstep)      
        sig1 = gk.HKS(evals1, evecs1, eigsToUse, timeSample1)
                
        timeSample2 = gk.logTimeScale(evals2, eigsToUse, nstep)    
        sig2 = gk.HKS(evals2, evecs2, eigsToUse, timeSample2)
            
        nodeMapping = gk.signatureBasedMapping(sig1, sig2, kernel, topKs[-1])    
        resK = gk.mappingAccuracy(nodeMapping, forwardMap, topKs, equivsG)
        print resK





def testing_Standard_HKS(exampleG, eigsToLoad, eigsToUse, nstep, smallest=True, plot=False):
    ###EXTRA_PARAMETERS###
    save = True
    gcc = True
    kernel = gk.heatKernel
    topKs = [1,2]     
    #####-#####-#####    
    
    Lg, Lgp, forwardMap = graph_analysis.IO.compute_or_load(graph_analysis.IO.makeISOLaplacians, "../Data/Matrices/Laplacians/"+exampleG+"_GCC_&_Isomorfic_Laplacians_1", False, exampleG, gcc)
    equivsG             = graph_analysis.IO.load_data("../Data/Graphs/"+exampleG+"_GCC_equiv_Classes").next()

    if smallest:
        at = "SM"; eigType = "Smallest"
    else:
        at = "LM"; eigType = "Largest"
             
    eigsWhere ="../Data/Spectrums/"+at+"/"+exampleG+"_GCC_"+str(eigsToLoad+1)+"_Eigs"    
    evals1, evecs1 = graph_analysis.IO.compute_or_load(gk.computeSpectrum, eigsWhere, save, Lg, eigsToLoad+1, smallest)
 
    eigsWhere ="../Data/Spectrums/"+at+"/"+exampleG+"_ISO_GCC_"+str(eigsToLoad+1)+"_Eigs"    
    evals2, evecs2 = graph_analysis.IO.compute_or_load(gk.computeSpectrum, eigsWhere, save, Lgp, eigsToLoad+1, smallest)        

    res = []    
    if type(eigsToUse) == list: #many eigenvalues to test_myRolx
        for k in eigsToUse:
            timeSample1 = gk.logTimeScale(evals1, k, nstep)      
            sig1 = gk.HKS(evals1, evecs1, k, timeSample1)
                    
            timeSample2 = gk.logTimeScale(evals2, k, nstep)    
            sig2 = gk.HKS(evals2, evecs2, k, timeSample2)
                
            nodeMapping = gk.signatureBasedMapping(sig1, sig2, kernel, topKs[-1])    
            resK = gk.mappingAccuracy(nodeMapping, forwardMap, topKs, equivsG)
            print "Eigenvalues used: "+ str(k)
            print resK
            res.append(resK)
            
        
        accForAllTopK = []
        for i, topK in enumerate(topKs):
            accForAllTopK.append([exp[i][3] for exp in res])
        label = ["TopK="+str(topK) for topK in topKs]
        inText = "Heat Kernel\nLogTimeSample |t|="+str(nstep)
        xlabel = eigType + " Eigenvalues"
        myPLT.xQuantityVsAccuracy(exampleG, eigsToUse, accForAllTopK, xlabel, label, inText, "-", save)                    
    else:
        timeSample1 = gk.logTimeScale(evals1, eigsToUse, nstep)      
        sig1 = gk.HKS(evals1, evecs1, eigsToUse, timeSample1)
                
        timeSample2 = gk.logTimeScale(evals2, eigsToUse, nstep)    
        sig2 = gk.HKS(evals2, evecs2, eigsToUse, timeSample2)
            
        nodeMapping = gk.signatureBasedMapping(sig1, sig2, kernel, topKs[-1])    
        resK = gk.mappingAccuracy(nodeMapping, forwardMap, topKs, equivsG)
        print resK




def computeSpectrums(M1, M2, eigsToUse, smallest, matrixName=None, M1eigsFile=None, M2eigsFile=None):
                         
    print "Eigenvalues that will be computed:" + str(eigsToUse+1)
    save = True


    if M1eigsFile == None:
        assert(M2eigsFile == None and matrixName != None)
        if smallest:    
            M1eigsFile = "../Data/Spectrums/SM/"+matrixName+"_Eigs_"+str(eigsToUse+1)
            M2eigsFile = "../Data/Spectrums/SM/"+matrixName+"_Eigs_"+str(eigsToUse+1)+"_ISO"
        else:
            M1eigsFile = "../Data/Spectrums/LM/"+matrixName+"_Eigs_"+str(eigsToUse+1)
            M2eigsFile = "../Data/Spectrums/LM/"+matrixName+"_Eigs_"+str(eigsToUse+1)+"_ISO"
    
    evals1, evecs1 = graph_analysis.IO.compute_or_load(gk.computeSpectrum, M1eigsFile, save, M1, eigsToUse+1, smallest)
    evals2, evecs2 = graph_analysis.IO.compute_or_load(gk.computeSpectrum, M2eigsFile, save, M2, eigsToUse+1, smallest)
    
    return evals1, evecs1, evals2, evecs2
    
        
                                  
def makeISOPoisson(numberOfNodes, averDegree, gcc, save=False, saveAt=None):
    def degreeSample():
        return np.random.poisson(averDegree)
     
    g = generation.random_graph(numberOfNodes, degreeSample, directed=False)
    if gcc:
        l = topology.label_largest_component(g)  #Keep Largest Connected Component
        g.set_vertex_filter(l); g.purge_vertices()
        
    Ag = spectral.adjacency(g).todense()  
    gp, forwardMap = graph_analysis.IO.createIsoGraph(Ag)
         
    if save:
        assert(saveAt != None)
        graph_analysis.IO.save_data(saveAt, g, gp, forwardMap)
    
    return g, gp, forwardMap
    



#TODO: see and put in right places
# 
# 1. calculate the signature first
# sigsPCA = PCA(sig) 
# print "Fraction of Variance in the first three PCA dimensions: " +str(sum(sigsPCA.fracs[0:3]))    
# myPalet = palets("intToColor")
# nodesColored = [myPalet[i] for i in bm.a]
# myPLT.plot_PCA_Res(sigsPCA, dim=2, legendLabel=None, colorMap=[nodesColored, palets("intToColor")], save=False, saveAt=None)
    
# 2.
# graphFile   = "../Data/Graphs/web_snapShot/web.gml"
# exampleName = "webSnapShot"
# g           = IO.loadGraph(graphFile, gcc=True, removeSL=True, removePE=True)
# kcore       = kcore_decomposition(g, deg="total")
# numCores    = max(kcore.a) + 1
# print "Number of Cores = ", numCores
# # core_hist = stats.vertex_hist(g, kcore)
# # print core_hist
# # plt.title("K-Shell Distribution of Nodes, WebSnapShot")
# # plt.xlabel("K-shell")
# # plt.hist(kcore.a.tolist(), bins=maxCore, normed=True)
# # plt.show()
# # print sum(core_hist[0][1:4]) / float(g.num_vertices()) #fraction of 3 first roles among all roles

#3
# from sklearn.neighbors import KNeighborsClassifier
# labels = kcore.a
# knn = KNeighborsClassifier()
# knn.fit(sig[0:18000,:], labels[0:18000])
# KNeighborsClassifier(algorithm='auto', leaf_size=30, metric='minkowski', n_neighbors=1, p=4, weights='uniform')
# yP = knn.predict(sig[18001:])
# print sum (yP)
# print len(yP)
# print sum(yP == labels[18001:]) / float(len(yP))






if __name__ == "__main__":

    pass