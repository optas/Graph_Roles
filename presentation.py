
'''
Created on Jul 1, 2014

Contact: pachlioptas@gmail.com

Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''


from myUtils import *
import graph_tool as gt
from graph_tool import topology, spectral
from role_discovery import graph_kernels as gk
from graph_analysis import IO
import sys
import os
import myPloting as myPlt
from role_discovery import graphs_with_roles as roles
from __builtin__ import raw_input
from clustering.distances import jaccard_vectorial
from scipy.sparse import linalg, identity, diags






    
        
def plotResBest(graphName, allMethods, all_res, trueTaxaNum, comment=None):  #rewrite
    vmea = []; amut = []; aran = [];
        
    for r in all_res:
        if r[0] != 1:           #not a single cluster
            vmea.append(r[3]); aran.append(r[4]); amut.append(r[5])
        else:
            vmea.append(0); aran.append(0); amut.append(0);

    vmea = relative_to_worst(vmea)
    amut = relative_to_worst(amut)
    aran = relative_to_worst(aran)

    x     = allMethods
    y     = [vmea, amut, aran]
    
    
    legendLabel=["V","M","R"]; save = True; saveAt = "../OutPut/" + graphName + comment  
    xAxisLabel=None; yAxisLabel='Improvement from worst (percent)', 
    title="Clustering Efficiency: " + graphName; commentText=None; stamp="*" ; xTick=True; yNorm = False
    
    myPlt.x_vs_y(graphName, x, y, xAxisLabel, yAxisLabel, title, legendLabel, commentText, stamp, xTick, yNorm, save, saveAt)