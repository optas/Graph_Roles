'''
Created on Feb 27, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

class nodeClusteringResults(object):
    '''
    Keep the results of different experiments measuring graph node-clusterings.
    '''
        
    def __init__(self, methodName, graphName, methodParams=None):
        '''
        Constructor
        '''
        self.methodName   = methodName
        self.graphName    = graphName
        self.methodParams = methodParams
        self.clusterings  = {}               # Keeps (methodName - param) -> node clustering produced
        self.evaluations  = {}                # Keep the statistics of different evaluations 
        
    
    def add_clustering_labels(self, clusteringMethod, clusteringParams, producedLabels):        
        self.clusterings[(clusteringMethod, clusteringParams)] = producedLabels

    def get_clustering_labels(self, clusteringMethod, clusteringParams):
        return self.clusterings[(clusteringMethod, clusteringParams)]
    
    def add_clustering_evaluation(self, clusteringMethod, clusteringParams,  measureName, measureResult=None):
        if type(measureName) == str:
            self.evaluations[(clusteringMethod, clusteringParams, measureName)] = measureResult
        else:
            assert(type(measureName) == dict)
            for key in measureName.keys():            
                self.evaluations[(clusteringMethod, clusteringParams, key)] = measureName[key]

    
    def get_clustering_evaluation(self, clusteringMethod, clusteringParams,  measureName):
        return self.evaluations[(clusteringMethod, clusteringParams, measureName)]
    
    