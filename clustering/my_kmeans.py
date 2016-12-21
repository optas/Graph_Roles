'''
Created on Feb 10, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''



estimator = KMeans(init='k-means++', n_clusters=n_digits, n_init=10)
estimator.fit(data)