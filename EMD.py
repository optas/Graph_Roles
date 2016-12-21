'''
Created on Dec 31, 2014

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2014, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

# 
# class MyClass(object):
#     '''
#     classdocs
#     '''
# 
# 
#     def __init__(self, params):
#         '''
#         Constructor
#         '''
#         pass
#     
#     def Static lala(self):
        
try:
    from cv2 import *
except:
    pass
 
import numpy as np


def clear_repeated_converged_values2(M, thres):
    '''   M - (m x n) np.array of variables x observations.
    '''
    print M
    m, n = M.shape
    assert(n > 2)
    res  = []
    for row in xrange(m):
        stopped = False
        diffOld = abs(M[row,1] - M[row,0])               
        for col in xrange(1, n-1):
            diffNew = abs(M[row, col] - M[row, col+1])
#             print diffOld, diffNew
#             raw_input()        
#             if abs(diffNew - diffOld) < thres * diffOld:
            print abs(diffNew - diffOld)
            raw_input()
            if abs(diffNew - diffOld) < 0.001:
                res.append(np.array(M[row, 0:col]))
                stopped = True
                print col
                break;
            else:
                diffOld = diffNew                        
        if stopped == False:
            res.append(np.array(M[row,:]))
    return res


def calc_emd(x, y, wx=None, wy=None, distance=None):
    xs = np.size(x); 
    ys = np.size(y);

    # Set the default weight of each coordinate to uniform.
    if wx == None:
        wx = np.ones((xs, 1)) * (1. / xs)
    if wy == None:
        wy = np.ones((ys, 1)) * (1. / ys)
    
    
    x = np.hstack((wx, x.reshape(xs, 1)))
    y = np.hstack((wy, y.reshape(ys, 1)))
    

    # Convert from numpy array to CV_32FC1 Mat
    x64 = cv.fromarray(x)
    x32 = cv.CreateMat(x64.rows, x64.cols, cv.CV_32FC1)
    cv.Convert(x64, x32)

    y64 = cv.fromarray(y)
    y32 = cv.CreateMat(y64.rows, y64.cols, cv.CV_32FC1)
    cv.Convert(y64, y32)

    # Calculate Earth Mover's
    return cv.CalcEMD2(x32, y32, cv.CV_DIST_L2)