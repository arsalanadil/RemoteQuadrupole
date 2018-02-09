#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 19:10:50 2017
Testing different parallelism methods.....

@author: arsalanadil
"""

from multiprocessing import Pool
import RemoteQuad as Remote
import time
import numpy as np
from functools import partial
#import RemoteQuad_Harmonic as Remote2
import matplotlib.pyplot as plt
def f(v1,v2):
    return Remote.wignerRotation(v1,v2)[0][4][4]
    
#
"""
ta = time.time()    
print([f(i) for i in range(1,N)])
tb = time.time()
print("Time taken in second process:", tb-ta)
"""
"""
N = 10#number of processes
clusters = np.load("VoxList.npy")
if __name__ == '__main__':
    t0 = time.time()
    with Pool() as p:
        for j in range(0,N):
            print("CLUSTER NUMBER,",j)
            print()
            print()
            print(p.starmap(f, [(clusters[i], clusters[j]) for i in range(j,N)] ))
    t1 = time.time()
    print("TIme Taken:", t1-t0)
"""
#def g(v1,v2):
#    return Remote2.CovRel(v1,v2)[0]

"""
----------------------
The following code can be used to generate the variance matrix by using multiple processors.
"""
N = 30#number of processes
clusters = np.load("VoxList.npy")
#clusters = np.load("testClusters.npy")

varMatOld = np.zeros((N,N)).astype(complex)
varMatNew = np.zeros((N,N)).astype(complex)
"""
zlist = np.random.uniform(0.1,1,N)#sort this list, thats the one we want
z_testList = np.sort(zlist)
phi_testList = np.random.uniform(0,np.pi/2,N)
theta_testList = np.random.uniform(0,np.pi/2,N)
"""
#clusters = np.vstack((z_testList, theta_testList, phi_testList)).T
import sys                      
if __name__ == '__main__':
    t0 = time.time()
    with Pool(3) as p:
        temp = (p.starmap(f, [(clusters[i], clusters[j]) for i in range(0,N) for j in range(i,N)] ))
        #temp = (p.starmap(g, [(i,j) for i in range(0,N) for j in range(i,N)] ))
    t1 = time.time()
    print("TIme Taken:", t1-t0)
    k = 0
    for i in range(0,N):
        for j in range(i,N):
            varMatOld[i][j] = temp[k]
            varMatOld[j][i] = np.conjugate(temp[k])
            k+=1
    p.close()
    p.join()
    p.terminate()
    sys.exit()
    
        
    """    
    t0 = time.time()
    with Pool(8) as p:
        temp2 = (p.starmap(g, [(clusters[i], clusters[j]) for i in range(0,N) for j in range(i,N)] ))
        #temp = (p.starmap(g, [(i,j) for i in range(0,N) for j in range(i,N)] ))
    t1 = time.time()
    print("TIme Taken:", t1-t0)
    k = 0
    for i in range(0,N):
        for j in range(i,N):
            varMatNew[i][j] = temp2[k]
            varMatNew[j][i] = np.conjugate(temp2[k])
            k+=1
    
    vMatNewReal = np.abs(varMatNew + np.conjugate(varMatNew))/2
    vMatNewIm =  np.abs(varMatNew - np.conjugate(varMatNew))/2

    #vMat = (varMatNew) - np.conjugate(varMatOld)
    #vMat = np.abs(np.imag(varMatOld)) - vMatNewIm
    vMat = np.abs(np.real(varMatOld)) - vMatNewReal
    for i in range(0,N):
        for j in range(0,N):
            if(np.abs(vMat[i][j]) > 0.1):
                print("Masla at", i,j,"Val:",vMat[i][j])
    """
    
"""
temp = [Remote2.CovRel([1,0.3,0.9],[1,i,0.4])[0] for i in np.arange(0,np.pi/2,0.2)]
plt.plot((temp + np.conjugate(temp))/2)

plt.plot([Remote.wignerRotation([1,0.3,0.9],[1,i,0.4])[0][0][0] for i in np.arange(0,np.pi/2,0.2)])
"""
"""
t0 = time.time()
[f(clusters[i], clusters[j]) for i in range(0,N) for j in range(i,N)]
t1 = time.time()
print("TIME USING OLD CODE", t1-t0)
"""