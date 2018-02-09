t#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 16:22:00 2017

@author: arsalanadil
"""

import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import RemoteQuad_Harmonic as Remote
#import AngPS as Local
import CrossCor as Cross
import matplotlib.pyplot as plt
import time


def LocalCMB(l):    
    length = l**2 + 2*l -3
    store = np.zeros(l-1)
    clmat = np.zeros((length,length)).astype(complex)
    clmat2 = np.zeros((length,length)).astype(complex)
    
    Cls = np.load("locallmax7.npy")#stored values for Cls upto lmax=7. 
    
    for i in range(2, l+1):
        #store[i-2] = Local.Integral(i)
        store[i-2] = Cls[i-2]
        #print(store[i-2], i)
        
    clmat[0][0] = store[0]
    clmat2[0][4] = store[0]
    p=0 
    n=0
    y=0
    for j in range(2,l+1):
        i = 0
        for x in range(n, n + 2*j + 1):
            #print(y+1)
            clmat[y+1][y+1] = store[p]
            clmat2[y+1][n + 2*j - i] = -(-1)**y * store[p]            
            
            y = x
            i += 1
            
        n = y+1
        p+= 1
    
    clmat2[1][4] = 0 #the algorithm only gets this co-ordinate wrong! So fixing it manually
        
    return clmat, clmat2
    
def CrossCMB(clusters,l , n):
    length = l**2 + 2*l -3 
    CrossCovMat = np.zeros((n, length)).astype(complex)
    CrossRelMat = np.zeros((n, length)).astype(complex)
    
    for j in range(0,n):
        print(" Cluster Numnber ", j)        
        crosscorel = Cross.wignerRotation(clusters[j,:], 2)[0]         
        Row = crosscorel[4][:]#we want only one row for the covariance matrix
        Row2 = crosscorel[0][:]#and only one row for the relation matrix
        
        for i in range(3,l+1):
            crosscorel2 = Cross.wignerRotation(clusters[j,:], i)[0]            
            temp = crosscorel2[4][:]
            temp2 = crosscorel2[0][:]
            
            Row = np.concatenate((Row, temp))
            Row2 = np.concatenate((Row2, temp2))
            
        CrossCovMat[j] = np.conjugate(Row)
        #CrossCovMat[j] = Row
        CrossRelMat[j] = Row2
        
    return CrossCovMat, CrossRelMat

def RemoteCMB(clusters, n):

    varMatrix = np.zeros((n,n)).astype(complex)
    relMatrix = np.zeros((n,n)).astype(complex)
    
    for i in range(n):
        for j in range(i+1):
            #print("For Remote: Cluster number", i, "and",j)
            temp = Remote.CovRel(clusters[i,:],clusters[j,:])
            print("For Remote: Cluster number", i, "and",j,".Covariance",temp[0])
            #print("wignerRotation() returned:", temp[0])        
            varMatrix[i][j] = temp[0]
            varMatrix[j][i] = np.conjugate(temp[0])
            relMatrix[i][j] = temp[1]
            relMatrix[j][i] = temp[1]
    
    return np.conjugate(varMatrix), relMatrix

 
   
def Baap(l, clusters, n):
    
    Local = LocalCMB(l)
    #Cross = np.load("Cross2ClustersR1.npy")# or call CrossCMB(clusters,l,n)    
    Cross = CrossCMB(clusters,l,n)
    Remote = RemoteCMB(clusters,n)    
    
    A = Remote[0]
    B = Cross[0]  
    C = np.transpose(np.conjugate(B))
    D = Local[0]
    
    #print("Remote is:", A.shape, "Cross is:", B.shape, "Local is:", D.shape)
    
    Covariance = np.bmat([[A , B],[C, D]])
    
    A1 = Remote[1]
    B1 = Cross[1]
    C1 = np.transpose(B1)
    D1 = Local[1]
    
    Relation = np.bmat([[A1, B1],[C1, D1]])
    
    return Covariance, Relation

def powerSpec(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**3
"""
randClusters = np.load("rand20Clusters.npy")#cluster list
n = 10#number of clusters we want to consider
t0 = time.time()
crossCode = CrossCMB(randClusters, 5, n)
t1 = time.time()

print("TIme taken for remote for",n,"clusters:", (t1-t0))
"""

"""
randClusters = np.load("rand20Clusters.npy")#cluster list
n = 4#number of clusters we want to consider
l = 2#lmax     
CovRel = Baap(l, randClusters,n)
"""
#cov = Baap(l, np.array([[0,0,0],[0,0,0]]), 1)

#a = LocalCMB(5)
#print(a)   
"""
q = CrossCMB(randClusters, l, n) 
print(q)
np.save("Cross10ClustersL8", q)
"""