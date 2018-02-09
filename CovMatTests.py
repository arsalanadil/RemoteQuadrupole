#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:35:42 2017
Testing CovandMat

@author: arsalanadil
"""


import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import VarMatrix_Generalized as Remote
import LocalPS as Local
#import VarMatrix_Gen_Local as Cross
import matplotlib.pyplot as plt

def LocalCMB(l):    
    length = l**2 + 2*l -3
    store = np.zeros(l-1)
    clmat = np.zeros((length,length)).astype(complex)
    clmat2 = np.zeros((length,length)).astype(complex)
    
    for i in range(2, l+1):
        store[i-2] = Local.Integral(i,powerSpec)
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
        crosscorel = Remote.wignerRotation(clusters[j,:],np.array([0,0,0]))[0]         
        Row = crosscorel[4][:]#we want only one row for the covariance matrix
        Row2 = crosscorel[0][:]#and only one row for the relation matrix
        
            
        CrossCovMat[j] = np.conjugate(Row)
        CrossRelMat[j] = Row2
        
    return CrossCovMat, CrossRelMat

def RemoteCMB(clusters, n):

    varMatrix = np.zeros((n,n)).astype(complex)
    relMatrix = np.zeros((n,n)).astype(complex)
    
    for i in range(n):
        for j in range(i+1):
            temp = Remote.wignerRotation(clusters[i,:],clusters[j,:])[0]
            #print("wignerRotation() returned:", temp[0])        
            varMatrix[i][j] = temp[0][0]
            varMatrix[j][i] = np.conjugate(temp[0][0])
            relMatrix[i][j] = temp[0][4]
            relMatrix[j][i] = temp[0][4]
    
    return varMatrix, relMatrix

 
   
def Baap(l, clusters, n):
    
    Local = LocalCMB(l)
    #Cross = np.load("Cross2ClustersR1.npy")# or call CrossCMB(clusters,l,n)    
    Cross = CrossCMB(clusters,l,n)
    Remote = RemoteCMB(clusters,n)    
    
    A = Remote[0]
    B = Cross[0]  
    C = np.transpose(np.conjugate(B))
    D = Local[0]
    
    print("Remote is:", A.shape, "Cross is:", B.shape, "Local is:", D.shape)
    
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
n = 2#number of clusters we want to consider
l = 3#lmax     
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