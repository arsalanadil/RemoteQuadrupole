#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 08:56:27 2017

This routine is meant to test the result of real valued covariance matrix generated using "RealValuedGmat.py". Currently, the problem I am facing
is that for lmax>2 for N ~ 20, the covariance matrix has (~1) negative eigenvalue, which is quite perplexing since if I reduced N~10 that problem
goes away. 
Thanks to optimizations, we can now carry out these tests fairly quickly which is really helpful. 


@author: arsalanadil
"""

import numpy as np
import COVandREL_Opt as CV
#import COVandREL_Harmonic as CV
import time 

def CovRel(l, clusters, n):
    return CV.Baap(l,clusters,n)


def GMat(lmax, clusters,n): 
    CovAndRel = CovRel(lmax, clusters, n)#l,clusters,n

    varMatrix = CovAndRel[0]
    relMatrix = CovAndRel[1]
    
    GammaXX = 1/2 * (np.real(varMatrix) + np.real(relMatrix))
    GammaXY = 1/2 * ((np.imag(relMatrix) - np.imag(varMatrix)))
    GammaYX = 1/2 * ((np.imag(relMatrix) + np.imag(varMatrix)))
    GammaYY = -1/2 *((np.real(relMatrix) - np.real(varMatrix)))
            
    gammaR = np.bmat([[GammaXX , GammaXY],[GammaYX, GammaYY]])
    
    return gammaR



randClusters = np.load("rand100Clusters.npy")#cluster list
n = 20
lmax = 3
"""
for k in range(0,n):
    for j in range(0,n):
        print("Cluster", k, "and", (j))
        v1 = randClusters[k]
        v2 = randClusters[j]
        gammaR = GMat(lmax, np.array([v1,v2]),2)
        
        a = np.linalg.eigvals(gammaR)#check for positive definiteness
        for i in range(0, gammaR.shape[0]):
            if(a[i]<-0.001 or np.abs(a[i].imag)>0.001):
                print("GammaR. Tere chudday. Not Positive Definite")
        
        
        c = np.sort(np.diag(gammaR))[(gammaR.shape[0]-1)]
        for i in range(0,gammaR.shape[0]):#check that largest value occurs on diagonal
            for j in range(0,gammaR.shape[0]):
                if(gammaR[i,j] > c and i != j):
                    print("Ye masla hai bc,", "(i,j)",i,j,"value",gammaR[i,j])
            
        b = gammaR - np.transpose(gammaR)#check for symmetricity
        for i in range(0,b.shape[0]):
            for j in range(0,b.shape[1]):
                if(b[i,j] > 0.000001):
                    print("GammaR. BC tere tou chudday. Not Symmetric")
"""        
                      
#v1 = randClusters[i]
#print("Cluster location:", v1)
#clusters = np.vstack([randClusters[j,:] for j in range(10,20)])
#v1 = np.array([2.94,0.75, 0.172])
#v2 = np.array([2.09,0.3,0.213])
#v3 = np.array([0.09,1.1,0.57])
#CovAndRel = np.load("CovRel2.npy")
#v1 = randClusters[19]
#v2 = randClusters[0]
t0 = time.time()
gammaR = GMat(lmax, randClusters,n)
t1 = time.time()
print("TIME TAKEN IN GENERATING GAMMAR for Lmax =5 and N=50:", (t1-t0))

a = np.linalg.eigvals(gammaR)#check for positive definiteness
for i in range(0, a.shape[0]):
    if(a[i]<-0.001 or np.abs(a[i].imag)>0.001):
        print("GammaR. Tere chudday. Not Positive Definite")


c = np.sort(np.diag(gammaR))[(gammaR.shape[0]-1)]
for i in range(0,gammaR.shape[0]):#check that largest value occurs on diagonal
    for j in range(0,gammaR.shape[0]):
        if(gammaR[i,j] > c and i != j):
            print("Ye masla hai bc,", "(i,j)",i,j,"value",gammaR[i,j])
    
b = gammaR - np.transpose(gammaR)#check for symmetricity
for i in range(0,b.shape[0]):
    for j in range(0,b.shape[1]):
        if(b[i,j] > 0.000001):
            print("GammaR. BC tere tou chudday. Not Symmetric")

