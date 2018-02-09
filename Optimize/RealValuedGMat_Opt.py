# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 21:42:40 2017

Generating the Real Valued Variance Matrix using the Covariance and Relation Matrices obtained from COVandREL.py

@author: arsalanadil
"""
import numpy as np
import COVandREL_Opt as CV
import time 

def CovRel(l, clusters, n):
    return CV.Baap(l,clusters,n)

#randClusters = np.load("rand20Clusters.npy")#cluster list
randClusters = np.load("rand30Clusters.npy")
#v1 = randClusters[i]
#print("Cluster location:", v1)
#clusters = np.vstack([randClusters[j,:] for j in range(470,500)])
#1 = np.array([2.94,0.75, 0.172])
#v2 = np.array([2.09,0.3,0.213])
#v3 = np.array([0.09,1.1,0.57])
#CovAndRel = np.load("CovRel2.npy")

v1 = randClusters[0]
v2 = randClusters[9]
t0 = time.time()
#CovAndRel = CovRel(3, np.array([v1,v2]), 2)#l,clusters,n
CovAndRel = CovRel(3, randClusters,20)
#last run at n=17
varMatrix = CovAndRel[0]
relMatrix = CovAndRel[1]

GammaXX = 1/2 * (np.real(varMatrix) + np.real(relMatrix))
GammaXY = 1/2 * ((np.imag(relMatrix) - np.imag(varMatrix)))
GammaYX = 1/2 * ((np.imag(relMatrix) + np.imag(varMatrix)))
GammaYY = -1/2 *((np.real(relMatrix) - np.real(varMatrix)))
        
gammaR = np.bmat([[GammaXX , GammaXY],[GammaYX, GammaYY]])
t1 = time.time()
print("TIME TAKEN IN GENERATING GAMMAR:", (t1-t0))

a = np.linalg.eigvals(gammaR)#check for positive definiteness
for i in range(0, a.shape[0]):
    if(a[i]<-0.001 or np.abs(a[i].imag)>0.001):
        print("GammaR. Tere chudday. Not Positive Definite. Eigval:", a[i])


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
n = 2
for j in range(0,n):
    for i in range(n,14):
        temp = np.abs(varMatrix[j,i])<= np.sqrt(np.abs(varMatrix[i,i] * varMatrix[j,j]))
        if(temp == False):
            print("Equivalence relation not true!")
"""
     
       
"""       
def covMat():#trying to reconstruct covariance matrix
     cov = GammaXX + GammaYY + 1j*(GammaYX - GammaXY)
     a = np.linalg.eigvals(cov)#check for positive definiteness
     for i in range(0, a.shape[0]):
        if(a[i]<-0.001 or np.abs(a[i].imag)>0.001):
            print("CovMatrix. Tere chudday. Not Positive Definite")
        
     b = cov - np.conjugate(np.transpose(cov))#check for Hermicity
     for i in range(0,b.shape[0]):
        for j in range(0,b.shape[1]):
            if(np.abs(b[i,j]) > 0.000001):
                print("CovMatrix. BC tere tou chudday. Not Hermitian")
     return cov
 
def relMat():#reconstructing the relation matrix
    rel = GammaXX - GammaYY + 1j*(GammaYX + GammaXY)
    
    b = rel - (np.transpose(rel))#check for Symmetricity
    for i in range(0,b.shape[0]):
        for j in range(0,b.shape[1]):
            if(np.abs(b[i,j]) > 0.000001):
                print("RelMatrix. BC tere tou chudday. Not Symmetric")
    return rel
"""