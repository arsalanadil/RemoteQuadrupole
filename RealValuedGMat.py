# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 21:42:40 2017

Generating the Real Valued Variance Matrix using the Covariance and Relation Matrices obtained from COVandREL.py

@author: arsalanadil
"""
import numpy as np
import COVandREL as CV


def CovRel(l, clusters, n):
    return CV.Baap(l,clusters,n)


randClusters = np.load("rand20Clusters.npy")#cluster list
#v1 = np.array([2.94,0.75, 0.172])
#v2 = np.array([2.09,0.3,0.213])
#v3 = np.array([0.09,1.1,0.57])
#CovAndRel = np.load("CovRel2.npy")
CovAndRel = CovRel(3, randClusters, 10)#l,clusters,n

varMatrix = CovAndRel[0]
relMatrix = CovAndRel[1]

GammaXX = 1/2 * (np.real(varMatrix) + np.real(relMatrix))
GammaXY = 1/2 * ((np.imag(relMatrix) - np.imag(varMatrix)))
GammaYX = 1/2 * ((np.imag(relMatrix) + np.imag(varMatrix)))
GammaYY = -1/2 *((np.real(relMatrix) - np.real(varMatrix)))
        
gammaR = np.bmat([[GammaXX , GammaXY],[GammaYX, GammaYY]])


a = np.linalg.eigvals(gammaR)#check for positive definiteness
for i in range(0, a.shape[0]):
    if(a[i]<-0.001 or np.abs(a[i].imag)>0.001):
        print("GammaR. Tere chudday. Not Positive Definite")
    
b = gammaR - np.transpose(gammaR)#check for symmetricity
for i in range(0,b.shape[0]):
    for j in range(0,b.shape[1]):
        if(b[i,j] > 0.000001):
            print("GammaR. BC tere tou chudday. Not Symmetric")
print("N=10 and Lmax=3 with increased precision")            
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