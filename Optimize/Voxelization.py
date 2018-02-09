#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 14:50:39 2017

This script is meant to determine the Voxelization parameters.

@author: arsalanadil
"""
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import RemoteQuad as Remote
import matplotlib.pyplot as plt
import time
#import RemoteQuad_Harmonic as Remote2

def RemoteCMB(clusters, n):

    varMatrix = np.zeros((n,n)).astype(complex)
    relMatrix = np.zeros((n,n)).astype(complex)
    
    for i in range(n):
        print(i)
        for j in range(i+1):
            #print("For Remote: Cluster number", i, "and",j)
            temp = Remote.wignerRotation(clusters[i,:],clusters[j,:])[0]
            #print("For Remote: Cluster number", i, "and",j,".Covariance",temp[0][0])
            #print("wignerRotation() returned:", temp[0])        
            varMatrix[i][j] = temp[0][0]
            varMatrix[j][i] = np.conjugate(temp[0][0])
            relMatrix[i][j] = temp[0][4]
            relMatrix[j][i] = temp[0][4]
    
    return varMatrix, relMatrix

"""
The following piece of code gives you the voxelization in z (using the old code - VarMatrixGeneralized.py)
"""
"""
zVoxel = np.zeros(20)#don't know how much we'll need. keeping 20 to be safe.
n = 200
zvals = np.linspace(0,2,n)
clusters = np.array([[zvals[i],0,0] for i in range(0,n)])
#varMat = RemoteCMB(clusters,n)[0]
varMat = np.load("VoxelizeCovarianceMat.npy")

cut = 0.96
m = 0
k = 0
i = 0
print("RESULT USING OLD CODE:")
while(i < n):
    for j in range(k,n):
        corel = varMat[i][j]/np.sqrt(varMat[i][i] * varMat[j][j])
        
        if(i != j and corel.real < cut):
            #print(i,j)
            #cut = varMat[i][j]/np.sqrt(varMat[i][i] * varMat[j][j])
            print("Correlation Coefficient:", corel)
            print("z",j,clusters[j])
            zVoxel[m] = zvals[j]#store this value since we need it to figure out theta and phi voxelization.
            k = j
            i = k
            m+=1
            break        

"""








"""
------------------------------------------------------------------
Use the following if you want to use/compare results with the new code in Harmonic Space:
"""
            
"""
def RemoteHarm(clusters,n):
    varMatrix = np.zeros((n,n)).astype(complex)
    for i in range(n):
        print(i)
        for j in range(i+1):
            #print("For Remote: Cluster number", i, "and",j)
            temp = Remote2.UnrotCovRel(0,clusters[i][0],clusters[j][0])[0]#for things on the z-axis we only care about UnrotCovRel
            #print("For Remote: Cluster number", i, "and",j,".Covariance",temp)
            varMatrix[i][j] = temp
            varMatrix[j][i] = np.conjugate(temp)
    return varMatrix

def RempoteHarmComp(clusters,n):
    varMatrix = np.zeros((n,n)).astype(complex)
    for i in range(n):
        print(i)
        for j in range(i+1):
            #print("For Remote: Cluster number", i, "and",j)
            temp = Remote2.CovRel(clusters[i],clusters[j])[0]
            #print("For Remote: Cluster number", i, "and",j,".Covariance",temp)
            varMatrix[i][j] = temp
            varMatrix[j][i] = np.conjugate(temp)
    return varMatrix
"""
"""
n = 10
zvals = np.linspace(0,0.5,n)
cut = 0.95

varMat2 = RemoteHarm(zvals,n)
print("RESULT USING Harmonic Space CODE:")

k = 0
i = 0
while(i < n):
    for j in range(k,n):
        corel = varMat2[i][j]/np.sqrt(varMat2[i][i] * varMat2[j][j])
        
        if(corel.real < cut):
            #print(i,j)
            #cut = varMat[i][j]/np.sqrt(varMat[i][i] * varMat[j][j])
            print("Correlation Coefficient:", corel)
            print("z",j,zvals[j])
            k = j
            i = k
            break
"""



"""
---------------------------------------------------------------------------
The following piece is for determining Nside for each "z-shell". 
"""
"""
zVox = np.load("zVox.npy")
zVox = np.concatenate((zVox,[2]))
Nstore = np.zeros(len(zVox)+1)
j = 0
for i in range(0,len(zVox)):
    print('\n')
    print("At shell number:", i, "corresponding to z=",zVox[i])
    corel = 1
    Nside = 32
    while(corel > 0.96):
        z = zVox[i]
        Np = 12 * Nside**2
        theta = np.sqrt(4*np.pi/Np)
        corel = np.abs(
            Remote.wignerRotation([z,0,0],[z,theta,0])[0][0][0]/np.sqrt(Remote.wignerRotation([z,0,0],[z,0,0])[0][0][0] 
            * Remote.wignerRotation([z,theta,0],[z,theta,0])[0][0][0]))
        corel2 = np.abs(
            Remote.wignerRotation([z,0,theta],[z,theta,theta/np.sin(theta)])[0][0][0]/np.sqrt(Remote.wignerRotation([z,theta,theta/np.sin(theta)],[z,theta,theta/np.sin(theta)])[0][0][0] 
            * Remote.wignerRotation([z,theta,0],[z,theta,0])[0][0][0]))
        print(" Theta Corelation:", corel, "Phi Corel:",corel2, "Nside:",Nside, "delta Theta:", theta, "delta phi:", theta/np.sin(theta))
        if(corel < 0.96 or corel2 <0.96):
            Nstore[j] = Nside * 2#stores the Nside that I want to use for each zshell (for ease of checking in the end)
            j+=1
        Nside = Nside/2#Nside needs to be in powers of two...

    #print("Corel:", corel,"theta:",theta, "Corresponding to Nside=",Nside)

"""




"""
------------------------------------------------------------------
The following piece is for verifying that \delta theta is uniform for any given z-shell
"""
"""
Nside = 4
z = 0.070351758794
Np = 12 * Nside**2
tempTheta = np.sqrt(4*np.pi/Np)
delTheta = 0
theta = tempTheta
while(delTheta < np.pi):
    corel = np.abs(
    Remote.wignerRotation([z,0,delTheta],[z,0,theta])[0][0][0]/np.sqrt(Remote.wignerRotation([z,0,delTheta],[z,0,delTheta])[0][0][0] 
    * Remote.wignerRotation([z,0,theta],[z,0,theta])[0][0][0]))
    
    #At each theta value, we want delta phi = delta t*heta/ sin(theta) and this should return the same correlation as above.
    
    delPhi = tempTheta/np.sin(theta)
    corel2 = np.abs(
            Remote.wignerRotation([z,0,theta],[z,delPhi,theta])[0][0][0]/np.sqrt(Remote.wignerRotation([z,delPhi,theta],
                                  [z,delPhi,theta])[0][0][0] 
            * Remote.wignerRotation([z,0,theta],[z,0,theta])[0][0][0]))
    print("theta:",theta,"theta2:",delTheta,"Corelation",corel, "Phi Correlation:", corel2, "delPhi:", delPhi)
    print("Number of Phi Voxels at this theta val:", 2*np.pi/delPhi)
    theta = delTheta
    delTheta += tempTheta
""" 

"""
--------------------------------------------------------------
The following code checks for discrepencies in the covariance matrices between 
the old (real space) code and the new (harmonic space) code 
"""
"""
n = 10
zvals = np.linspace(1,2,n)
clusters = np.array([[zvals[i],0,0] for i in range(0,n)])
clusters = np.load("rand20Clusters.npy")
vMatOld = RemoteCMB(clusters, n)[0]
vMatNew = RemoteHarm(clusters,n) 

print(vMatNew - vMatOld)
"""       


