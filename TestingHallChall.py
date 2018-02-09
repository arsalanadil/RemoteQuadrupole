#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 13:01:22 2017

This is to try and reproduce results from Hall and Challinor. Specifcally, the thing they call "zeta" in equations 16, and subsequently 
plot in Fig.2 For more details on how to equate what we have with they do (they did their math in Harmonic space and we did it in real space)
see pg. 65 in notes.
In short, zeta_{lm}(theta phi) = <p(r, theta phi) a_{lm}(r,theta,phi)>/2Ylm(theta phi)
essentially removing the constraint on the distance to the cluster. This means that the only thing we care about is the SEPERATION between two clusters. Dividing by spin-2
spherical harmonic removes that constraint as well (polarization field has the same "shape" essentially). This means in the end, for A FIXED REDSHIFT, zeta is a constant.


I'm using a routine I found online for calculating the spin-2 spherical harmonics. Note that the arguments are reversed from our notation. i.e., what we call theta they call phi
and vice versa. 

@author: arsalanadil
"""

import VarMatrix_Gen_Local as Cross
import Spin2Harmonics as SWSH
import numpy as np, matplotlib.pyplot as plt

def clusterList(z,n):
    """
    Generates a random list of clusters at fixed red-shift
    """
    clusters = np.zeros((n,3))
    phiList = np.random.uniform(0.0, 2*np.pi,n)
    thetaList = np.arccos(np.random.uniform(-1,1,n))

    for i in range(0,n):
        clusters[i] = np.array([z,thetaList[i],phiList[i]])
        
    return clusters
"""
First, we want to test whether we get directional independence. The result should also be independent of the "m" we choose,
as long as we choose the same one for both the numerator (i.e. the cross covariance) and the denominator (i.e. the spin-2 harmonic).

Yes we do. This code also agrees with the results (for l=2) with that from VarMatrixGeneralized.py.

"""
"""        
z = 0.425#red-shift (needs to be fixed for this procedure)
n = 2#numebr of clusters to compare (each should yield the final same value of zeta for constant z). Need for "clusterList"
l = 3#local multipole of interest
m = -2#corresponding m-value. Note that zeta should be independent of this term, but that this is only for testing purposes.
#mLength = 2*l+1

#clusters = clusterList(z,n)

V1 = np.array([z,0.1,0.3])
V2 = np.array([z,1.1,2.2])
#V3 = np.array([z,2.3,0.69])
clusters = np.vstack((V1,V2))
result = np.zeros(n).astype(complex)
#mlist = np.zeros(mLength).astype(complex)
for i in range(0,n):
    flm =Cross.wignerRotation(clusters[i],l,Cross.powerSpec)[0][0][l+m]
    ylm = SWSH.sYlm(2,l,m,clusters[i][2],clusters[i][1])
    result[i] = flm/ylm
    print("flm", flm,"ylm",ylm,"....","flm/ylm", result[i])
"""






"""
We can now test this against Fig. 2. Note that if this works, that's really really good news because we will then have an external check on our work!

"""
z=0.025#0.025,0.425,0.975
lmax = 2
m = -2

result = np.zeros(lmax-1).astype(complex)
toPlot = np.zeros(lmax-1)
for j in range(2,lmax+1):
    flm = Cross.wignerRotation(np.array([z,1,1]),j, Cross.powerSpec)[0][0][j+m]
    print("FLM:", flm)
    ylm = SWSH.sYlm(2,j,m,1,1)
    result[j-2] = (flm/ylm)
    toPlot[j-2] = np.abs(flm/ylm)
""" 
         
plt.plot(np.arange(2,lmax+1), toPlot, label="z = 0.025")

z=0.425#0.025,0.425,0.975
lmax = 5
m = -2

result1 = np.zeros(lmax-1).astype(complex)
toPlot1 = np.zeros(lmax-1)
for j in range(2,lmax+1):
    flm = Cross.wignerRotation(np.array([z,1,1]),j, Cross.powerSpec)[0][0][j+m]
    ylm = SWSH.sYlm(2,j,m,1,1)
    result1[j-2] = (flm/ylm)
    toPlot1[j-2] = np.abs(flm/ylm)
 
         
plt.plot(np.arange(2,lmax+1), toPlot1, label="z=0.425")

z=0.975#0.025,0.425,0.975
lmax = 5
m = -2

result2 = np.zeros(lmax-1).astype(complex)
toPlot2 = np.zeros(lmax-1)
for j in range(2,lmax+1):
    flm = Cross.wignerRotation(np.array([z,1,1]),j, Cross.powerSpec)[0][0][j+m]
    ylm = SWSH.sYlm(2,j,m,1,1)
    result2[j-2] = (flm/ylm)
    toPlot2[j-2] = np.abs(flm/ylm)
 
         
plt.plot(np.arange(2,lmax+1), toPlot2,label= "z =0.975")
"""

