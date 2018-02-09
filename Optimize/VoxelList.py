#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 17:47:29 2017

Creates an array of voxels with the respective theta, phi and z values.
Note: for Nside =8, there is a grid of 768 voxels for each z-shell.
Thus, for 13 z-shells, the total number of voxes is 13*768
(but actually there's 14 z-shells, because there's also the last one at z = 2)

@author: arsalanadil
"""
import healpy as hp
import numpy as np

Nside = 8
zVals = np.concatenate((np.load("zVox.npy"),[2]))
temp = hp.pix2ang(Nside,np.arange(12*(Nside**2)))
thetaVals = temp[0]
phiVals = temp[1]
clusters = np.zeros((len(thetaVals)*len(zVals),3))


tpList = np.zeros((len(phiVals),2))
for i in range(len(phiVals)):
    tpList[i] = [thetaVals[i],phiVals[i]]

i = 0
for z in range(0,len(zVals)):
    print("At Z =", zVals[z])
    for j in range(0,len(thetaVals)):    
        clusters[i] = np.concatenate([[zVals[z]], tpList[j]])
        i+=1
