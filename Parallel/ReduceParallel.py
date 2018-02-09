#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:49:57 2017

@author: arsalanadil
"""

import numpy as np

#n = 6912
#segments = 4
n = 12
segments = 3

import EqualArea as EA
nVals = EA.divisions(n,segments)#generated using EqualArea.py to divide covariance matrix into equal area segments
def parallelScatter(number):
    nValss = np.concatenate([nVals,[n]])
    return nValss[number]

varMatrix = np.zeros((n,n)).astype(complex)
relMatrix = np.zeros((n,n)).astype(complex)

#temp = np.load("results0.npz")
#tempVar = temp['arr_0']
#tempRel = temp['arr_1']
for i in range(0,segments+1):
    number = i 
    nStart = int(parallelScatter(number-1))
    nEnd = int(parallelScatter(number))
    
    if(number == 0):
        nEnd = int(parallelScatter(number))
        nStart = 0
    
    tempVar = np.load("results%s"%i+".npz")['arr_0']
    tempRel = np.load("results%s"%i+".npz")['arr_1']
    
    #print(tempVar)
    p = 0
    for j in range(nStart,nEnd):
        for k in range(0,j+1):
            varMatrix[j][k] = tempVar[p]
            varMatrix[k][j] = np.conjugate(tempVar[p])
            p+=1
    
    #tempVar = tempVar + np.load("results%s"%i+".npz")['arr_0']
    #tempRel = tempRel + np.load("results%s"%i+".npz")['arr_1']    
