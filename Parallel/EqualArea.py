#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 21:00:09 2017
Iz nayce; calculates equal area segments of a right triangle. 
@author: arsalanadil
"""

import numpy as np

def divisions(N,segments):
    #N = The total number of elements
    #segments = #number of segments we want to divide the triangle into
    
    nVals = np.zeros(segments)
    
    TArea =N * N#Total area of triangle (ignore factor of 2-- that get's cancelled)
    #dA = segments/N * 2 * TArea
    dA = TArea/segments#area of each segment
    
    n0 = np.sqrt(dA)
    #print("n0 =",int(n0))
    
    temp =int(n0)
    nVals[0] = temp
    for i in range(1,segments):
        n = np.sqrt(dA + temp**2)
        n = int(n)
        nVals[i] = n
        Area = 1/2 * (n + temp)* (n - temp)
        #print("n%s"%i,"=",n, "Area=", Area)
        temp = n
    return nVals