# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 13:39:17 2016

This is a modified version of the python file, VarianceMatrix. This routine adds the integrated Sachs-Wolfe term
(i.e. the entire transfer function) in calculating the variance matrix. Moreover, it takes as input redshift values
for galaxy clusters rather than conformal times.
As before, it is still restricted to the case of galaxy clusters on the "z-axis" or "poles".

@author: arsalanadil
"""
import scipy.special as special
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import TransferFunction as tf

def F(kq):
    
    if(abs(kq)<0.001):
        return 1            
    
    f =  (-15*(3*kq*np.cos(kq) + (-3 + (kq)**2)*np.sin(kq)))/(kq)**5 
    return f

def Gamma2(zj, zk, powerSpec):
    
    Nj = tf.ZtoN(zj)#These are corresponding conformal time values
    Nk = tf.ZtoN(zk) 
    
    #print(Nj)
    #print(Nk)
    q = (Nj-Nk)
    
    func = lambda k: ((k**2)*powerSpec(k)*F(k*q)*np.real((tf.delta(2,k,zj))*np.real(tf.delta(2,k,zk))))
        
  
    temp =  quad(func, 0, 20, limit=30)
    #technically, we want the integral from k=0 to infinity which is computationally very intensive, but plotting
    #the integrand(func) reveals that it falls to zero fairly quickly.    
    
    return temp[0]
    #temp[1] corresponds to the uncertainty in the integral


zj = 1
zk = 3
#Zrec = 1100
#R = 1
#p = lambda k: 1/k
def p(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**3
        
print(Gamma2(zj,zk,p))

"""
n=10
Variance = np.zeros((n,n))
for x in range(n):
    for y in range(n):
        Variance[x][y] = np.real(Gamma2(x,y,p))

#print(Variance)

A = np.linalg.cholesky(Variance)
print(A) 
"""
#print(Gamma2(zj,zk,p))