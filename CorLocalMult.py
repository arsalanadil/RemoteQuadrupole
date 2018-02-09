# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 19:25:31 2016

Previously, in the code VarMatrixGneneralized.py, we created elements of the Variance (and Corelation) matrices
which symbolized the corelation between remote quadrupoles from two different galaxy clusters. VarAndCorMat.py is a routine
that used this process to generate the entire Variance and Corelation matrices. 

Now we're interested in finding the corelation between our local MULTIPOLE (multiple l values) -which is well known-
with a remote quadrupole from a distant galaxy cluster. This routine is an attempt to do so. 

Alot of the code here is similar to VarMatrix_Genereralized.py - but alot of it isn't.


@author: arsalanadil
"""

import scipy.special as special
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import TransFuncGen as tf
import ClebGord as cg
from sympy.physics.quantum.spin import Rotation
import sympy as sp
import time

def bessel(l,q):
    """
    Parameters:
    q : product of k.del(r) 
    l : order of the spherical bessel function; 
    
    ToDo: Instead of closed form (form Mathematica), think about just using special.sph_jn(l,q) from python
    """
    if l == 0:
        if q == 0:
            return 1
            
        return np.sin(q)/q
    elif l == 2:
        if q == 0:
            return 0
        return (-3*np.cos(q))/q**2 + ((3 - q**2)*np.sin(q))/q**3
    
    elif l == 4:
        if q == 0:
            return 0
        return (5*(-21 + 2*q**2)*np.cos(q))/q**4 + ((105 - 45*q**2 + q**4)*np.sin(q))/q**5
    else:
        jn = special.sph_jn(l,q)[0][2] #this routine in python computes the spehrical bessel function and its first
        #derivates. I only want one value. See documentation
        return jn
        
        
        
        
def KIntegral(l,z,L,ps):
    """
    Calculates the "K-integral". See notes on June 30th '16    
    
    Parameters:
    l order of local multipole
    z: red-shift to cluster
    L: sum over these values
    """
    
    R = -tf.ZtoN(z)    
    
    Klm = lambda k: (   (k**2 * ps(k) * tf.delta(2,k,z) * 
                         np.conjugate(tf.delta(l,k,0)) * bessel(L,(k*R)))    )
    """
    Question: What are my parameters for the bessel function in the above?
    the second parameter doesn't make sense... I think. Should I just be multiplying them???
    Where am I using the theta and phi values of the cluster. It seems like this is invariant under
    rotation - which shouldn't be the case.
    """
    
    temp = quad(Klm,0,10, limit=20, epsrel = 1e-5,epsabs = 1e-5)
    #choose an upper limit for integration carefully. This function dies down really fast
    #so k=10 is not bad.
    return temp[0]
    
def KBar(L,l,m1,m2,z,ps):
    """
    Calculates a matrix element for <a2m'(r) a*lm(0)>
    
    l,m1:    local multipole
    m2: remote quadrupole (2,m2)
    
    
    """
    M = m2 - m1 
    
    cleb = cg.J(2,l,L,m2,m1,M)
    
    print("l", l, "L", L, "m2", m2, "m1", m1, "M", M, "clebsh", cleb)
    KlmBar = (4*np.pi*KIntegral(l,z,L,ps)*cleb)*(1j)**L
    
    return KlmBar

def  KBarCor(l,m1,m2,z,ps):

    temp = 0
    
    for x in range(-2,3,2):
        L = l + x
        if(np.abs(m2-m1) <= L):  
            temp += KBar(L,l,m1,m2,z,ps)
    
    return temp
    
def CorMatrix(l,z,ps):
    """
    The multipole upto which we want our corelation
    """
    
    KMat = np.zeros((5,(2*l + 1))).astype(complex)
    
    for m2 in range(-2,3,1):
        for m1 in range(-l,l+1,1):
            KMat[m2][m1] = KBarCor(l,m1,m2,z,ps)
            
    return KMat
    
def WignerRotation(cluster,l,ps):
    
    z = cluster[0]
    t1 = cluster[1]
    p1 = cluster[2]
    
    elDm = lambda m1, m2: Rotation.D(l,m1, m2,-t1,p1,0).doit()
    MatDm = np.zeros((2*l+1, 5)).astype(complex)
    
    for m2 in range(-2,3):
        for m1 in range(-l,l+1):
            MatDm[m1+l][m2+2] = elDm(m1,m2)
            
    temp = CorMatrix(l,z,ps)
    corelation = np.linalg.multi_dot([MatDm,temp])
    
    return corelation
    
def powerSpec(k):
    if(k==0):
        return 0
    else:
        return 1/k**3
        
v1 = np.array([0,0,0])
temp = WignerRotation(v1,2,powerSpec)
print(temp)