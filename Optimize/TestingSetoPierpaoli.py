#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:06:57 2017

Testing Figure 2 from Seto and Pierpaoli using our original code that was written in real space. 
To use pre-optimized code (VarMatrix_Generalized.py), blank out "import Remote Quad as remote" and instead use 
"import VarMatrix_GenOld as remote"

@author: arsalanadil
"""

import RemoteQuad as remote
#import VarMatrix_GenOld as remote
import CrossCor as cross
import numpy as np, matplotlib.pyplot as plt
import Spin2Harmonics as s2
from sympy.abc import x
import sympy as sp
from sympy.physics.quantum.spin import Rotation
import scipy.integrate as integrate
"""
def quadrupole(z):
    temp = remote.wignerRotation(np.array([z,0,0]),np.array([z,0,0]))[0]
    #return np.abs(np.sum(np.diag(temp)))
    return temp[0][0]
"""
"""
def multipole(z,l):

    integral = 0
    
    if(z == 0):
        blist = np.linspace(0,np.pi,200)
        temp = np.zeros(200).astype(complex)
        val = remote.wignerRotation(np.array([z,0,0]),np.array([z,0,0]))[0][0][0]
        for i in blist:
            #temp[i] = val * np.abs(sp.N(Rotation.d(l,-2,-2,i).doit())) * np.sin(i)
            temp[i] = val * np.cos(i/2)**4 * (-2 + 3*np.cos(i))*np.sin(i)
        integral = integrate.simps(temp,blist)
        #print(integral)
    
    else:
        #z = 2 
        #l = 2
        blist = np.linspace(0,np.pi,100)
        temp = np.zeros(100).astype(complex)
        for i in blist:
            #temp[i] = remote.wignerRotation(np.array([z,1,0]),np.array([z,1,i]))[0][0][0] * sp.N(Rotation.d(l,2,2,i).doit()) * np.sin(i)
            temp[i] = remote.wignerRotation(np.array([z,0,0]),np.array([z,0,i]))[0][0][0] * np.cos(i/2)**4 * (-2 + 3*np.cos(i))*np.sin(i)
        integral = integrate.simps(temp,blist)
    
    print(integral)
    return integral
"""

def multipole(z,l):

    integral = 0
    
    if (z == 0 and l ==2):
        val = np.real(remote.wignerRotation(np.array([z,0,0]),np.array([z,0,0]))[0][0][0])
        
        func = lambda b: val * np.cos(b/2)**4 *np.sin(b)
        
        integral = integrate.quad(func,0,np.pi)[0]
        print(integral)
        return integral
    
    elif(z == 0 and l==3):
        
        val = np.real(remote.wignerRotation(np.array([z,0,0]),np.array([z,0,0]))[0][0][0])
        
        func = lambda b: val * (np.cos(b/2))**4 * (-2 + 3*np.cos(b))*np.sin(b)
        #func = lambda b: val * sp.N(Rotation.d(l,2,2,b).doit())* np.sin(b)
        
        integral = integrate.quad(func,0,np.pi)[0]
        print(integral)
        return integral
    
    else:
        func = lambda b: np.real(remote.wignerRotation(np.array([z,0,0]),np.array([z,b,0]))[0][0][0]) * (np.cos(b/2))**4 * (-2 + 3*np.cos(b))*np.sin(b)
        #func = lambda b: np.real(remote.wignerRotation(np.array([z,0,0]),np.array([z,b,0]))[0][0][0]) * sp.N(Rotation.d(l,2,2,b).doit()) * np.sin(b)
        integral = integrate.quad(func,0,np.pi)[0]
        print(z,integral)
        return integral
    


"""
def multipole2(z,l):
         
    func = lambda i: np.real(remote.wignerRotation(np.array([z,0,0]),np.array([z,0,i]))[1][0][0]) * Rotation.d(l,2,2,i).doit() * sp.sin(i)
    
    temp = sp.mpmath.quad(func,[0,np.pi])
    
    return temp[0]
"""

"""
n = 10
l = 3
zmax = 2
#zmax = 1.5
C2 = np.zeros(n)
#C20 = quadrupole(0)
C20 = multipole(0,2)
 
i=0
zlist = np.linspace(0.2,zmax,n)

for z in zlist:
    C2[i] =  np.real(multipole(z,l) ) 
    #C2[i] = np.abs(multipole(z,l) )/C20 
    i+=1

#plt.figure(2)
#plt.plot(np.concatenate([[0],zlist]), np.concatenate([[1],C2]))
plt.plot(zlist,C2/C20)
"""
"""
n = 5
zmax = 2
C21= np.zeros(n)
i=0
zlist1 = np.arange(1,zmax,zmax/n)

for z in np.arange(1,zmax,zmax/n):
    #C2[i] = 7/5 * np.real(multipole(z,l) )/C20 
    C21[i] = np.real(multipole(z,l) )/C20 
    i+=1

plt.plot(zlist1, C21)
"""
"""
thetalist = np.linspace(0,np.pi,20)
plt.plot(thetalist, [remote.wignerRotation(np.array([1,0,0]), np.array([1,0,i]))[0][0][0] for i in thetalist])
"""
""" 
n = 11
l = 2
xhilist = np.zeros(n)
C20 = 5 / (3*np.pi) * quadrupole(0)
i=0
for z in np.arange(0,2.2,0.2):
    xhilist[i] = (2*l+1)/(3*np.pi) * quadrupole(z)/C20 
    #xhilist[i] = multipole(z,l)
    i+=1

plt.plot(np.arange(0,2.2,0.2), xhilist)
"""
