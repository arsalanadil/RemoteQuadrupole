# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:08:55 2017

Making this routine to find the angular power spectrum, the Cl's that is. 

Cl = <alm(0) alm*(0)> = Integrate[|delta|^2 Pk j_l^2(k,R) k^2 dk]

This is only valid for diagonal elements. That's because all other elements are zero.
i.e. the corelation between different multipoles is the zero. 

@author: Arsalan
"""

import scipy.special as special
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import TransFuncGen as tf
import sympy as sp
import time

#R: distance to SLS (in conformal time)


def bessel(l,q):
    """
    Parameters:
    q : product of k.del(r) 
    l : order of the spherical bessel function; L=l, l-1, l+2
    """    
    sphbessel = special.sph_jn(l,q)[0][l]
    return sphbessel


def Integral(l, ps):
    #R = -3.11293
    #R = 2.867681345810709 
    #func = lambda k: ( k**2 * ps(k)  * np.abs(tf.delta(l,k,0))**2 * ( bessel(l,(k*R))**2 ) )#remove Bessel function.
    #the final normalization factor is: 1.89767320555/0.00740471 = 256.2792068224144 (where numerator is generated using VarMatrix_Generalized.py)
    func = lambda k: ( k**2 * ps(k)  * np.abs(tf.delta(l,k,0)[0])**2   )
    #should bessel have k*R??

    #plt.plot([func(k) for k in range(0,40)])
    #plt.show()
        
    
    #temp =  1/(9*(np.pi)**3) *  quad(func,0,300, limit=50, epsrel = 1e-5,epsabs = 1e-5)[0]
    temp =  quad(func,0,300, limit=50, epsrel = 1e-5,epsabs = 1e-5)[0]
    #for high l vals, use limit = 200.
    #print(temp)
    return temp


def ClMatrix(l,ps):
    
    clmat = np.zeros((l-1,l-1))
    
    for i in range(2,l+1):
        clmat[i-2][i-2] = Integral(i,ps)
    
    return clmat


def powerSpec(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**3

def Dl(l):
    cll = Integral(l,powerSpec)
    return 100*l*(l+1)*cll/(2*np.pi)

"""
temp = ClMatrix(2, powerSpec)
print(temp)
"""
"""
lmax = 60
lList = np.arange(2,lmax+1)
temp = np.zeros(lmax-1)
for l in range(2,lmax+1):
    temp[l-2] =  l*(l+1)*(Integral(l, powerSpec))

print(temp)
plt.figure(4)
plt.plot(lList,temp)
"""
#plt.ylabel(r'l(l+1)C_l)