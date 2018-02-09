#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 16:17:01 2017

@author: arsalanadil
"""


import numpy as np, scipy.special as special, scipy.interpolate as interpolate, scipy.integrate as integrate, time
import matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad


qlist = np.loadtxt("qList1.csv", delimiter = ',')

zq = qlist[:,0]
qq = qlist[:,1]

nlist = np.loadtxt("nList1.csv", delimiter = ',')
zn = nlist[:,0]
nn = nlist[:,1]#array of eta values; zn is corresponding z vals

rec = -3.11293 #Conformal time at recombination

#Spherical bessel function of second order
def j2(k):
    #bessel, der = special.sph_jn(2, k)#There's nothing wrong with this, 
    #except Python has a stack overflow when integrating this to infinity
    #return bessel[2]
    
    if(k==0):
        return 0
    else:
        j2 = -((3*k*np.cos(k) + (-3 + k**2)*np.sin(k))/k**3)
        return j2
"""
Converts redshift(z) value to it's corresponding conformal time(eta or N) value
"""
def ZtoN(z):
    n = np.interp(z, zn, nn)
    return n

def delta(k,z):
    
    """    
    l = 50*k
    if(l<50):
        l = 50
    
    zp = np.arange(l).astype(float) * (100 -z)/(l)+z#this needs to confirmed
    """
    
    n1 = ZtoN(z)
    SW = 3/10*(j2(k*(n1-rec)))#This is the first term in the transfer function
    
    zp1 = np.arange(z,11,0.01)
    zp2 = np.arange(10,100,1)
    zp = np.concatenate([zp1,zp2])
    nprime = np.interp(zp,zn,nn)
    temp = np.zeros(len(nprime))
    
    for x in range(len(nprime)):
        temp[x] = j2(k*(n1-nprime[x])) 
    
    temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
    #of eta_prime and the second multiplicative term is derivative of 
    #the "matter perturbation growth" (as a function of z) w.r.t the scale factor
    
    ISW = -9/5*integrate.simps(temp1, zp)#minus sign is there because our qlist is positive.
    
    transFunc = SW + ISW
    
    return transFunc, SW, ISW

def C2(z1):
    """
    Parameters:
    k : Wavenumber
    Vr1,Vr2: position vectors of galaxy clusters
    """    
    #print("l is", l, "delR is,", delR)
    
    func = lambda k: (k**2 * powerSpec(k) * (delta(k,z1)[0]**2) )
    #func = lambda k: (k**2 * powerSpec(k) * (tf.SWTerm(k,z1)**2)  * bessel(l,(k*deltR)))
    
    #plt.plot([func(k) for k in range(0,30)])
    #plt.show()
    
    
    temp = quad(func,0,300, limit=30, epsrel = 1e-6,epsabs = 1e-6)
    #for high z vals, use upper limit of 60.....
    print("C2 at z=,",z1,"is:", temp[0])   
    print("Error in integral is,", temp[1])
    return temp[0]

def powerSpec(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**(3)

"""    
n=40
zlist = np.logspace(-3,2,n)
clusterList = np.zeros((n,3))
c2 = np.zeros(n)
for i in range(0,n):
    c2[i] = C2(zlist[i])

plt.figure(1)
print("This window plots C2(z) vs Z (total n=1)")
plt.semilogx(zlist, c2/0.0074979107580118816)
plt.title('C2(z) vs Z total(n=1)')
#plt.plot(zlist,c2)
"""
