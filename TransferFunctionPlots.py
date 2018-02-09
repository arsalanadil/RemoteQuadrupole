#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 13:04:00 2017
This routine is intended for creating plots of the tranfer function
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

def ISWTerm(k,z):
        
    l = 50*k
    if(l<50):
        l = 50
    
    zp = np.arange(l).astype(float) * (100 -z)/(l)+z#this needs to confirmed
    
    n1 = ZtoN(z)
    
    nprime = np.interp(zp,zn,nn)
    temp = np.zeros(len(nprime))
    
    for x in range(len(nprime)):
        temp[x] = j2(k*(n1-nprime[x])) 
    
    temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
    #of eta_prime and the second multiplicative term is derivative of 
    #the "matter perturbation growth" (as a function of z) w.r.t the scale factor
    
    ISW = 6 * integrate.simps(temp1, zp)
    
    return ISW

def IntegrateISW(z):
    ISWIntegrand = lambda k : ISWTerm(k,z)
    temp = quad(ISWIntegrand,0,60, limit=30, epsrel = 1e-5,epsabs = 1e-5)[0]
    
    return  temp

def SWTerm(k,z):
    n1 = ZtoN(z)
    SW = (j2(k*(n1-rec)))#This is the first term in the transfer function
    
    return SW

def IntegrateSW(z):
    SWIntegrand = lambda k : SWTerm(k,z)
    temp = quad(SWIntegrand,0,40, limit=30, epsrel = 1e-5,epsabs = 1e-5)[0]
    
    return  temp

"""
n=200
c2 = np.zeros(n)
zlist = np.logspace(0,1,n)
for i in range(0,n):
    c2[i] = (IntegrateISW(zlist[i]))**2
    print("At i value",i)
    
print("This window plots ISW-ISW vs Z")
plt.figure(3)
plt.title("ISW-ISW")
plt.semilogx(zlist, c2)
"""

"""

def delta(k,z):
        
    l = 50*k
    if(l<50):
        l = 50
    
    zp = np.arange(l).astype(float) * (100 -z)/(l)+z#this needs to confirmed
    
    n1 = ZtoN(z)
    SW = (j2(k*(n1-rec)))#This is the first term in the transfer function
    
    
    nprime = np.interp(zp,zn,nn)
    temp = np.zeros(len(nprime))
    
    for x in range(len(nprime)):
        temp[x] = j2(k*(n1-nprime[x])) 
    
    temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
    #of eta_prime and the second multiplicative term is derivative of 
    #the "matter perturbation growth" (as a function of z) w.r.t the scale factor
    
    ISW = 6 * integrate.simps(temp1, zp)
    
    transFunc = (-4*np.pi/3)*(SW - ISW)
    
    return transFunc

plt.figure(4)
plt.plot([delta(k,1) for k in np.arange(0,100,0.1)])
plt.title("Transfer Function vs K at z=1")
"""