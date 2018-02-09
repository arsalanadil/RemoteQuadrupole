#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:59:30 2017

This class is meant to write, and test, the idea of creating interpolation tables for tabulating the transfer funciton.

There are a few things we can do:
    interpolate only the ISW term 
    interpolate the entire transfer function
    interpolate the KIntegral (highly unlikely -- too many variables)

In anycase, this may or may not turn out to be a very good idea considering the bessel functino is very
oscillatory. 

@author: arsalanadil
"""

import numpy as np, scipy.special as special, scipy.interpolate as interpolate, scipy.integrate as integrate, time
import matplotlib.pyplot as plt
import numba as nb
from numba import float64, int32, jit
from mpl_toolkits.mplot3d import Axes3D
import time


qlist = np.loadtxt("qList1.csv", delimiter = ',')

zq = qlist[:,0]
qq = qlist[:,1]

nlist = np.loadtxt("nList1.csv", delimiter = ',')
zn = nlist[:,0]
nn = nlist[:,1]#array of eta values; zn is corresponding z vals

rec = -3.11293 #Conformal time at recombination

@jit(float64(int32,float64))
def jl(l,k):
    q = k
    if l == 0:
        if q == 0:
            return 1
            
        return np.sin(q)/q
    
    elif l == 2:
        if q == 0:
            return 0
        return (-3*np.cos(q))/q**2 + ((3 - q**2)*np.sin(q))/q**3
    
    elif l ==3:
        if q ==0:
            return 0
        return (   ((-15+q**2) * np.cos(q))/q**3 - (3* (-5+2*q**2) * np.sin(q) )/q**4   )
    
    elif l == 4:
        if q == 0:
            return 0
        return (5*(-21 + 2*q**2)*np.cos(q))/q**4 + ((105 - 45*q**2 + q**4)*np.sin(q))/q**5
    
    else:
        #print("L must be 0, 2 or 4")
        sphbessel = special.sph_jn(l,k)[0][l]
        return sphbessel
    #sphbessel = special.sph_jn(l,k)[0][l]
    #return sphbessel

@jit(float64(float64))
def j2(k):
    #bessel, der = special.sph_jn(2, k)#There's nothing wrong with this, 
    #return bessel[2]
    
    if(k==0):
        return 0
    else:
        j2 = -((3*k*np.cos(k) + (-3 + k**2)*np.sin(k))/k**3)
        return j2

"""
Converts redshift(z) value to it's corresponding conformal time(eta or N) value
"""
@jit(float64(float64))
def ZtoN(z):
    n = np.interp(z, zn, nn)
    return n

@jit(float64(float64,float64))
def delta2(k,z):
    n1 = ZtoN(z)
    
    zp1 = np.arange(z,11,0.01)#create a list of zprimes to sample the functino to use Simpson's rule (numerical integration)
    zp2 = np.arange(10.9,100,1)#lower redshifts require "finer" sampling compared to higher ones
    zp = np.concatenate([zp1,zp2])
    #note that we are actually supposed to integrate upto z=1100 but the function drops rather rapidly so that we can 
    #ignore those. This needs more testing. 
    
    nprime = np.interp(zp,zn,nn)#niterpolate "nlist" which converts Z to Eta.
    temp = np.zeros(len(nprime))
    
    for x in range(len(nprime)):
        temp[x] = j2(k*(n1-nprime[x])) 
    
    temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
    #of eta_prime and the second multiplicative term is derivative of 
    #the "matter perturbation growth", i.e. the growth function D, divded by scale factor a (as a function of z)
    
    ISW = -6*integrate.simps(temp1, zp)#minus sign is there because our qlist is positive.
    
    #transFunc = -4*np.pi/3 * (SW + ISW)#note that the -ve sign is due to (-1j)**2 
    
    #return transFunc, SW, ISW
    return ISW

def delta2ISW(ks,zs):
    
    
    #ISWtable = np.zeros((kk.shape[1],zz.shape[0]))
    ISWtable = np.zeros((len(ks),len(zs)))
    
    i = 0
    for k in ks:
        print("i =,", i)

        j = 0
        for z in zs:
            n1 = ZtoN(z)
            zp1 = np.arange(z,11,0.01)#create a list of zprimes to sample the functino to use Simpson's rule (numerical integration)
            zp2 = np.arange(10.9,100,1)#lower redshifts require "finer" sampling compared to higher ones
            zp = np.concatenate([zp1,zp2])
            #note that we are actually supposed to integrate upto z=1100 but the function drops rather rapidly so that we can 
            #ignore those. This needs more testing. 
            
            nprime = np.interp(zp,zn,nn)#niterpolate "nlist" which converts Z to Eta.
            temp = np.zeros(len(nprime))
            
            for x in range(len(nprime)):
                temp[x] = j2(k*(n1-nprime[x])) 
            
            temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
            #of eta_prime and the second multiplicative term is derivative of 
            #the "matter perturbation growth", i.e. the growth function D, divded by scale factor a (as a function of z)
            
            ISW = -6*integrate.simps(temp1, zp)#minus sign is there because our qlist is positive.
            ISWtable[i][j] = ISW
            
            j+=1
        i+=1        
    #return np.transpose(ISWtable)
    return ISWtable

def deltaL(n,k,z):
    """
    Parameters:
    n : order of transfer function (affects bessel function)
    k : wave number
    z : red-shift to cluster
    """
    if(n<2):
        print("l<2 in Transfer Function not allowed")
        
    if(n==2):
        #print("In if statement")
        return -delta2(k,z)#-ve sign: since for general l, we have taken care of the factor (-i)**l within the
        #cross correlation code. Note including this would be double counting.
    
    zp1 = np.arange(z,11,0.01)#create a list of zprimes to sample the functino to use Simpson's rule (numerical integration)
    zp2 = np.arange(10.9,100,1)#lower redshifts require "finer" sampling compared to higher ones
    zp = np.concatenate([zp1,zp2])
    
    n1 = ZtoN(z)
    
    
    nprime = np.interp(zp,zn,nn)
    temp = np.zeros(len(nprime))
    
    for x in range(len(nprime)):
        temp[x] = jl(n,k*(n1-nprime[x])) 
    
    temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
    #of eta_prime and the second multiplicative term is derivative of 
    #the "matter perturbation growth" (as a function of z) w.r.t the scale factor
    
    ISW = -6 * integrate.simps(temp1, zp)
    #print(ISW)    
    
    #return transFunc, SW, ISW
    return ISW

def deltaLISW(n,ks,zs):#note that ks and zs are arrays where we want to sample the function
    ISWtable = np.zeros((len(ks),len(zs)))
    
    i = 0
    for k in ks:
        print("i =,", i)

        j = 0
        for z in zs:
            zp1 = np.arange(z,11,0.01)#create a list of zprimes to sample the functino to use Simpson's rule (numerical integration)
            zp2 = np.arange(10.9,100,1)#lower redshifts require "finer" sampling compared to higher ones
            zp = np.concatenate([zp1,zp2])
            
            n1 = ZtoN(z)
            
            
            nprime = np.interp(zp,zn,nn)
            temp = np.zeros(len(nprime))
            
            for x in range(len(nprime)):
                temp[x] = jl(n,k*(n1-nprime[x])) 
            
            temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
            #of eta_prime and the second multiplicative term is derivative of 
            #the "matter perturbation growth" (as a function of z) w.r.t the scale factor
            
            ISW = -6 * integrate.simps(temp1, zp)
            ISWtable[i][j] = ISW
            
            j+=1
        i+=1
    return ISWtable        


def deltaLZero(n,ks,z):#for calculating local multipole at fixed redshift
    iswtable = np.zeros(len(ks))
    i = 0
    for k in ks:
        print("i =,", i)
        zp1 = np.arange(z,11,0.01)#create a list of zprimes to sample the functino to use Simpson's rule (numerical integration)
        zp2 = np.arange(10.9,100,1)#lower redshifts require "finer" sampling compared to higher ones
        zp = np.concatenate([zp1,zp2])
            
        n1 = ZtoN(z)
            
            
        nprime = np.interp(zp,zn,nn)
        temp = np.zeros(len(nprime))
            
        for x in range(len(nprime)):
            temp[x] = jl(n,k*(n1-nprime[x])) 
            
        temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
            #of eta_prime and the second multiplicative term is derivative of 
            #the "matter perturbation growth" (as a function of z) w.r.t the scale factor
            
        ISW = -6 * integrate.simps(temp1, zp)
        iswtable[i] = ISW
        i+=1
    
    return iswtable

"""
#tests for delta2
kk1 =  np.arange(0,50,0.1)
kk2 = np.arange(50,300,2)
kkk = np.concatenate([kk1,kk2])
zz = np.arange(0,4,0.01)
#kk, zz = np.meshgrid(k,z)
zVals = delta2ISW(kkk,zz)
t0 = time.time()
ISWinterp = interpolate.RectBivariateSpline(kkk,zz,zVals)
t1 = time.time()
print("Time taken in interpolation:", t1-t0)
print("Concatenated array wala hai ye")

k0 =50 
k1 = 60
z1 = 0.33
plt.plot([k for k in range(k0,k1)], [delta2(k,z1) for k in range(k0,k1)], 'ro-', 
          [k for k in range(k0,k1)], [ISWinterp(k,z1)[0] for k in range(k0,k1)] )

#np.savez_compressed("ISWL2",kVals = kkk,zvals = zz,isw = zVals)

#t0 = time.time()
#ISWinterp = interpolate.interp2d(kk,zz,zVals)
#t1 = time.time()
#print("Time taken in interpolation:", t1-t0)
"""
"""
#tests for general l
l=5
k1 =  np.arange(0,50,0.1)#need finer sampling in the region k<50
k2 = np.arange(50,300,2)
k = np.concatenate([k1,k2])
z = np.arange(0,4,0.01)
zVals = deltaLISW(l,k,z)
t0 = time.time()
ISWinterp = interpolate.RectBivariateSpline(k,z,(zVals))
t1 = time.time()
print("Time taken in interpolation:", t1-t0)


k0 = 0 
k1 = 300
z1 = 0.33
plt.plot([k for k in range(k0,k1)], [deltaL(l,k,z1) for k in range(k0,k1)], 'ro-', 
          [k for k in range(k0,k1)], [ISWinterp(k,z1)[0][0] for k in range(k0,k1)] )

#np.savez_compressed("ISWL5",kVals = k,zvals = z,isw = zVals)
"""

#but we'll only ever need l>2 at z=0. So we can use normal 1D interpolation for that purpose.
l=2
k1 =  np.arange(0,50,0.1)#need finer sampling in the region k<50
k2 = np.arange(50,300,4)
k = np.concatenate([k1,k2])
#k = np.arange(0,10,0.1)
zVals = deltaLZero(l,k,0)
t0 = time.time()
#ISWinterp = interpolate.UnivariateSpline(k, np.transpose(zVals))
ISWinterp = interpolate.interp1d(k, np.transpose(zVals))
t1 = time.time()
print("Time taken in interpolation:", t1-t0)

np.savez_compressed("ISWL2Z0",kVals = k,isw = zVals)

"""
k0 = 290 
k1 = 295
z1 = 0
plt.plot([k for k in range(k0,k1)], [deltaL(l,k,z1) for k in range(k0,k1)], 'ro-', 
          [k for k in range(k0,k1)], [ISWinterp(k) for k in range(k0,k1)] )
"""

"""
deltaMatrix = np.zeros((len(k), len(z)))
j = 0
for k1 in range(0,len(k)):
    i=0
    for z1 in range(0, len(z)):
        deltaMatrix[j][i] = delta2(k1,z1)
        i+=1
    j+= 1

fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(kk,zz,np.transpose(deltaMatrix))

interpMatrix = np.zeros((len(k), len(z)))
j = 0
for k1 in range(0,len(k)):
    i=0
    for z1 in range(0, len(z)):
        interpMatrix[j][i] = ISWinterp(k1,z1)
        i+=1
    j+= 1

ax.plot_surface(kk,zz, np.transpose(interpMatrix))
#plt.plot([k for k in range(0,10)], [ISWinterp(k,1)[0] for k in range(0,10)])
"""
