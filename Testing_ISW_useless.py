# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:39:19 2016

These are just test cases I'm running with the file "VarianceMatrix_w_ISW" in order to find
an efficient way of doing it. As it stands, the code takes VERY long to generate even one
element of the variance matrix(and we want it to generate an entrie matrix of such elements).

To Try:
1. One solution that would most definitely work is to make a 2D array consisting of different
values of k and z in the delta() method. Maybe I could even use a dictionary. Then, I could
just interpolate it and see how fast that is.

2. We've noticed that delta(k,z) slows down significantly at higher value of k. It happens to be
the case that higher values of k contribute negligibly to the overall system. As of now I'm integrating
from k=0 to infinity. One obvious choice is to plot the integrand as a funciton of k and see where it approaches
zero. Then I could just integrate upto that value

3. We could break up our integral into two parts: for low values of k we sample at a high number of points
but for higher values we can change "dx" (python lets us do that)

@author: arsalanadil
"""

import scipy.special as special
import numpy as np, matplotlib.pyplot as plt, time
from scipy.integrate import dblquad, quad
import scipy.integrate as integrate
#import TransferFunction as tf
i = 0

qlist = np.loadtxt("qlist.csv", delimiter = ',')

zq = qlist[:,0]
qq = qlist[:,1]

nlist = np.loadtxt("nlist.csv", delimiter = ',')
zn = nlist[:,0]
nn = nlist[:,1]#array of eta values; zn is corresponding z vals

rec = -3.11293 #Conformal time at recombination

#Spherical bessel function of second order
def j2(k):
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
    global i    
    
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
    #print(transFunc)    
    i += 1
    return transFunc


def F(kq):
    
    if(kq<0.1):
        return 1            
    
    f =  (-15*(3*kq*np.cos(kq) + (-3 + (kq)**2)*np.sin(kq)))/(kq)**5 
    return f

def Gamma2(zj, zk, powerSpec):
    
    Nj = ZtoN(zj)#These are corresponding conformal time values
    Nk = ZtoN(zk) 
    
    #print(Nj)
    #print(Nk)
    q = (Nj-Nk)
    
    func = lambda k: ((k**2)*powerSpec(k)*F(k*q)*delta(k,zj)*delta(k,zk))
    
    """
    def func(k):
        fff =  ((k**2)*powerSpec(k)*F(k*q)*tf.delta(k,zj)*tf.delta(k,zk))
        #print(k,fff)
        return fff
    """
    
    #temp =  quad(func, 0, np.inf, limit = 2000, epsrel = 1e-4, epsabs = 0)
    #experiment with: limits of integration, limit of subdivisions, transferFunction    
    #temp =  quad(func, 200,np.inf, limit = 20)
    #k = 0 to 200, limit = default(50): call to delta() = 4158; time=100; temp = 1.153 +- 0.0336;
    #k = 0 to 200, limit = 20: call to delta() = 1638; time=35.8; temp =1.151 +- 0.0690;
    #k = 0 to 200, limit = 10: call to delta() = 798; time=23.33; temp =1.155 +- 2.935; 
    # k = 0 to inf, limit = default(50): call to delta() = 2970; time = 108; temp = 1.141 +- 0.0149;
    #k = 200 to inf, limit = 20: call to delta() = 1170; time = 132.9; temp =0.00574+-0.00531  
    
    #temp =  quad(func, 200, np.inf, limit = 50)[0]
    
    temp = integrate.trapz([func(k) for k in range(200,2000,10)])  
    return temp

#Zrec = 1100
#R = 1

zj = 3
zk = 1
p = lambda k: 1/k**2
t0 = time.time()
print(Gamma2(zj,zk,p))
t1 = time.time()
print("Time it took = ", (t1-t0))
print("Number of time delta() is called: ", i)
"""
n=2
Variance = np.zeros((n,n))
for x in range(n):
    for y in range(n):
        Variance[x][y] = np.real(Gamma2(x,y,p))

print(Variance)

A = np.linalg.cholesky(Variance)
print(A) 
"""
