# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:33:28 2016
This python routine is for determinig a transfer function. For more details on transfer function,
see Bunn(2006)

THis is the generalized version and to be usde in COVandREL.py and its constituent classes. 
@author: arsalanadil
"""
import numpy as np, scipy.special as special, scipy.interpolate as interpolate, scipy.integrate as integrate, time
import matplotlib.pyplot as plt

qlist = np.loadtxt("qList1.csv", delimiter = ',')

zq = qlist[:,0]
qq = qlist[:,1]

nlist = np.loadtxt("nList1.csv", delimiter = ',')
zn = nlist[:,0]
nn = nlist[:,1]#array of eta values; zn is corresponding z vals

rec = -3.11293 #Conformal time at recombination

def jl(l,k):
    
    sphbessel = special.sph_jn(l,k)[0][l]
    return sphbessel
    
"""
Converts redshift(z) value to it's corresponding conformal time(eta or N) value
"""
def ZtoN(z):
    n = np.interp(z, zn, nn)
    return n

def delta(n,k,z):
    """
    Parameters:
    n : order of transfer function (affects bessel function)
    k : wave number
    z : red-shift to cluster
    """    
    zp1 = np.arange(z,11,0.01)#create a list of zprimes to sample the functino to use Simpson's rule (numerical integration)
    zp2 = np.arange(10.9,100,1)#lower redshifts require "finer" sampling compared to higher ones
    zp = np.concatenate([zp1,zp2])
    
    n1 = ZtoN(z)
    SW = (jl(n,k*(n1-rec)))#This is the first term in the transfer function
    
    nprime = np.interp(zp,zn,nn)
    temp = np.zeros(len(nprime))
    
    for x in range(len(nprime)):
        temp[x] = jl(n,k*(n1-nprime[x])) 
    
    temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
    #of eta_prime and the second multiplicative term is derivative of 
    #the "matter perturbation growth" (as a function of z) w.r.t the scale factor
    
    ISW = -6 * integrate.simps(temp1, zp)
    #print(ISW)    
    
    """
    Why am I adding the (-1j)**n factor here???! This is very confusing. 
    What were you thinking when you were originally doing this?!
    
    transFunc = (-1j)**n * (4*np.pi/3)*(SW - ISW)
    
    if(n == 2):
        return transFunc.real
    """
    transFunc = (4*np.pi/3)*(SW + ISW)
    
    return transFunc, SW, ISW

#temp = lambda k: delta(k,0)
#a = integrate.quad(temp, 0, 10)
#print(a)
"""
def powSp(k):
    if k==0:
        return 0
    else:
        return 1/k**3

plt.figure()
plt.title("Transfer Func")    
plt.plot([(delta(2,k,0) * k**2 * powSp(k)) for k in np.arange(0,10,.1)])
#plt.plot([delta(3,k,0).imag for k in np.arange(0,10,0.1)])
plt.show()
"""

"""
def powSp(k):
    if k==0:
        return 0
    else:
        return 1/k**3

plt.figure()
plt.title("Transfer Func")    
plt.plot([(np.abs(delta(2,k,0))**2 * k**2 * powSp(k)) for k in np.arange(0,10,.1)])
#plt.plot([delta(3,k,0).imag for k in np.arange(0,10,0.1)])
plt.show()
"""
"""
temp1 = lambda k: delta(k,2)[0]

#plt.figure()
#plt.title("ISW")
plt.plot([temp1(k) for k in range(0,50)])
plt.show() 
"""
"""
temp = lambda k: delta(k,1)[0]
temp1 = lambda k: delta(k,1)[1]
temp2 = lambda k: delta(k,1)[2]
plt.figure()
plt.plot([temp(k) for k in range(100)], label = 'SW')
plt.plot([temp1(k) for k in range(100)], label = 'ISW')
plt.plot([temp2(k) for k in range(100)], label = 'Transfer Function (w/o normalization)')
plt.title("SW, ISW, and Transfer Function")
plt.legend()
plt.show()
"""
"""    
d = np.zeros(100)
for j in range(100):
    d[j] = delta(j,1)
plt.figure(10)
plt.plot(d)
plt.show()
"""
"""
n = 20
table = np.zeros(n)
for t in range(1,n):
    t0 = time.time()
    [delta(k,1) for k in range(10*t,(n+1)*t)]   
    t1 = time.time()
    table[t] = t1-t0
plt.figure(2)
plt.plot(np.arange(0,10*n,10),table)
plt.show()
"""
#print("Computational time = ", (t1-t0))