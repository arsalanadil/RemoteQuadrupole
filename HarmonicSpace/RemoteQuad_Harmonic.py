# -*- coding: utf-8 -*-
"""
Spyder Editor

June 15th 2017: 
    Till now we've done everything in real space; we're now going to transfer over to harmonic space. 
This is the remote quadrupole code. For details, see your notes and the paper by Hall and Challinor; Seto and Pierpaoli;
Hu&White; and E.Bunn

The methods "interpTF, j2, delta, ZtoN" are used for calculating the transfer function.

"""

import numpy as np, scipy.special as special, scipy.interpolate as interpolate
import scipy.integrate as integrate, sympy as sp
from sympy.physics.quantum.spin import Rotation
import matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import scipy
from scipy.misc import factorial
import numba as nb
from numba import float64, int32, jit, complex64, cfunc
from mpl_toolkits.mplot3d import Axes3D

qlist = np.loadtxt("qList1.csv", delimiter = ',')

zq = qlist[:,0]
qq = qlist[:,1]

nlist = np.loadtxt("nList1.csv", delimiter = ',')
zn = nlist[:,0]
nn = nlist[:,1]#array of eta values; zn is corresponding z vals

rec = -3.11293 #Conformal time at recombination

Interpolate = 'ON'#OFF or ON

interpTable = np.load("ISWL2.npz")
kk = interpTable['kVals']
zz = interpTable['zvals']
isw = interpTable['isw']
ISWinterp = interpolate.RectBivariateSpline(kk,zz,isw)

@jit('float64(float64,float64)')
def interpTF(k,z):
    ISW = ISWinterp(k,z)[0][0]
    n1 = ZtoN(z)
    SW = (bessel(2, k*(n1-rec)))#This is the Sachs-Wolfe term
    
    transFunc = 1/3 * (SW + ISW)#note that the -ve sign is due to (-1j)**2 
    return transFunc 


#Spherical bessel function of second order
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
@jit('float64(float64)')
def ZtoN(z):
    n = np.interp(z, zn, nn)
    return n
@jit('float64(float64,float64)')
def delta(k,z):
    
    """    
    l = 50*k
    if(l<50):
        l = 50
    
    zp = np.arange(l).astype(float) * (100 -z)/(l)+z#this needs to confirmed
    """
    
    n1 = ZtoN(z)
    SW = (j2(k*(n1-rec)))#This is the Sachs-Wolfe term
    
    zp1 = np.arange(z,11,0.01)#create a list of zprimes to sample the functino to use Simpson's rule (numerical integration)
    zp2 = np.arange(10.9,100,1)#lower redshifts require "finer" sampling compared to higher ones
    zp = np.concatenate([zp1,zp2])
    
    nprime = np.interp(zp,zn,nn)#niterpolate "nlist" which converts Z to Eta.
    temp = np.zeros(len(nprime))
    
    for x in range(len(nprime)):
        temp[x] = j2(k*(n1-nprime[x])) 
    
    temp1 =  temp*np.interp(zp,zq,qq) #Here, temp is the bessel funciton (evaluated at different values
    #of eta_prime and the second multiplicative term is derivative of 
    #the "matter perturbation growth", i.e. the growth function D, divded by scale factor a (as a function of z)
    
    #ISW = -6*integrate.simps(temp1, zp)#minus sign is there because our qlist is positive.
    
    #transFunc = -4*np.pi/3 * (SW + ISW)#note that the -ve sign is due to (-1j)**2 
    ISW = -1*integrate.simps(temp1, zp)
    transFunc = 1/3 * SW + 2 * ISW #following Hall and Challinor convention here
                              
    return transFunc

def powerSpec(k,ns):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**(4-ns)

def bessel(l,q):
    """
    Parameters:
    q : product of k.del(r) 
    l : order of the spherical bessel function; sensible values = 0,2,4
    """
    
    """
    if l == 0:
        if q == 0:
            return 1
            
        return np.sin(q)/q
    elif l == 2:
        if q == 0:
            return 0
        return (-3*np.cos(q))/q**2 + ((3 - q**2)*np.sin(q))/q**3
    
    else:
        sphbessel = special.sph_jn(l,q)[0][l]
        return sphbessel
    """
    sphbessel = special.sph_jn(l,q)[0][l]
    return sphbessel
@jit('float64(float64, int32, float64, float64)')
def func(k,l,z1,z2):
    r1 = -ZtoN(z1)#we want "magnitude" of r. hence need positive value
    r2 = -ZtoN(z2)
    if (Interpolate == 'OFF'):
        return 10**8 * 1/k**2 * bessel(l, k*r1)/r1**2 * bessel(l,k*r2)/r2**2 * delta(k,z1) * delta(k,z2) * powerSpec(k,1)
    return 10**8 * 1/k**2 * bessel(l, k*r1)/r1**2 * bessel(l,k*r2)/r2**2 * interpTF(k,z1) * interpTF(k,z2) * powerSpec(k,1)

@jit('float64(float64, float64, int32)')
def KIntegral(z1,z2,l):#this returns the function "\xi_l(z,z')
    
    #const = 9/(8 * np.pi) * factorial(l+2)/factorial(l-2)
    const= 9 * np.pi / 4 * factorial(l+2)/factorial(l-2)
    
    #func = lambda k: 1/k**2 * bessel(l, k*r1)/r1**2 * bessel(l,k*r2)/r2**2 * delta(k,z1)[0] * delta(k,z2)[0] * powerSpec(k,1)
    
    #return func
    
    integral = integrate.quad(func,0,300, limit = 50, args=(l,z1,z2))[0]/10**8
    
    return const * integral

def UnrotCovRel(beta,z1,z2):
    lmax = 8
    
    xi = np.zeros(lmax-1).astype(complex)
    for l in range(2,lmax+1):
        xi[l-2] = KIntegral(z1,z2,l)
        #print(xi[l-2])
    
    rel = 0
    for l in range(2,lmax+1):
        rel += (2*l + 1)/(4*np.pi) * xi[l-2] * sp.N(Rotation.d(l,2,-2,beta).doit())
        
    var = 0
    for l in range(2,lmax+1):
        var += (2*l + 1)/(4*np.pi) * xi[l-2] * sp.N(Rotation.d(l,2,2,beta).doit())
        
    return var,rel

from numpy import sin, cos, arccos, arctan2, arcsin
def Rotate(t1,p1,t2,p2):
    l = np.abs(p1-p2)
    #l = p2-p1
    beta = arccos(cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2))
    if(beta == 0):
        return np.array([0,0,0])
    alpha = arcsin(sin(l) * sin(t1)/sin(beta))
    #gamma = arctan2(np.tan(alpha), cos(beta))
    meem = arcsin(sin(l)*sin(t2)/sin(beta))#note that gamma = np.pi-meem (see notes)
    
    #return np.array([alpha,beta,gamma, np.pi-meem, meem])
    #return np.array([alpha,beta,meem])
    return np.array([alpha,beta,np.pi-meem])

def CovRel(v1,v2):
    euler = Rotate(v1[1],v1[2],v2[1],v2[2])
    alpha = euler[0]
    beta = euler[1]
    gamma = euler[2]
    print("Angles:", euler)
    
    unrot = UnrotCovRel(beta,v1[0],v2[0])
    #print("unrot", unrot)
    VarRot = sp.N(np.e**(2 * 1j* (alpha - gamma))*unrot[0])
    #VarRot = sp.N(np.e**(2 * 1j* (alpha+gamma))*unrot[0])
    RelRot = sp.N(np.e**(2 * 1j* (-alpha - gamma))*unrot[1])
    
    return VarRot,RelRot

def XiInterpolate(l,spacing):
    #zvals = np.arange(0,2.1,spacing)
    #zvals = np.concatenate([[0],2*np.logspace(-3,0,30)])

    z1 = np.arange(0,0.12,0.02)
    z2 = np.arange(0.12,0.50,0.08)
    z3 = np.arange(0.50,2.1,0.1)
    zvals = np.concatenate([z1,z2,z3])
    temp = [KIntegral(z,z,l) for z in zvals]
    interp = interpolate.interp1d(zvals, np.transpose(temp))
    return interp

def XiInterpolate2(z1,l):
    #zvals = np.arange(0,2.1,spacing)
    #zvals = np.concatenate([[0],2*np.logspace(-3,0,30)])

    zz1 = np.arange(0,0.12,0.02)
    zz2 = np.arange(0.12,0.50,0.08)
    zz3 = np.arange(0.50,2.1,0.1)
    zvals = np.concatenate([zz1,zz2,zz3])
    temp = [KIntegral(z1,z,l) for z in zvals]
    interp = interpolate.interp1d(zvals, np.transpose(temp))
    return interp
"""    
philist = np.arange(0,2*np.pi+0.2,0.2)
temp = np.zeros(len(philist))
j = 0
for i in philist:
    cov = CovRel([0.5,0,0],[0.5,1,i])[0]
    #temp[j] = np.abs(1j*(cov - np.conjugate(cov))/2)#to plot only the imaginary part
    temp[j] = (cov + np.conjugate(cov))/2            #plot only real part
    j+=1
plt.plot(philist,temp)        
"""
#temp = KIntegral(2,2,4)
#plt.plot([temp(k) for k in np.arange(0.01,30,0.1)])

#plt.semilogy(np.arange(0.1,2.1,0.1),[10**4*KIntegral(i,i,4) for i in np.arange(0.1,2.1,0.1)])

#l = 4
#H20 = 5*KIntegral(0,0,2)
#plt.plot(np.arange(0,2.1,0.05),[(2*l + 1)*KIntegral(i,i,l)for i in np.arange(0,2.1,0.05)]/H20)

#plt.semilogy(np.arange(2,100),[10**4 * KIntegral(2,2,l) for l in range(2,100)])
"""    
l = 8
zlist = np.arange(0,2,0.01)
interp = XiInterpolate(l,0.08)
interPlot = [interp(z) for z in zlist]
realPlot = [KIntegral(z,z,l) for z in zlist]
plt.figure()
plt.plot(zlist, interPlot)
plt.plot(zlist, realPlot)
"""
"""
C20 = KIntegral(0,0,2)
l=7
z1 = 0.5
zlist = np.arange(0,2,0.01)
plt.figure()
interp = XiInterpolate2(z1,l)
plt.plot(zlist,[KIntegral(z,z1,l) for z in zlist]/C20)
plt.plot(zlist, [interp(z) for z in zlist]/C20)
"""
"""
#Let's make some 3D plots now;
l = 2        
zlist = np.arange(0,2.2,0.2)
T = np.zeros((len(zlist),len(zlist)))
plotList = np.zeros((len(zlist),len(zlist)))
for i in range(0,len(zlist)):
    for j in range(0,len(zlist)):
        plotList[i][j] = KIntegral(zlist[i],zlist[j],l)
        T[i][j] = np.abs(zlist[j] - zlist[i])
        print(zlist[i],zlist[j],plotList[i][j])

X,Y = np.meshgrid(zlist,zlist)
plt.figure()
cp = plt.contourf(X,-T,plotList)            
plt.colorbar(cp)
#ax = Axes3D(fig)
#ax.plot_surface(X,Y,plotList)
"""