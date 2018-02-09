# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:46:01 2016

This is one of the more useful (in the sense that it is an active part of the simulations) pieces of code.

Previously, (in VarianceMatrix_w_ISW.py) we had added an additional layer of depth by incorporating
the entire Transfer funciton(w/ the ISW term). However, that piece of code was still entirely built 
on the assumption that we consider galaxy clusters only on the z-axis. This is an attempt
to generalise that argument to the entire sphere.

This routine generates elements of the covariance matrix <a2m1 a2m2*> and builds upon code written previously;
most notably on TransferFunction and ClebGord(routine that calculates the integral over three Ylm coefficients).

@author: arsalanadil
"""

import scipy.special as special, scipy.interpolate as interpolate
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
#import TransferFunctionPlots as tf
import TransFuncOpt as tf
import ClebGord as cg
from sympy.physics.quantum.spin import Rotation
import sympy as sp
import time
#import Spin2Harmonics as SWSH 

import numba as nb
from numba import float64, int32, jit, complex64, cfunc

Interpolate = 'ON'#OFF or ON

rec = -3.11293 #Conformal time at recombination

interpTable = np.load("ISWL2.npz")
kk = interpTable['kVals']
zz = interpTable['zvals']
isw = interpTable['isw']
ISWinterp = interpolate.RectBivariateSpline(kk,zz,isw)

@jit('float64(float64,float64)')
def interpTF(k,z):
    ISW = ISWinterp(k,z)[0][0]
    n1 = tf.ZtoN(z)
    SW = (bessel(2, k*(n1-rec)))#This is the Sachs-Wolfe term
    
    transFunc = -4*np.pi/3 * (SW + ISW)#note that the -ve sign is due to (-1j)**2 
    return transFunc 

@jit("int32(int32,float64)")
def bessel(l,q):
    """
    Parameters:
    q : product of k.del(r) 
    l : order of the spherical bessel function; sensible values = 0,2,4
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
        print("L must be 0, 2 or 4")
    

    #sphbessel = special.sph_jn(l,q)[0][l]
    #return sphbessel

@jit("float64(float64,int32,float64,float64,float64)")
def func(k,l,deltR,z1,z2):
    if(Interpolate == 'OFF'):
        return (   k**2 * powerSpec(k) * (tf.delta2(k,z1) * tf.delta2(k,z2))  * bessel(l,(k*deltR))   )
    else:
        return (  k**2 * powerSpec(k) * (interpTF(k,z1) * interpTF(k,z2))  * bessel(l,(k*deltR))   )

@jit("float64(int32,float64,float64,float64)")
def I(l,delR,z1,z2):
    """
    Parameters:
    k : Wavenumber
    l : order of spherical bessel function
    Vr1,Vr2: position vectors of galaxy clusters
    """    
    #print("l is", l, "delR is,", delR)
    
    deltR = delR[0]
    #if(deltR == 0):
    #    print('DeltaR is Zero')
    
    #func = lambda k: (k**2 * powerSpec(k) * (tf.delta2(k,z1) * tf.delta2(k,z2))  * bessel(l,(k*deltR)))
    #func = lambda k: (k**2 * powerSpec(k) * (tf.SWTerm(k,z1)**2)  * bessel(l,(k*deltR)))
    
    #plt.plot([func(k) for k in range(0,30)])
    #plt.show()
    
    
    temp = quad(func,0,500, args=(l,deltR,z1,z2), epsrel = 1e-5,epsabs = 1e-5)
    #for high z vals, use upper limit of 60.....
    #print("THE INTEGRAL IS,", temp[0])    
    return temp[0]

@jit("complex128(float64,int32,int32)")    
def spherHarm(delR,m,l):
    """
    Vr1,Vr2 : vectorzied spherical co-ordiantes
    """   
    theta = delR[1]
    phi = delR[2]     
    
    Y = special.sph_harm(m,l,phi,theta)
    
    return Y
    
def deltaR(Vr1,Vr2):
    """
    Calculates the resultant vector from Vr1-Vr2
    Note: This is currently sensitive to the relative ordering of Vr1 and Vr2. That's an unnecessary constraint
    and we are only interested in the relative distance (I think). TODO: Change this. 
    Paramters:
    Vr1, Vr2: Spherical Coordinates (vectors) of two galaxy glusters
    """
    #print("Vr1 and Vr2:", Vr1, Vr2)
    temp = SpherToCart(Vr1) - SpherToCart(Vr2)
    #resultant = CartToSpher(np.abs(temp))
    resultant = CartToSpher(temp)
    return resultant
    
def SpherToCart(position):
    """
    Parameters:
    position = Vectorized spherical co-ordinates, {r1,theta (polar angle), phi(coltaidunal angle)}
    
    Converts spherical to cartesian
    """
    r = position[0]
    theta = position[1]
    phi = position[2]
    
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    return np.array([x,y,z])

def CartToSpher(position):
    
    x = position[0]
    y = position[1]
    z = position[2]

    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y,x)
    theta = np.arctan2(np.sqrt(y**2 + x**2),z)

    return np.array([r,theta,phi])    
    
def varElem(m1, m2,delR,Il):#computes a2m-bar matrix elements
    
    M = m1-m2
        
    l = 0
    if(np.abs(M) > 2):
        l = 4
    
    elif(np.abs(M) > 0):
        l = 2
    
    
    a2m = 0    
    while(l<=4):#Sums the a2m's at l=0,2 and 4.
        Ylm =  np.conjugate(spherHarm(delR,M,l))#conjugated
        Ill = Il[l]
        j = cg.J(2,2,l,m1,m2,M)
        a2m  += (4*np.pi*((1j)**l)*Ill*Ylm*j)
        #print("m1:", m1, "m2:", m2, "M:",M,"L:",l,"Il:",Ill,"j:",j,"a2m:",a2m)        
        l+=2
    
    return a2m
    
    
def powerSpec(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**(3)


#Generates the <aBar2m1 aBar2m2> correlation matrix
def aBarCor(delR,Integrand):
    var = np.zeros((5,5)).astype(complex)
    for mj in range(-2,3):
        for mk in range(-2,3):
            #print("mj", mj)
            #print("mk",mk)
            var[mj+2][mk+2] = varElem(mj,mk,delR, Integrand)
            
    return var

def wignerRotation(V1,V2):
    z1 = V1[0]
    t1 = V1[1]
    p1 = V1[2]
    
    z2 = V2[0]    
    t2 = V2[1]
    p2 = V2[2]
    
    Vr1 = np.array([-tf.ZtoN(z1),t1,p1])    
    Vr2 =  np.array([-tf.ZtoN(z2),t2,p2])
    
    elDm1 = lambda m1, m2: sp.N(Rotation.D(2,m2,m1,-p1,t1,0).doit())#using physics convention: phi is the azimuth (in x-y plane)
    elDm2 = lambda m1, m2: sp.N(Rotation.D(2,m2,m1,-p2,t2,0).doit())
        
    delR = deltaR(Vr1,Vr2)
        
    Integrand = np.zeros(5)#store the values of Ils to save computational time
    for i in range(0,len(Integrand),2):
        Integrand[i] = I(i,delR,z1,z2)#toDo: pass in delR rather than Vr1 and Vr2;
        #same thing for aBarCor(Vr1 and Vr2 are only used in spherHarm)
    
    
    corelation = np.zeros((5,5)).astype(complex)    
    
    
    MatDm1 = np.zeros((5,5)).astype(complex)#The two wigner matrices
    MatDm2 = np.zeros((5,5)).astype(complex)
    
    for i in range(-2,3):
        for j in range(-2,3):
            MatDm1[i+2][j+2] = elDm1(i,j)
            MatDm2[i+2][j+2] = np.conjugate(elDm2(i,j))
            #MatDm1[i+2][j+2] = np.conjugate(elDm1(i,j))
            #MatDm2[i+2][j+2] = elDm2(i,j)
            
    var = aBarCor(delR,Integrand)#var is the <aBar2m1 aBar*2m2> matrix
    
    #print("var:\n", var, "\n MatDm1: \n", MatDm1,"\n MatDm2: \n", MatDm2)    
    
    #corelation = np.linalg.multi_dot([np.transpose(MatDm1),var,(MatDm2)])
    corelation = np.linalg.multi_dot([(MatDm1),var,np.transpose(MatDm2)])
    
    #print("Before rotation (A2m bar):\n:", var, "\n After Wigner rotation(A2m):\n", corelation)
    
    return corelation, var, MatDm1, MatDm2

"""
#clusters = np.load("rand100Clusters.npy")
clusters=np.load("VoxList2.npy")
n = 10
for i in range(0,n):
    for j in range(0,i+1):
      temp = wignerRotation(clusters[i],clusters[j])
      print(temp[0][0][0])
"""
"""
thetalist = np.arange(0,np.pi,0.2)
philist = np.arange(0,2*np.pi,0.4)

temp = [wignerRotation([1,0,0],[1,1,i])[0][0][0] for i in philist]
"""
"""
V1 = np.array([1,0.1,0.3])    
V0 = np.array([0,0,0])

temp = wignerRotation(V1, V0)[0][0][0]
print(temp)
"""

"""
#v1 = np.array([1,1.1,6.2])
#v2 = np.array([3,0.45,0.97])
v3 = np.array([1,4.5,0.69])
v4 = np.array([1,0,0])
v0 = np.array([0,0,0])
t0 = time.time()
temp=wignerRotation(v0,v0)[0][0][0]
t1 = time.time()
print(temp)
print("Time this process took:", t1-t0)
"""
"""

z = 0.025
m = -2
Vl = np.array([0,0,0])
V1 = np.array([z,0,0])
V2 = np.array([z,1.1,2.2])
V3 = np.array([z,2.3,0.69])
#Vr1 = np.array([-tf.ZtoN(V1[0]),V1[1],V1[2]])    
t0 = time.time()
#temp = I(2,2,Vr1, 0, powerSpec)
#temp = aBarCor(Vr1,3,powerSpec,V1[0])
f1 = wignerRotation(V1, Vl)[0][0][m+2]
s1 = SWSH.sYlm(2,2,m,V1[2],V1[1])
f2 = wignerRotation(V2, Vl)[0][0][m+2]
s2 = SWSH.sYlm(2,2,m,V2[2],V2[1])
f3 = wignerRotation(V3, Vl)[0][0][m+2]
s3 = SWSH.sYlm(2,2,m,V3[2],V3[1])
print("VarMatGeneralized.py.....","zeta1:", f1/s1,"zeta2:", f2/s2, "zeta3:", f3/s3)

t1 = time.time()
print("Time this process took: ", (t1-t0))
"""

"""
Vl = np.array([0,0,0])
V1 = np.array([1,0.1,0.3])
V2 = np.array([1,1.1,2.2])
#Vr1 = np.array([-tf.ZtoN(V1[0]),V1[1],V1[2]])    
t0 = time.time()
#temp = I(2,2,Vr1, 0, powerSpec)
#temp = aBarCor(Vr1,3,powerSpec,V1[0])
f1 = wignerRotation(V1,Vl)[0][0][0]
s1 = SWSH.sYlm(2,2,-2,0.3,0.1)
f2 = wignerRotation(V2, Vl)[0][0][0]
s2 = SWSH.sYlm(2,2,-2,2.2,1.1)
t1 = time.time()
print("Time this process took: ", (t1-t0))
print("zeta1:", f1/s1,"zeta2:", f2/s2)
"""
"""
n=30
zlist = np.logspace(-3,1,n)
clusterList = np.zeros((n,3))
for z in range(0,n):
    clusterList[z] = np.array([zlist[z],0,0])

c2 = np.zeros(n)
for i in range(0,n):
    c2[i] =wignerRotation(clusterList[i,:],clusterList[i,:])[0][0][0]#this is wignerRotation from VarMatrix_Generalized.py

plt.figure(2)
print("This window plots C2(z) vs Z")
plt.semilogx(zlist, c2)
plt.title('total(n=1)')
#plt.plot(zlist,c2)
"""
"""    
n=100
clusterList2 = np.zeros((n,3))
for z in range(0,n):
    clusterList2[z] = np.array([z/1000,0,0])
    
c22 = np.zeros(n)
for i in range(0,n):
    c22[i] = wignerRotation(clusterList2[i,:],clusterList2[i,:])[0][0][0]#this is wignerRotation from VarMatrix_Generalized.py

print("This window plots C2(z) vs Z")
plt.figure(6)
plt.semilogx(np.arange(0,n/1000,0.001), c22,'r')
"""
"""
V1 = np.array([0,0,0])
#V2 = np.array([2,0,0])
#V3 = np.array([0.1,1.2,0.7])
V3 = np.array([2,1.1,2.3])
t0 = time.time()
temp = wignerRotation(V1,V1)
t1 = time.time()
print("Time this process took: ", (t1-t0))
print(temp[0][0][0])
"""

"""
thetalist = np.arange(50)/10
toPlot = [wignerRotation( np.array([0,0,0] ),np.array([thetalist[i],0,0]))[0][0][0]/1.89 for i in range(thetalist.shape[0])]
#1.43 is the normalization!!!!
plt.figure()
plt.plot(thetalist,toPlot)
plt.show()
"""
"""
#Let's make some 3D plots now;

plotList = np.zeros((20,20))
zlist = np.arange(20)*0.1
thetalist = np.arange(20)*np.pi/20
for i in range(0,20):
    for j in range(0,20):
        plotList[i][j] = wignerRotation(np.array([2,0,0]),np.array([zlist[i],0,thetalist[j]]))[0][0][0]/1.43
        print(plotList[i][j])
        
            
X, Y = np.meshgrid(zlist, thetalist)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X,Y,plotList)
"""


"""
var = np.zeros((5,5)).astype(complex)
for mj in range(-2,3):
    for mk in range(-2,3):
        var[mj+2][mk+2] = varElem(mj,mk,p,2,2,0,0,0,0, Integrand)
        
        #print("Var element: ", var[mk][mj])
    #print(var[mj][:])
    
print("Variance Matrix: \n", var)
"""

