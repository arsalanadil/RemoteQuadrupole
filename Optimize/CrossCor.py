# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:46:01 2016

This routine generates cross corelation (between local multipoles and remote quadrupoles) covariance and relation matrices. It is 
mostly based off VarMatrix_Generalized.py with the necessary mathematical tweaks. Meant to be called from COVandREL.py.
Note that we are using the generalized Transfer Funciton routine for this.

@author: arsalanadil
"""


import scipy.special as special, scipy.interpolate as interpolate
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import TransFuncOpt as tf
import ClebGord as cg
from sympy.physics.quantum.spin import Rotation
import sympy as sp
import time
import Spin2Harmonics as SWSH
from numba import float64, int32, jit, complex64, float32, cfunc

#import FisherInfo_toy_mdoel as fi

Interpolate = 'on' #interpolate the transfunction? Off if Interpolate = 'OFF'. Else, it is ON.

rec = -3.11293 #Conformal time at recombination

interpTable2 = np.load("ISWL2.npz")
kk = interpTable2['kVals']
zz = interpTable2['zvals']
isw2 = interpTable2['isw']

interpTable20 = np.load("ISWL2Z0.npz")
interpTable3 = np.load("ISWL3Z0.npz")
interpTable4 = np.load("ISWL4Z0.npz")
interpTable5 = np.load("ISWL5Z0.npz")
isw20 = interpTable20['isw']
isw3 = interpTable3['isw']
isw4 = interpTable4['isw']
isw5 = interpTable5['isw']
kk0 = interpTable3['kVals']

ISWinterp2 = interpolate.RectBivariateSpline(kk,zz,isw2)
ISWinterp20 = interpolate.interp1d(kk0,isw20)
ISWinterp3 = interpolate.interp1d(kk0,isw3)
ISWinterp4 = interpolate.interp1d(kk0,isw4)
ISWinterp5 = interpolate.interp1d(kk0,isw5)

@jit(float64(int32,float64,float64))
def interpolation(l,k,z):
    if(l==2):
        if(z==0):
            return ISWinterp20(k)
        else:
            return ISWinterp2(k,z)[0][0]

    if(l==3):
        return ISWinterp3(k)
    
    if(l==4):
        return ISWinterp4(k)
    
    if(l==5):
        return ISWinterp5(k)

    else:
        print("Transfer Function only defined L=2,3,4,5!")

@jit(float64(int32,float64,float64))
def interpTF(l,k,z):
    ISW = interpolation(l,k,z)
    n1 = tf.ZtoN(z)
    SW = (tf.jl(l, k*(n1-rec)))#This is the Sachs-Wolfe term
    
    transFunc = 4*np.pi/3 * (SW + ISW)#note that the -ve sign is due to (-1j)**2 
    
    #if(l==2 and z!=0):
    #    return -transFunc
    
                            
    return transFunc 

    

@jit(float64(int32,float64))
def bessel(l,q):
    """
    Parameters:
    q : product of k.del(r) 
    l : order of the spherical bessel function; L=l, l-1, l+2
    """
    sphbessel = special.sph_jn(l,q)[0][l]
    return sphbessel

@jit(float64(float64,int32,int32,float64,float64))
def func(k, L,l,vr1,z1):
    if(Interpolate == 'OFF'):
        return (   (k**2 * powerSpec(k) * tf.delta2(k,z1) * 
                        tf.delta(l,k,0) * bessel(L,(k*vr1)))    )
    
    return (   (k**2 * powerSpec(k) * -1*interpTF(2,k,z1) * 
                        interpTF(l,k,0) * bessel(L,(k*vr1)))    )

@jit(float64(int32,int32,float64,float64))       
def KIntegral(L,l,Vr1,z1):
    """
    Parameters:
    delR: co-ordinates (r,theta,phi) of the difference between cluster1 and cluster2
    L : order of spherical bessel function (that needs to be summed over)
    l : local multipole
    z1, z2: red shift of clusters
    """    
        
    vr1 = Vr1[0]   
    
    #nint = lambda k,L,l,vr1,z1: cfunc('float32(int32,int32,float32,float32)')(func)
    
    temp = quad(func,0,250, args=(L,l,vr1,z1), epsrel = 1e-5,epsabs = 1e-5)
    
    #print("l is", l, "L is, ", L, "Integrand is, ", temp[0])    
    
    return temp[0]
    
def spherHarm(Vr1,m,l):
    """
    V1: Co-ordinates of the difference between two clusters in SHERICAL co-ordinates
    m,l: whatever m and l are for the spherical harmonics (Ylm)
    """   
    theta = Vr1[1]
    phi = Vr1[2]     
    
    #Y = special.sph_harm(m,l,theta,phi)
    Y = special.sph_harm(m,l,phi,theta)#m,l,phi,theta
    return Y

    
def varElem(m1, m2,Vr1,l,z1,Kint):#computes a2m-bar matrix elements
    
    M = m1-m2
    #M = m2-m1    
    L = l-2
    j = 0
    """  
    l = 0
    if(np.abs(M) > 2):
        l = 4
    
    elif(np.abs(M) > 0):
        l = 2
    """
    i = 0   
    plm = 0    
    while(L<=l+2):#Sums the a2m's at L=0,2 and 4.        
        if (np.abs(M)>L):
            plm += 0
            L+=2
            i+=1
            continue
        
        j = cg.J(2,l,L,m1,m2,M)
        #print("Clebsch Gordon,",j)
        #note that the parameters are passed differently than in VarMatrix_Generalized. This is to take into account 
        #conjugation of the REMOTE part, i.e. <bar(a2m1) alm2(0)*>. See pg 36 in notes.
        
                                                  
        Ylm =  np.conjugate(spherHarm(Vr1,M,L))#conjugated               
        
        #Ill = (-1j)**(-l) * KIntegral(L,l,Vr1,z1)
        Ill = Kint[i]
        #print("L=",L,"Ill=",Ill)
         
        plm  += (4*np.pi*((1j)**L)*Ylm*Ill*j)
        
        L+=2
        i+=1
        
    return plm
    
    


#Generates the <aBar2m1 alm> correlation matrix
def aBarCor(Vr1,l,z1):
    length = 2*l + 1
    var = np.zeros((5,length)).astype(complex)
    
    Kint = np.zeros(3).astype(complex)
    L = l-2
    i = 0
    while(L<=l+2):
        Kint[i] = (-1j)**(-l) * KIntegral(L,l,Vr1,z1)
        i+=1
        L+=2
        
    for mj in range(-2,3):
        for mk in range(-l,l+ 1):
            #print("mj",mj)
            #print("mk", mk)
            var[mj+2][mk+l] =  varElem(mj,mk,Vr1, l,z1, Kint)
            #print(var)
    #print("unrotated var matrix for l =", l, "at cluster location", Vr1, var)
    return var



def wignerRotation(V1, l):
    z1 = V1[0]
    t1 = V1[1]
    p1 = V1[2]
    
    Vr1 = np.array([-tf.ZtoN(z1),t1,p1])    
    
    elDm1 = lambda m1, m2: Rotation.D(2,m2,m1,-p1,t1,0).doit()
    
    
 
    var = aBarCor(Vr1,l,V1[0])#var is the <aBar2m1 aBar*2m2> matrix


    
    MatDm1 = np.zeros((5,5)).astype(complex)#The wigner matrix
    
    for i in range(-2,3):
        for j in range(-2,3):
            MatDm1[i+2][j+2] = elDm1(i,j)
    
    
    
    corelation = np.matmul(MatDm1, var)
    return corelation, var, MatDm1
    



def powerSpec(k):
    if(k==0):
        return 0
    else:
        return 1/k**3
    
"""  
V1 = np.array([0,0,0])    
temp = wignerRotation(V1, 2)[0][0][0]
print(temp)
"""
"""    
    
#V1 = np.array([1,0.1,0.3])
#V2 = np.array([1,1.1,2.2])
#V3 = np.array([1,4.5,0.69])
z = 1
l=5
m=-2#m must be less than l
V1 = np.array([z,0,0])
V2 = np.array([z,1.1,2.2])
V3 = np.array([z,2.3,0.69])
t0 = time.time()
#temp = I(2,2,Vr1, 0, powerSpec)
#temp = aBarCor(Vr1,3,powerSpec,V1[0])
#f1 = wignerRotation(V1, l)[0][0][(m+l)]
#s1 = SWSH.sYlm(2,l,m,V1[2],V1[1])
f2 = wignerRotation(V2, l)[0][0][(m+l)]
print(f2)
s2 = SWSH.sYlm(2,l,m,V2[2],V2[1])
f3 = wignerRotation(V3, l)[0][0][(m+l)]
s3 = SWSH.sYlm(2,l,m,V3[2],V3[1])
print("VarMatrixGenLocal.py.....","zeta1:", "zeta2:", f2/s2, "zeta3:",f3/s3)
#print(f2/s2)
t1 = time.time()
print("Time this process took: ", (t1-t0))
#print(temp[0])
#print(f1)
"""


"""
v1 = np.array([0.4,0,0])
toPlot = [np.abs(wignerRotation(v1,i,powerSpec)[0][0][0])/1.43 for i in range(2,5)]
plt.plot(np.array([2,3,4]), toPlot)
print(toPlot)
"""
"""
zlist = np.arange(20)/10
toPlot = [wignerRotation(np.array([thetalist[i],0,0]),2,powerSpec )[0][0][0]/1.85 for i in range(zlist.shape[0])]
plt.figure()
plt.plot(zlist,toPlot)
plt.show()
"""

"""
thetalist = np.arange(20)*np.pi/10
toPlot = [wignerRotation(np.array([2,0,np.pi/4]),np.array([2,0,thetalist[i]]))[0][0][0] for i in range(20)]
plt.figure()
plt.plot(thetalist,toPlot)
plt.show()
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
