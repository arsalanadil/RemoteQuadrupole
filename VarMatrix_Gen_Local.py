# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:46:01 2016

This routine generates cross corelation (between local multipoles and remote quadrupoles) covariance and relation matrices. It is 
mostly based off VarMatrix_Generalized.py with the necessary mathematical tweaks. Meant to be called from COVandREL.py.
Note that we are using the generalized Transfer Funciton routine for this.

@author: arsalanadil
"""


import scipy.special as special
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad
import TransFuncGen as tf
import TransferFunction as tf2
import ClebGord as cg
from sympy.physics.quantum.spin import Rotation
import sympy as sp
import time
import Spin2Harmonics as SWSH 

#import FisherInfo_toy_mdoel as fi

def bessel(l,q):
    """
    Parameters:
    q : product of k.del(r) 
    l : order of the spherical bessel function; L=l, l-1, l+2
    """
    sphbessel = special.sph_jn(l,q)[0][l]
    return sphbessel
        
def KIntegral(L,l,Vr1,z1,ps):
    """
    Parameters:
    delR: co-ordinates (r,theta,phi) of the difference between cluster1 and cluster2
    L : order of spherical bessel function (that needs to be summed over)
    l : local multipole
    z1, z2: red shift of clusters
    """    
        
    #deltR = delR[0]
    
    func = lambda k: (   (k**2 * ps(k) * tf2.delta(k,z1)[0] * 
                         tf.delta(l,k,0)[0] * bessel(L,(k*Vr1[0])))    )
    
    #func = lambda k: (   k**2 * ps(k) * bessel(L,(k*Vr1[0])) )    
    
    #tf.delta is the transfer function
    
    #for k in range(1,20):
    #    print(bessel(L,(k*Vr1[0])) )
    #plt.plot(np.arange(0,30),[func(k) for k in range(0,30)])
    #plt.show()
    
    
    
    temp = quad(func,0,100, limit=10, epsrel = 1e-5,epsabs = 1e-5)
    
    #print("l is", l, "L is, ", L, "Integrand is, ", temp[0])
    #choose an upper limit for integration carefully. This function dies down really fast
    #so k=10 is not bad.
    print(temp[0])
    return temp
    
def spherHarm(Vr1,m,l):
    """
    V1: Co-ordinates of the difference between two clusters in SHERICAL co-ordinates
    m,l: whatever m and l are for the spherical harmonics (Ylm)
    """   
    theta = Vr1[1]
    phi = Vr1[2]     
    
    Y = special.sph_harm(m,l,theta,phi)
    
    return Y

    
def varElem(m1, m2,Vr1,l,z1,ps):#computes a2m-bar matrix elements
    
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
   
    #this part is a little tricky.... need to ensure correctness of the factors of i
    
    plm = 0    
    while(L<=l+2):#Sums the a2m's at L=0,2 and 4.        
        if (np.abs(M)>L):
            plm += 0
            L+=2
            continue
        
        #j = cg.J(l,2,L,m2,m1,M)
        j = cg.J(2,l,L,m1,m2,M)
        #note that the parameters are passed differently than in VarMatrix_Generalized. This is to take into account 
        #conjugation of the REMOTE part, i.e. <bar(a2m1) alm2(0)*>. See pg 36 in notes.
        
                                                  
        Ylm =  np.conjugate(spherHarm(Vr1,M,L))#conjugated               
        #Ylm =  spherHarm(Vr1,M,L)
        #Ill = np.conjugate(((1j)**l) * KIntegral(L,l,Vr1,z1,ps)[0]   ) #the (-1j)**l term is correction for transfer function
        Ill = (-1j)**(-l) * KIntegral(L,l,Vr1,z1,ps)[0]
        
        
        #Ill = (-1j)**(-l) * KIntegral(L,l,Vr1,z1,ps)[0]
        #Ill =  ((1j)**l) * KIntegral(L,l,Vr1,z1,ps)[0]              
        #plm  += (-1)**m1 * (4*np.pi*((1j)**L)*Ylm*Ill*j)
        #print("m1:", m1, "m2:", m2, "M:",M,"L:",L,"Il:",Ill,"j:",j,"a2m:",plm)        
        plm  += (4*np.pi*((1j)**L)*Ylm*Ill*j)
        
        L+=2
    
    return plm
    
    


#Generates the <aBar2m1 alm> correlation matrix
def aBarCor(Vr1,l,ps,z1):
    length = 2*l + 1
    var = np.zeros((5,length)).astype(complex)
    for mj in range(-2,3):
        for mk in range(-l,l+ 1):
            #print("mj",mj)
            #print("mk", mk)
            var[mj+2][mk+l] = varElem(mj,mk,Vr1, l,z1,ps)
            #print(var)
    #print("unrotated var matrix for l =", l, "at cluster location", Vr1, var)
    return var



def wignerRotation(V1, l, powerSpec):
    z1 = V1[0]
    t1 = V1[1]
    p1 = V1[2]
    
    Vr1 = np.array([-tf.ZtoN(z1),t1,p1])    
    
    elDm1 = lambda m1, m2: Rotation.D(2,m2,m1,-t1,p1,0).doit()
    #elDm2 = lambda m1,m2: Rotation.D(2,m2,m1,-t2,p2,0).doit()
    
        
    
    #corelation = np.zeros((5,5)).astype(complex)    
 
    var = aBarCor(Vr1,l,powerSpec,V1[0])#var is the <aBar2m1 aBar*2m2> matrix


    
    MatDm1 = np.zeros((5,5)).astype(complex)#The wigner matrix
    #MatDm2 = np.zeros((5,5)).astype(complex)
    
    for i in range(-2,3):
        for j in range(-2,3):
            MatDm1[i+2][j+2] = elDm1(i,j)
            #MatDm2[i+2][j+2] = np.conjugate(elDm2(i,j))
    
    
    #print("var:\n", var, "\n MatDm1: \n", MatDm1,"\n MatDm2: \n", MatDm2)    
    
    #corelation = np.linalg.multi_dot([MatDm1,var,np.transpose((MatDm2))])#matrix after rotation by the respective 
    #wignerD matrices
    
    #print("Before rotation (A2m bar):\n:", var, "\n After Wigner rotation(A2m):\n", corelation)
    corelation = np.matmul(MatDm1, var)
    return corelation, var, MatDm1
    



def powerSpec(k):
    if(k==0):
        return 0
    else:
        return 1/k**3
"""
V3 = np.array([1,4.5,0.69])
#Vr1 = np.array([-tf.ZtoN(V1[0]),V1[1],V1[2]])    
t0 = time.time()
#temp = I(2,2,Vr1, 0, powerSpec)
#temp = aBarCor(Vr1,3,powerSpec,V1[0])
f1 = wignerRotation(V3, 2, powerSpec)[0][0][0]

#s1 = SWSH.sYlm(2,2,-2,V1[1],V1[0])
#f2 = wignerRotation(V2, 2, powerSpec)[0][0][0]
#s2 = SWSH.sYlm(2,2,-2,V2[1],V2[0])
#f3 = wignerRotation(V3, 2, powerSpec)[0][0][0]
#s3 = SWSH.sYlm(2,2,-2,V3[1],V3[0])
#print("VarMatrixGenLocal.py.....","zeta1:", f1/s1,"zeta2:", f2/s2, "zeta3:",f3/s3)
 
t1 = time.time()
print("Time this process took: ", (t1-t0))
#print(temp[0])
print(f1)
"""

    
"""    
#V1 = np.array([0.4,0,0])
    
#temp = wignerRotation( p, 2,2,np.pi/3,np.pi/3+0.3,np.pi/4+.3,np.pi/2)

#V1 = np.array([2,1.1,3.2])
#V2 = np.array([3,1.7,0.7])
#V1 = np.array([0.1,1.2,0.7])
V1 = np.array([1,0.1,0.3])
V2 = np.array([1,1.1,2.2])
V3 = np.array([1,4.5,0.69])
#Vr1 = np.array([-tf.ZtoN(V1[0]),V1[1],V1[2]])    
t0 = time.time()
#temp = I(2,2,Vr1, 0, powerSpec)
#temp = aBarCor(Vr1,3,powerSpec,V1[0])
f1 = wignerRotation(V1, 2, powerSpec)[0][0][0]
s1 = SWSH.sYlm(2,2,-2,V1[1],V1[0])
f2 = wignerRotation(V2, 2, powerSpec)[0][0][0]
s2 = SWSH.sYlm(2,2,-2,V2[1],V2[0])
f3 = wignerRotation(V3, 2, powerSpec)[0][0][0]
s3 = SWSH.sYlm(2,2,-2,V3[1],V3[0])
print("VarMatrixGenLocal.py.....","zeta1:", f1/s1,"zeta2:", f2/s2, "zeta3:",f3/s3)
 
t1 = time.time()
print("Time this process took: ", (t1-t0))
#print(temp[0])

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