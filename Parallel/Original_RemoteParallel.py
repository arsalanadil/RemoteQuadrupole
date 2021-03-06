#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 12:25:31 2017
For archiving puprposes. This version has been depricated and is no longer in use. Use the new verion, RemoteParallel.py
@author: arsalanadil
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Created on  Jul  17 2017
    
    @author: arsalanadil
"""

from multiprocessing import Pool
import time
import numpy as np
from scipy.integrate import quad
import scipy.special as special, scipy.interpolate as interpolate
from sympy.physics.quantum.spin import Rotation
import sympy as sp
from scipy.special import binom
import sys
import numba as nb
from numba import float64, int32, jit, complex64, cfunc

from pympler.tracker import SummaryTracker
tracker = SummaryTracker()

rec = -3.11293 #Conformal time at recombination

interpTable = np.load("ISWL2.npz")
kk = interpTable['kVals']
zz = interpTable['zvals']
isw = interpTable['isw']
ISWinterp = interpolate.RectBivariateSpline(kk,zz,isw)

qlist = np.loadtxt("qList1.csv", delimiter = ',')

zq = qlist[:,0]
qq = qlist[:,1]

nlist = np.loadtxt("nList1.csv", delimiter = ',')
zn = nlist[:,0]
nn = nlist[:,1]#array of eta values; zn is corresponding z vals

#@jit("float64(float64)")              
def ZtoN(z):
    n = np.interp(z, zn, nn)
    return n

#@jit("float64(int32, int32, int32, int32, int32, int32)")              
def clebgor(j1, j2, j, m1, m2, m):
    """
        Parameters:    j1, j2, j: the angular momenta input
        m1, m2, m: the z components of angular momenta input
        Returns:       The numerical value of the Clebsch-Gordan coeffcient.
        Remarks:       Note that in the sum none of the binomial coeffcients
        can have negative values.  Thus, zmin is there to make
        sure that the sums have a cut-off.
    """
    zmin = int(min([j1 - m1, j2 + m2]))
    J = j1 + j2 + j
    return (int(m1 + m2 == m) *
            int(np.abs(j1 - j2) <= j <= (j1 + j2)) *
            int(np.abs(m1) <= j1) *
            int(np.abs(m2) <= j2) *
            int(np.abs(m) <= j) *
            int((j1 + m1) >= 0.0) *
            int((j2 + m2) >= 0.0) *
            int((j + m) >= 0.0) *
            int(J >= 0) *
            np.sqrt(binom(2 * j1, J - 2 * j) *
                    binom(2 * j2, J - 2 * j) /
                    (binom(J + 1, J - 2 * j) *
                     binom(2 * j1, j1 - m1) *
                     binom(2 * j2, j2 - m2) *
                     binom(2 * j, j - m))) *
            np.sum([(-1) ** z * binom(J - 2 * j, z) *
                    binom(J - 2 * j2, j1 - m1 - z) *
                    binom(J - 2 * j1, j2 + m2 - z) for z in range(zmin + 1)]))

#@jit("float64(int32, int32, int32, int32, int32, int32)")              
def cleb(j1,j2,j,m1,m2,m):
    if((j1 + j2 + j)%2 != 0):
        return 0

    temp = j1#Switch {j1,m1} with {j,m} to follow convention
    temp1 = m1
    j1 = j
    m1 = m
    j = temp
    m = temp1
    
    """
        Now {j1,m1} are the parameters of the first Y2m^* term and {j,m} are conventionally
        defined to correspond to the YLM term (the so called "total angular momentum")
    """

    cg = clebgor(j1, j2, j, m1, m2, m)
    cg0 = clebgor(j1,j2,j,0,0,0)
    val = (np.sqrt(((1 + 2*j1)*(1 + 2*j2))/((1 + 2*j)*4*np.pi)))*cg*cg0
    return val

#@jit('float64(float64,float64)')
def interpTF(k,z):
    ISW = ISWinterp(k,z)[0][0]
    n1 = ZtoN(z)
    SW = (bessel(2, k*(n1-rec)))#This is the Sachs-Wolfe term
    
    transFunc = -4*np.pi/3 * (SW + ISW)#note that the -ve sign is due to (-1j)**2
    return transFunc

#@jit("int32(int32,float64)")
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


@jit("float64(float64,int32,float64,float64,float64)")
def func(k,l,deltR,z1,z2):
    return (   k**2 * powerSpec(k) * (interpTF(k,z1) * interpTF(k,z2))  * bessel(l,(k*deltR))   )

@jit("float64(int32,float64,float64,float64)")
def I(l,delR,z1,z2):
    """
        Parameters:
        k : Wavenumber
        l : order of spherical bessel function
        Vr1,Vr2: position vectors of galaxy clusters
        """
    
    deltR = delR[0]
    
    
    temp = quad(func,0,300, args=(l,deltR,z1,z2), epsrel = 1e-5,epsabs = 1e-5)
    return temp[0]

#@jit("complex128(float64,int32,int32)")    
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
    temp = SpherToCart(Vr1) - SpherToCart(Vr2)
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
        j = cleb(2,2,l,m1,m2,M)
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
    
    Vr1 = np.array([-ZtoN(z1),t1,p1])
    Vr2 =  np.array([-ZtoN(z2),t2,p2])
    
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

    var = aBarCor(delR,Integrand)#var is the <aBar2m1 aBar*2m2> matrix
    
    
    #corelation = np.linalg.multi_dot([np.transpose(MatDm1),var,(MatDm2)])
    corelation = np.linalg.multi_dot([(MatDm1),var,np.transpose(MatDm2)])
    
    #print("Before rotation (A2m bar):\n:", var, "\n After Wigner rotation(A2m):\n", corelation)
    
    return corelation, var, MatDm1, MatDm2


def paralFunc(v1,v2):
    #print(v1,v2)
    return wignerRotation(v1,v2)[0]
"""
n =200
clusters = np.load("VoxList2.npy")

t0 = time.time()
temp = RemoteCMB(clusters,n)
t1 = time.time()
print("Time taken:", t1-t0)
segments = 1
"""

n = 6912#number of processes
#n = 12
segments = int(sys.argv[2])
clusters = np.load("VoxList2.npy")
#clusters = np.load("rand100Clusters.npy")

import EqualArea as EA
nVals = EA.divisions(n,segments)#generated using EqualArea.py to divide covariance matrix into equal area segments
def parallelScatter(number):
    nValss = np.concatenate([nVals,[n]])
    return nValss[number]
    
"""
def RemoteCMB(clusters, n):

    varMatrix = np.zeros((n,n)).astype(complex)
    relMatrix = np.zeros((n,n)).astype(complex)
    
    for i in range(n):
        for j in range(i+1):
            #print("For Remote: Cluster number", i, "and",j)
            temp = wignerRotation(clusters[i,:],clusters[j,:])[0]
            #print("For Remote: Cluster number", i, "and",j,".Covariance",temp[0][0])
            #print("wignerRotation() returned:", temp[0])        
            varMatrix[i][j] = temp[0][0]
            varMatrix[j][i] = np.conjugate(temp[0][0])
            relMatrix[i][j] = temp[0][4]
            relMatrix[j][i] = temp[0][4]
    
    return varMatrix, relMatrix
"""

"""
segment = 'ON'
if __name__ == '__main__':
    if(segment == 'ON'):
        number = int(sys.argv[1])
        nStart = int(parallelScatter(number-1))
        nEnd = int(parallelScatter(number))
        if(number == 0):
            nEnd = int(parallelScatter(number))
            nStart = 0
        
        print(nStart, nEnd)
        varMatrix = np.zeros((n,n)).astype(complex)
        relMatrix = np.zeros((n,n)).astype(complex)
        
        t0 = time.time()
        with Pool(7) as p:
            temp = (p.starmap(paralFunc, [(clusters[i], clusters[j]) for i in range(nStart,nEnd) for j in range(0,i+1)] ))
        t1 = time.time()
        print("TIme Taken:", t1-t0)
        for k in range(0,len(temp)):
            print(temp[k][0][0])
        
        k = 0
        for i in range(nStart,nEnd):
            for j in range(0,i+1):
                varMatrix[i][j] = temp[k][0][0]#this should probably be [j][i] not [i][j]
                varMatrix[j][i] = np.conjugate(temp[k][0][0])
                #print(varMatrix[i][j])
                relMatrix[i][j] = temp[k][0][4]
                relMatrix[j][i] = temp[k][0][4]
                k+=1
                #np.savez("results%s"%number,varMatrix,relMatrix)
        #print(varMatrix, varMatrix.shape)
    else:
        print("NOT SEGMENTING")
        varMatrix = np.zeros((n,n)).astype(complex)
        relMatrix = np.zeros((n,n)).astype(complex)
        
        t0 = time.time()
        with Pool() as p:
            temp = (p.starmap(paralFunc, [(clusters[i], clusters[j]) for i in range(0,n) for j in range(0,i+1)] ))
        t1 = time.time()
        print("TIme Taken:", t1-t0)
        k = 0
        for i in range(0,n):
            for j in range(0,i+1):
                varMatrix[i][j] = temp[k][0][0]#this should probably be [j][i] not [i][j]
                varMatrix[j][i] = np.conjugate(temp[k][0][0])
                #print(varMatrix[i][j])
                relMatrix[i][j] = temp[k][0][4]
                relMatrix[j][i] = temp[k][0][4]
                k+=1
"""


if __name__ == '__main__':
    number = int(sys.argv[1])
    nStart = int(parallelScatter(number-1))
    nEnd = int(parallelScatter(number))
    if(number == 0):
        nEnd = int(parallelScatter(number))
        nStart = 0
    
    #if(number == segments):
    #    nEnd = n+1
        
    print(nStart, nEnd)
    
    #varMatrix = np.zeros((n,n)).astype(complex)
    #relMatrix = np.zeros((n,n)).astype(complex)
        
    t0 = time.time()
    with Pool() as p:
        temp = (p.starmap(paralFunc, [(clusters[i], clusters[j]) for i in range(nStart,nEnd) for j in range(0,i+1)] ))
    t1 = time.time()
    print("TIme Taken:", t1-t0)
    #k = 0
    #length = nEnd - nStart
    length = len(temp)
    varTemp = np.zeros(length).astype(complex)
    relTemp = np.zeros(length).astype(complex)
    for i in range(0,len(temp)):
        print(temp[i][0][0])
    
    for i in range(0,length):
        #print(k)
        varTemp[i] = temp[i][0][0]#this should probably be [j][i] not [i][j]
        #varMatrix[j][i] = np.conjugate(temp[k][0][0])
        #print(varMatrix[i][j])
        relTemp[i] = temp[i][0][4]
        #relMatrix[j][i] = temp[k][0][4]
        #k+=1
    np.savez("results%s"%number,varTemp,relTemp)
    print(varTemp)
        


tracker.print_diff()