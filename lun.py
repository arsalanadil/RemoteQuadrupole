# -*- coding: utf-8 -*-
"""
Created on Fri May 13 19:32:14 2016

@author: arsalanadil
"""
import numpy as np, matplotlib.pyplot as plt, spherepy as sp
from scipy.integrate import dblquad,quad
import scipy.special as spec
import parser
import scipy
from sympy.physics.quantum.spin import Rotation
from sympy import *
#import sympy as sp
#from sympy.abc import m
#import VarMatrix_Generalized as vm

"""
plt.figure(2)
x = np.arange(1, 200, 1)
y = np.exp(-(x/5)**2)
plt.plot(x,y)
plt.plot(x,y)
ax = plt.gca()
#ax.set_yscale('log')    
ax.set_xscale('log')
plt.title("1/exp(-k/5)^2")
plt.show()    
"""
"""
C = sp.zeros_coefs(5,5)
C[2,0] = 1
p = sp.ispht(C,50,50)
sp.plot_sphere_mag(p)
"""
"""
a = sp.zeros_coefs(3,3)
a[2,0] = 1

map = sp.ispht(a,256,256)
sp.plot_sphere_mag(map)
"""

"""
n=1000
g = np.zeros(n)

for x in range(n):
    g[x] = np.random.normal(0,1)
        
plt.hist(g)
plt.show()
"""
"""
mean = [0, 0]
cov = [[1, 0], [0, 100]]  # diagonal covariance
x, y = np.random.multivariate_normal(mean, cov, 5000).T
plt.plot(x, y, 'x')
plt.axis('equal')
plt.show()
"""
"""
A = np.matrix([[1,2],[3,4]])
print(A)

x = np.vstack(np.array([2,3]))
print(x)

#b = A*np.vstack(x)
c = np.dot(A,x)
print(A*x)
"""

"""
def covarianceMatrix(n):#Generates an nxn random positive definite (covariance) matrix    

    #A = np.zeros((n,n))
   
    #np.fill_diagonal(A, np.random.normal(0,1))    
    #np.fill_diagonal(A, np.random.randint(1,5))
        
    A = np.random.randint(1,100,size=(n,n))
    #A = A * np.matrix.transpose(A)
    A = np.matmul(A,np.matrix.transpose(A))    
    print(A)        
    print(np.linalg.eigvals(A))
    np.linalg.cholesky(A)    
        
covarianceMatrix(6)
"""
"""
def covarianceMatrix(n):#Generates an nxn random positive definite (covariance) matrix    

    A = np.zeros((n,n))
   
    #np.fill_diagonal(A, np.random.normal(0,1))    
    #np.fill_diagonal(A, np.random.randint(1,5))
        
    A = np.random.random(size=(n,n))
    A = np.matmul(A,np.matrix.transpose(A))    
    #print(A)        
    return A
    
cov = covarianceMatrix(6)
    
n = len(cov)
m=5
x= np.zeros((m,n))
   
for i in range(0,3):  
    x[i,:] = np.random.multivariate_normal(np.zeros(n), cov)

#print(np.average(x[1]))
print(x[:,0], "\n")
#print(np.average(x[2]))

#plt.hist(x[:,[1]])
#plt.show()
    
print(x)
"""

"""
def gaussVecGen(cov, m):#computes Cholesky decomposition and generates m n-dimenstional vectors 
    
    n = len(cov)

    A = np.linalg.cholesky(cov)
    A = np.matrix(A)
    x = np.zeros((m,n))
    
    for i in range(0, m):
        x[i,:] = np.dot(A,np.random.normal(0,1, size=n))

    return x

def covarianceMatrix(n):#Generates an nxn random positive definite (covariance) matrix    
       
    A = np.random.random(size=(n,n))
    A = np.matmul(A,np.matrix.transpose(A))#multiply matrix with its transpose to get positve definite cov. mat.    
    return A
    
def twoPointAvg(j,k, cov, m):#Finds ensebmle avergae at xj and xk
    temp = gaussVecGen(cov, m) #if we want to use gaussVecGen()
    temp2 = temp[:,j] * temp[:,k]
    
    avg = np.average(temp2)
    
    return avg
    #print("Average of xj and xk is = ",  avg)
    #print("Matrix element is = ", cov[j][k])
    

#cov = covarianceMatrix(90)
#twoPointAvg(30,30,cov,10000)

def playingAround(m,n):
    
    points = np.zeros(n)
    cov = covarianceMatrix(n)
    
    for i in range(0,n):
        points[i] = twoPointAvg(i,i,cov,m)
    
    X = np.arange(0,n)
    #print(len(X))
    #print(len(points))    
    plt.scatter(X,points)
    plt.show()

playingAround(100,100)
"""
"""
a = dblquad(special.sph_harm(2,2,lambda theta: theta,lambda phi: phi), 0, np.pi, lambda x:0, lambda x: 2*np.pi)
print(a)
"""

#Y = lambda l, m, theta, phi: special.sph_harm(m, l, phi, theta)
#Y22 = lambda theta, phi: Y(2,2,theta,phi)
#Y22 = lambda theta, phi: special.sph_harm(2,2,phi, theta)
#print(Y(2,2,np.linspace(0,np.pi,100),np.linspace(0,2*np.pi,100)))
#AbsY22 = lambda theta, phi: np.abs(Y22(phi, theta))**2
"""
AbsY22 = lambda theta, phi: np.abs(special.sph_harm(2,2,phi, theta))**2
plt.plot(np.linspace(0,100,1000),AbsY22(np.linspace(0,100,1000),np.linspace(0,100,1000)))
"""
"""
def complex_quadrature(func, a, b, **kwargs):
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])



E = lambda x : np.exp(10*x*1j)
#a = quad(np.real,0,1)[0]
#b = quad(np.imag(E),0,1)
#c = a +b*1j
d = complex_quadrature(E, 0, 1)
#realIntegral = quad(np.real(E(lambda(x)),0,1)
#imagIntegral = quad(np.imag(E(x)),0,1)
#g = realIntegral[0] +1j*imagIntegral[0]
print(d)
"""
"""
def bessel(q):
    bessel, der = special.sph_jn(2, q)
    return bessel[2]

a =quad(bessel, 0,10)

print(a)
"""
"""
f = lambda x: np.exp(-x) * np.sin(x)
a = quad(f,0, float('inf'))
print(a)
"""
"""
def j2(q):
    bessel, der = special.sph_jn(2, q)
    return bessel[2]
#plt.plot(j2(np.arange(1,100,1)))
temp = quad(j2, 0,float('inf'))
print(temp)
"""
"""
result = quad(lambda x: special.jv(4,x), 0, 20)
print(result) 
print( "Gaussian integral", np.sqrt(np.pi),quad(lambda x: np.exp(-x**2),-np.inf, np.inf))
"""
#f = lambda q: -((3*q*np.cos(q) + (-3 + q**2)*np.sin(q))/q**3)
#print(f(100))
#temp = quad(j2, 0, float('inf'))
"""
n=10
x = np.arange(n).astype(float) * n/100

print(x)
"""
"""
print(D.delta(100,1.5))
"""
"""
a = Rotation.D(4,3,4,1,1.5,2)
b = a.doit()
print(sp.N(b))
"""
"""
a = spec.sph_harm(2,2,np.pi/9,np.pi/3)
print(a)

b = sp.mpmath.spherharm(2,2,np.pi/3,np.pi/9)
print(b)

sp.mpmath.spherharm(2,k,np.pi/6, np.pi/9)
"""
"""
c = 0

b =  lambda k: Rotation.D(2,k,-2, np.pi/6, np.pi/2,0).doit() * spec.sph_harm(k,2,np.pi/6, np.pi/12)
print(spec.sph_harm(-2,2,0,-5*np.pi/12))
#b =  lambda k: Rotation.D(2,k,2, np.pi/6, np.pi/9,0).doit() * spec.sph_harm(k,2,np.pi/6, np.pi/9)
for i in range(-2,3,1):
    c += b(i)
print(sp.N(c)._eval_evalf(10))
"""
"""
n = 50#number of random clusters we want

z_testList = np.random.uniform(0.0,3,n)
phi_testList = np.random.uniform(0.0, 2*np.pi,n)
theta_testList = np.arccos(np.random.uniform(-1,1,n))

rand50Clusters = np.vstack((z_testList, theta_testList, phi_testList)).T
np.save("rand50Clusters", rand50Clusters)
"""
"""
x = Symbol('x')
y = Symbol('y')
a = diff(sin(x),x)
a.doit()
print(a(2))
"""
"""
def P(k,f):
    if(k==0):
        return 0
        
    elif(k<2):
        return (1-f)/k**3
    
    else:
        return 1/k**3

f = Symbol('f')
#dpdf = lambda k: diff(P(k,f),f).doit()
#print(dpdf(0.4))
plt.figure(8)
plt.plot(np.arange(0,2,0.1),[P(k,0.7)for k in np.arange(0,2,.1)])
plt.show()
"""
"""
n=0
a = spec.sph_jn(n,0)[0][n]
"""

def pToyModel(k):
    if k==0:
        return 0
    elif(k<2):
        return (1-0.7)/k**3
    return 1/k**3
    
plt.plot([pToyModel(k) for k in range(0,20)])
plt.show()
