# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 16:29:38 2016

This routine attempts to generate a variance matrix based on the average of two 
quadrupole spherical harmonic coefficientc, {a2m(r)}. This is outlined in eq.2 in Bunn (2006)

ToDo: Increase limit/mav no. of allowed subdivisions in quad/dblquad. (Then see if it converges at infinity)


@author: arsalanadil
"""

import scipy.special as special
import numpy as np, matplotlib.pyplot as plt
from scipy.integrate import dblquad, quad


"""
The methods doubleIntegral() and integral() are slightly modified scipy routines that
allow us to compute integrals of complex functions. todo: write a class and generalise this
to n-dimensions

to return uncertainty in integral calculation, append ", real_integral[1], imag_integral[1] within return statement"


"""
"""
def doubleIntegral(func, a, b, c, d, **kwargs):
    def real_func(k,x):
        return np.real(func(k,x))#real and imag parts need to be integrated seperately
    def imag_func(k,x):
        return np.imag(func(k,x))
    real_integral = dblquad(real_func, a, b, lambda y: c, lambda y: d, **kwargs)
    imag_integral = dblquad(imag_func, a, b, lambda y: c, lambda y: d, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0])

def integral(func, a, b, **kwargs):
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0])
"""


"""
this funciton needs to be calculated only once. We've calculated it and verified the general 
result using Mathematica. The output is labelled "f" in the funtion Gamma(..) that follows below, i.e.,

f = lambda k: (-15*(3*k*q*np.cos(k*q) + (-3 + (k*q)**2)*np.sin(k*q)))/(k*q)**5

where k is our running variable and q = rj -rk 
"""
def F(k,rj, rk):
    if(rj == rk):
        return 1
        
    E = lambda theta: np.exp(k*1j*(rj -rk)*np.cos(theta))
    
    AbsY22 = lambda theta, phi: np.abs(special.sph_harm(2,2,phi, theta))**2
        
    
    sinTheta = lambda theta: np.sin(theta)
    #return (E*AbsY22(lambda(theta,phi))*sinTheta(10))
    temp = lambda theta, phi: E(theta)*AbsY22(theta, phi)*sinTheta(theta)      
    return temp

"""
def j2(q):#There's nothing wrong with this, except Python has a stack overflow when integrating this to infinity
    bessel, der = special.sph_jn(2, q)
    return bessel[2]
"""
      

def Gamma(rj, rk, R, powerSpec):
    #f = lambda k: doubleIntegral(F(k,rj,rk),0,2*np.pi, 0, np.pi)#if we want to invoke our routin F(...) above     
    q = rj-rk
    q1 = R-rj
    q2 = R-rk
   
    if(rj == rk):
        f = lambda k: k**0
    else:

        f = lambda k: (-15*(3*k*q*np.cos(k*q) + (-3 + (k*q)**2)*np.sin(k*q)))/(k*q)**5 
    
    j1 = lambda q: -(-4*np.pi/3)*((3*q*np.cos(q) + (-3 + q**2)*np.sin(q))/q**3)
    
    
    func = lambda k: ((k**2)*powerSpec(k)*f(k)*j1(k*q1)*j1(k*q2))
    #print(f(0.000001))    
    #plt.figure(4)    
    #plt.plot([(j1(k*q1)*j1(k*q2)) for k in range(0,50)])
    #plt.plot([func(k) for k in range(10)])
   
   
   #temp = ((4*np.pi)**2)*quad(func, 0, np.inf, limit =2000)[0]
    temp = quad(func, 0, 20, limit =30)[0]
    #temp = quad(j2, 0, float('inf'))
    
    return temp

rj = 1.45585
rk = 0.76248
R= 3.11293
def p(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**3
#f = lambda k: doubleIntegral(F(k,rj,rk),0,2*np.pi, 0, np.pi)    
GAMMAjk = Gamma(rj,rk,R,p)
print(GAMMAjk)

"""
n = 5
Variance = np.zeros((n,n))
for x in range(n):
    for y in range(n):
        Variance[x][y] = np.real(Gamma(x/n,y/n,R,p))

#print(Variance)
"""
#A = np.linalg.cholesky(Variance)
#print(A) 
#print(Gamma(1.456,0.762,3.11,p))