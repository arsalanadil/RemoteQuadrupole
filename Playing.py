# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 16:37:15 2016

This routine generates a Variance matrix based on a given power spectrum and takes its cholesky decomposition
to give an mxn matrix of random gaussian vectors.

@author: arsalanadil
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 25 14:49:30 2016

This routine generates a random "real" relation matrix and uses that to generate a set of
n-dimensional, m-gaussian random vectors. The output matrix is hence mxn.
We use the Cholesky decomposition routine from numpy to acheive this. For more details, refer to
"Gaussian Random Vectors and a bunch of other stuff, Ted Bunn"

@author: arsalanadil
"""

import numpy as np, matplotlib.pyplot as plt, healpy as hp
import scipy.special as special
from scipy.integrate import dblquad, quad
import TransferFunction as tf
"""
Important note: Both multivarNormal() and gaussVecGen() accomplish the same purpose, namely to generate an
mxn matrix where each row is an N-dimensional vector "x". gaussVecGen() uses the Cholesky decomposition approach
while multiveNormal() uses a predefined routine in python to generate x form a multivariate distribution
"""
print(tf.delta(5,2))

def gaussVecGen(cov, m):
    
    n = len(cov)

    A = np.linalg.cholesky(cov)
    A = np.matrix(A)
    #print("Cholesky decomposition matrix A: \n", A)
    x = np.zeros((m,n))
    
    """
    One simple way to make this faster is to pass in values for "j" and "k" and
    just compute those vectors as they're the only ones we're interested in
    but the drawback is that we've lost our "matrix", i.e. we haven't computed any 
    other vector
    """
    for i in range(0, m):
        x[i,:] = np.dot(A,np.random.normal(0,1, size=n))
        
    
    #print("Avg value at one row = ", np.average(x[3]))
    #plt.hist(x[:,3])
    #plt.show()        
    
    #print("The cov matrix is \n",  cov)
    
    #print(x)
    return x

def multivarNormal(cov, m):
    n = len(cov)
    x= np.zeros((m,n))
    
    
    for i in range(0,m):  
        x[i,:] = np.random.multivariate_normal(np.zeros(n), cov)

    
    #x =  np.random.multivariate_normal(np.zeros(n), cov, size=(m,n))#This is fine, but just gives transpose of what I want
    return x

def covarianceMatrix(n):#Generates an nxn random positive definite (covariance) matrix    
    
      # -*- coding: utf-8 -*-
    """
    Created on Tue Jun  7 16:29:38 2016
    
    This routine attempts to generate a variance matrix based on the average of two 
    quadrupole spherical harmonic coefficientc, {a2,m(r)}. This is outlined in eq.2 in Bunn (2006)
    
    @author: arsalanadil
    """
    
    import scipy.special as special
    import numpy as np, matplotlib.pyplot as plt
    from scipy.integrate import dblquad, quad
    
    
    """
    phi = np.linspace(0,2*np.pi,100)
    theta = np.linspace(0,np.pi,100)
    #print(theta)
    
    Y22 = special.sph_harm(2,2,theta,phi)
    AbsY22 = (np.abs(Y22))**2
    
    #print(AbsY22)
    
    #y = (15*((np.sin(theta))**4))/(32*np.pi)
    ax = plt.gca()
    plt.axis([0,100,0,0.16])
    #ax.set_autoscale_on(False)
    plt.plot(theta,AbsY22)
    """
    
    """
    The methods doubleIntegral() and integral() are slightly modified scipy routines that
    allow us compute integrals of complex functions. todo: write a class and generalise this
    to n-dimensions
    
    to return uncertainty in integral calculation, append ", real_integral[1], imag_integral[1] within return statement"
    
    
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
        real_integral = quad(real_func, a, b, limit = 500, **kwargs)
        imag_integral = quad(imag_func, a, b,limit = 500, **kwargs)
        return (real_integral[0] + 1j*imag_integral[0])
    
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
        
        j1 = lambda q: -((3*q*np.cos(q) + (-3 + q**2)*np.sin(q))/q**3)
        
        
        func = lambda k: ((k**2)*powerSpec(k)*f(k)*j1(k*q1)*j1(k*q2))
        
        temp = ((4*np.pi)**2)*integral(func, 0, np.inf)
        #temp = quad(j2, 0, float('inf'))
        
        return temp
    
    #rj = 0.1
    #rk = 0.2
    R=1
    p = lambda k: 1/k**2
    #f = lambda k: doubleIntegral(F(k,rj,rk),0,2*np.pi, 0, np.pi)    
    #GAMMAjk = Gamma(rj,rk,R,p)
    #print(GAMMAjk)
    #n = 10
    Variance = np.zeros((n,n))
    for x in range(n):
        for y in range(n):
            Variance[x][y] = np.real(Gamma(x/n,y/n,R,p))
    
    return Variance
    
def twoPointAvg(j,k, cov, m):
    temp = gaussVecGen(cov, m) #if we want to use gaussVecGen()
    #temp = multivarNormal(cov, m) #if we want to use multivarNormal() routine    
    #print(temp)
    temp2 = temp[:,j] * temp[:,k]
    
    avg = np.average(temp2)
    
    print("Average of xj and xk is = ",  avg)
    print("Matrix element is = ", cov[j][k])
    

#inputCov = np.array([[1 , 0] , [0, 1]])
cov = covarianceMatrix(500)
#plt.imshow(gaussVecGen(cov,1000))
plt.figure(3)
#plt.plot(cov[0,:])
#z = np.ndarray.flatten(gaussVecGen(cov,55))
#cs = plt.pcolor(z)
#cb = plt.colorbar(cs, orientation = 'horizontal')
#plt.xlim(0,np.shape(z))
#plt.ylim(0,np.shape(z))
#plt.show()

#print(cov)
x = gaussVecGen(cov, 200)
plt.plot(x[2])
plt.show()
#twoPointAvg(1,1,cov,1000)
#print("Output matrix is: \n", x )                   
#print("Avg at one row = ", np.mean(x[500]))
#plt.hist(x[200],30)
#plt.show()  

#map = hp.sphtfunc.synfast(z,128)
#hp.mollview(map)
