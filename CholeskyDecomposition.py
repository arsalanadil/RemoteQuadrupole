# -*- coding: utf-8 -*-
"""
Created on Wed May 25 14:49:30 2016

This routine generates a random "real" relation matrix and uses that to generate a set of
n-dimensional, m-gaussian random vectors. The output matrix is hence mxn.
We use the Cholesky decomposition routine from numpy to acheive this. For more details, refer to
"Gaussian Random Vectors and a bunch of other stuff, Ted Bunn"

@author: arsalanadil
"""

import numpy as np, matplotlib.pyplot as plt

"""
Important note: Both multivarNormal() and gaussVecGen() accomplish the same purpose, namely to generate an
mxn matrix where each row is an N-dimensional vector "x". gaussVecGen() uses the Cholesky decomposition approach
while multiveNormal() uses a predefined routine in python to generate x form a multivariate distribution
"""

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

    A = np.zeros((n,n))
   
    #np.fill_diagonal(A, np.random.normal(0,1))    
    #np.fill_diagonal(A, np.random.randint(1,5))
        
    A = np.random.random(size=(n,n))
    A = np.matmul(A,np.matrix.transpose(A))#multiply matrix with its transpose to get positve definite cov. mat.    
    #print(A)        
    return A
    
def twoPointAvg(j,k, cov, m):
    temp = gaussVecGen(cov, m) #if we want to use gaussVecGen()
    #temp = multivarNormal(cov, m) #if we want to use multivarNormal() routine    
    #print(temp)
    temp2 = temp[:,j] * temp[:,k]
    
    avg = np.average(temp2)
    
    print("Average of xj and xk is = ",  avg)
    print("Matrix element is = ", cov[j][k])
    

#inputCov = np.array([[1 , 0] , [0, 1]])
cov = covarianceMatrix(9000)
#print(cov)
#x = gaussVecGen(cov, 10)
#twoPointAvg(1,89,cov,10000)
#print("Output matrix is: \n", x )                   
#print("Avg at one row = ", np.mean(x[500]))
#plt.hist(x[200],30)
#plt.show() 
plt.figure(2)
plt.imshow(gaussVecGen(cov,10000)) 
plt.show()