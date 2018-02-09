# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 15:28:05 2016

Generates the Variance and Relation matrices by extracting the required
elements from the matrix produced by VarMatrix_Generalized.py

Remember: Variance matrix must be Hermitian and Positive Definite
Relation matrix must be symmetric and 0 on the diagonals

If we're going to do stuff with Fisher Information in here, it is probably
worthwhile to (for testing purposes) generate and store
these matrices in a seperate file 

@author: arsalanadil
"""

import numpy as np
#import matplotlib.pyplot as plt
import VarMatrix_Generalized as vm
import time
#import FisherInfo_toy_mdoel as fi

randClusters = np.load("rand50Clusters.npy")

def pToyScore(k):
    k0 = 1.5    
    if(k==0):
        return 0
        
    elif(k>k0):
        return 0
    
    else:
        return -1/k**3

def pToyModel(k):
    if k==0:
        return 0
    elif(k<1.5):
        return (1-0.6)/k**3
        
    return 1/k**3
    
def powerSpec(k):
    if(k==0):
        return 0
    else:
        return 1/k**3
        
n = 2#number of galaxy cluster we want to consider
varMatrix = np.zeros((n,n)).astype(complex)
corMatrix = np.zeros((n,n)).astype(complex)
corMatrix2 = np.zeros((n,n)).astype(complex)#this should be the complex conjugate of the above corMatrix
t0 = time.time()
for i in range(n):
    for j in range(i+1):
        temp = vm.wignerRotation(randClusters[i,:],randClusters[j,:], powerSpec)[0]
        print("wignerRotation() returned:", temp[0])        
        varMatrix[i][j] = temp[0][0]
        varMatrix[j][i] = np.conjugate(temp[0][0])
        corMatrix[i][j] = temp[0][4]
        corMatrix[j][i] = temp[0][4]
    print(i)
        #corMatrix2[i][j] = temp[4][0]
t1 = time.time()
print("\n \n \n \n")    
print("VarianceMatrix: \n", varMatrix)
print("CorelationMatrix: \n", corMatrix)

#np.savez("VarAndCorMats", varMatrix = varMatrix, corMatrix=corMatrix)
"""
truncVar = np.zeros((n,n)).astype(complex)
truncCor = np.zeros((n,n)).astype(complex)

for x in range(n):
    for y in range(n):
        truncVar[x][y] = round(varMatrix[x][y],5)
        truncCor[x][y] = round(corMatrix[x][y],5)

print("Variance Matrix is hermitian? \n", (np.transpose(truncVar)==np.conjugate(truncVar)).all() )
print("Relation Matrix is symmetric? \n", (np.transpose(truncCor)==truncCor).all())
print("Time taken in generating all this is: ", (t1-t0))
"""




"""
We can now combine our variance and realtion matrices to form the 
block matrix, i.e. the (real valued) COVARIANCE matrix. 
"""
    
GammaXX = 1/2 * (np.real(varMatrix) + np.real(corMatrix))
GammaXY = 1/2 * ((np.imag(corMatrix) - np.imag(varMatrix)))
GammaYX = 1/2 * ((np.imag(corMatrix) + np.imag(varMatrix)))
GammaYY = -1/2 *((np.real(corMatrix) - np.real(varMatrix)))
        
gammaR = np.bmat([[GammaXX , GammaXY],[GammaYX, GammaYY]])
    

np.save("GammaReal_null_hyp_50Clusters", gammaR )