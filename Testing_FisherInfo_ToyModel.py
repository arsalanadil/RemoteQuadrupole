# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:04:18 2016

This routine generates a random vector (x1....xn,y1...yn) based on the "real valued" covariance matrix
for a given value of 'f' (the suppression factor in our toy model). From our fisher info, we have an expected error
(lower bound on standard devation) on f. 
Once we have all this we want to play the reverse game, i.e., estimate f given our random vector, d. We use a log likehood
function to do so (see notebook for mathematical details). When we do this very many times, we want the mean of our ESTIMATED
f to converge to the true value of f and the standard deviation to converge to the Cramer-Rao bound, i.e. F**-.5 

@author: arsalanadil
"""

import numpy as np, matplotlib.pyplot as plt
import FisherInfo_toy_mdoel as fi
#from FisherInfo_toy_mdoel import Score
from sympy import *
import scipy.optimize

#covMat = fi.GammaR()
#n = len(covMat)
"""
VarAndCor10Clusters = np.load("VarAndCorMats_10_Clusters.npz")

varMatrix = VarAndCor10Clusters["varMatrix"]
corMatrix = VarAndCor10Clusters["corMatrix"]
def Gamma0():
      

    GammaXX = 1/2 * (np.real(varMatrix) + np.real(corMatrix))
    GammaXY = 1/2 * ((np.imag(corMatrix) - np.imag(varMatrix)))
    GammaYX = 1/2 * ((np.imag(corMatrix) + np.imag(varMatrix)))
    GammaYY = -1/2 *((np.real(corMatrix) - np.real(varMatrix)))
        
    gamma = np.bmat([[GammaXX , GammaXY],[GammaYX, GammaYY]])
    
    return gamma
"""

GammaRNH =  np.load("GammaReal_null_hyp_20Clusters.npy")#unsuppressed covariance matrix; i.e. the null hypothesis
GammaScore = np.load("GammaScore.npy")#GammaPrime for suppressed powerSpec(f=0.6)
GammaRTM = np.load("GammaReal_ToyModel_20Clusters.npy")#suppressed power spec(f=0.6)

    
def vecGenerator():
    A = np.linalg.cholesky(GammaRTM)
    d = np.dot(A, np.random.normal(0,1,size=len(GammaRTM)))
    #either of these methods should be fine
    #d = np.random.multivariate_normal(np.zeros(40),GammaRNH)
    return d


def gammaF(f):
    
    Gamma = GammaRNH + f*GammaScore
    #Gamma = GammaRNH
    return Gamma

def logLH(f,d):
    Lambda = np.linalg.multi_dot([d,np.linalg.inv(gammaF(f)),np.transpose(d)]) + np.log(np.linalg.det(gammaF(f)))
    #print('*****',f,Lambda[0,0])
    #return Lambda[0,0]
    #return Lambda
    #print('********',f, '    ' , Lambda)
    return  Lambda
"""
def Lambda(x):
    return logLH(x[0])
"""    

runs = 100000
manyRuns = np.zeros(runs)
for i in range(runs):
    d = vecGenerator()
    temp = lambda f: logLH(f,d)
    #print(scipy.optimize.minimize_scalar(logLH)['x'])    
    manyRuns[i] = scipy.optimize.minimize_scalar(temp, bounds = (0,1.3), method = 'bounded')['x']
    print(i)
print("Average:", np.average(manyRuns), "\n Std:", np.std(manyRuns))
"""
flist = np.arange(-1,2,0.1)
d = vecGenerator()
tempFun = lambda f: logLH(f,d)
temp = [tempFun(f) for f in flist]
#plt.figure()   
plt.plot(flist,temp)
plt.show()
#print(scipy.optimize.minimize(tempFun, 0.02))
print(scipy.optimize.minimize_scalar(tempFun, bounds = (0,1.2), method = 'bounded'))
#print(np.argmin(temp)+2)
"""
"""
fMax = 10
manyRuns = np.zeros(100)
temp = np.zeros(fMax)
for i in range(len(manyRuns)):
    for f in range(len(temp)):
        temp[f] = logLH(f/10)
    fVal = np.argmin(temp)
    manyRuns[i] = fVal/10
    print(i)
print("Average:", np.average(manyRuns), "\n Std:", np.std(manyRuns))
"""