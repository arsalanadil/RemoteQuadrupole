# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:47:30 2016

This code computes the Fisher information for the toy model, i.e. the one where we 
exaplain the low quadrupole using a supression factor,f, in the power spectrum

@author: arsalanadil
"""
import numpy as np


#VarAndCor10Clusters = np.load("VarAndCorMats_10_Clusters.npz")
#VarAndCor10Clusters = np.load("VarAndCorMats_10_Clusters.npz")
#VarAndCor50Clusters = np.load("VarAndCorMats_50_Clusters.npz")
#GammaDerivMats = np.load("dgamma_term_fi.npz")


def powerSpec(k):
    k0 = 2    
    if(k==0):
        return 0
        
    elif(k>k0):
        return 0
    
    else:
        return -1/k**3
"""
def GammaR():
    
    varMatrix = VarAndCor10Clusters["varMatrix"]
    corMatrix = VarAndCor10Clusters["corMatrix"]

    GammaXX = 1/2 * (np.real(varMatrix) + np.real(corMatrix))
    GammaXY = 1/2 * ((np.imag(corMatrix) - np.imag(varMatrix)))
    GammaYX = 1/2 * ((np.imag(corMatrix) + np.imag(varMatrix)))
    GammaYY = -1/2 *((np.real(corMatrix) - np.real(varMatrix)))
    
    #GammaR = np.array(([GammaXX, GammaXY], [GammaYX, GammaYY])).reshape((20,20))
    
    GammaR = np.bmat([[GammaXX , GammaXY],[GammaYX, GammaYY]])
    
    return GammaR
    #Gamma = GammaXX + GammaYY + 1j*(GammaYX - GammaXY)
    #Relation = GammaXX - GammaYY + 1j*(GammaYX + GammaXY)
"""    
"""
def GammaR():
      
    varMatrix = VarAndCor10Clusters["varMatrix"]
    corMatrix = VarAndCor10Clusters["corMatrix"]
    
    GammaXX = 1/2 * (np.real(varMatrix) + np.real(corMatrix))
    GammaXY = 1/2 * ((np.imag(corMatrix) - np.imag(varMatrix)))
    GammaYX = 1/2 * ((np.imag(corMatrix) + np.imag(varMatrix)))
    GammaYY = -1/2 *((np.real(corMatrix) - np.real(varMatrix)))
        
    gamma = np.bmat([[GammaXX , GammaXY],[GammaYX, GammaYY]])
    
    return gamma

def Score():  
    dvarMatrix = GammaDerivMats["dvarMat"]
    dcorMatrix = GammaDerivMats["drelatMat"]
    
    
    dGammaXX = 1/2 * (np.real(dvarMatrix) + np.real(dcorMatrix))
    dGammaXY = 1/2 * ((np.imag(dcorMatrix) - np.imag(dvarMatrix)))
    dGammaYX = 1/2 * ((np.imag(dcorMatrix) + np.imag(dvarMatrix)))
    dGammaYY = -1/2 *((np.real(dcorMatrix) - np.real(dvarMatrix)))
    
    #dGamma = np.array(([dGammaXX, dGammaXY], [dGammaYX, dGammaYY])).reshape((20,20))
    dGamma = np.bmat([[dGammaXX, dGammaXY], [dGammaYX, dGammaYY]])
    
    return dGamma
"""
   
GammaScore = np.load("GammaScore.npy")#GammaPrime for suppressed powerSpec(f=0.6)
GammaRTM = np.load("GammaReal_ToyModel_20Clusters.npy")#suppressed power spec(f=0.6)

def FishyInfo():
        GammaRInv = np.linalg.inv(GammaRTM)
        dGamma = GammaScore        
        
        fishyInform = 1/2 * np.trace(np.linalg.multi_dot([GammaRInv,dGamma,GammaRInv,dGamma]))
        #fishyInform = 1/2 * np.trace(GammaRInv*dGamma*GammaRInv*dGamma)
        
        expError = (fishyInform)**(-.5)#what are the units of these numbers??
        
        return fishyInform, expError
        
print(FishyInfo())