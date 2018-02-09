# -*- coding: utf-8 -*-
"""
This is an attempt to translate Ted's Mathematica code, onepointfunction.nb, to python

@author: arsalanadil 
    
"""

import numpy as np, matplotlib.pyplot as plt

#User defined P(k_) funciton
def P(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**2
        
def kVal(j , npoints):

    if(j <= npoints/2):
        return (j - 1)
    
    return npoints + 1 -j
    
#Creates a "map" with npoints. The points are a random Gaussian distribution    
def randomFunction(npoints, P):
        
        ft = np.zeros(npoints).astype(complex)#Need to specify that array type is complex
       
        """
        for i in np.nditer(ft, op_flags=['readwrite']):
            i[...] = (np.random.normal(0,1) + 1j*np.random.normal(0,1))
        
        for j in np.nditer(ft, op_flags=['readwrite']):
            j[...] = j * np.sqrt(P(kVal(j,npoints) ))     
        
        """
        for x in range(0, npoints):
            ft[x] = (np.random.normal(0,1) + 1j*np.random.normal(0,1))
            ft[x] = ft[x]*np.sqrt(P(kVal(x,npoints) )) 
        
        complexFFT = np.fft.ifft(ft,axis = -1)#Take the inverse FFT         
        #Note: The default normalization factor is 1/N 
        
        realFFT = np.real(complexFFT)#Drop the complex parts in each element
                
        
        return realFFT

def pointList(P, npoints, nmaps, whichpoint):
       
    temp = np.zeros(nmaps)
    
    #Generates nmaps amount of maps (each containing npoints) and stores the value at "whichpoint" in an array
    for x in range(0, nmaps):
        temp[x] = npoints*randomFunction(npoints,P)[whichpoint]         
                    
    
    return temp    
    

#print (pointList(P,256,1000,30))
#plt.hist( pointList(P,250,1000,30),50)#The second parameter in plt.hist is the number of bins in the histogram

point =  pointList(P,256,1000,30)
print ("Mean is " , np.mean(point))
print ("Standard Deviation is " , np.std(point))
print ("Standard Deviation Squared is " , (np.std(point))**2)

