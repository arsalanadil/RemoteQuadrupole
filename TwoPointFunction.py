# -*- coding: utf-8 -*-
"""
Simulate a two point function; see onepointfunction for more details
@author: arsalanadil
"""

import numpy as np


def P(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**2


def P1(k):
    return np.exp(-(k/5)**2)

def P2(k):
    return np.exp(-(k/10)**2)
        
def kVal(j , npoints):

    if(j <= npoints/2):
        return (j - 1)
    
    return npoints + 1 -j
    
def randomFunction(npoints, P):
  
  ft = np.zeros(npoints).astype(complex)#Need to specify that array type is complex
    
  for x in range(0, npoints):
            ft[x] = (np.random.normal(0,1) + 1j*np.random.normal(0,1))
            ft[x] = ft[x]*np.sqrt(P(kVal(x,npoints) )) 
        
  complexFFT = np.fft.ifft(ft,axis = -1)#Take the inverse FFT         
        #Note: The default normalization factor is 1/N 
        
  realFFT = np.real(complexFFT)#Drop the complex parts in each element
                
        
  return realFFT
  
def twoPointAvg(P, npoints, nmaps, point1, point2):
    
    temp2 = np.zeros(nmaps)    
    
    for x in range(0, nmaps):
        temp = npoints*randomFunction(npoints, P)
        temp2[x] = temp[point1]*temp[point2]
    
    average = np.mean(temp2)
    
    print("The average of point", point1, "and", point2, "is", average)
    

    
    return average

#twoPointAvg(P, 256, 1000, 30,30)
#twoPointAvg(P, 256, 1000, 1,100)
#twoPointAvg(P, 256, 1000, 1,10)


twoPointAvg(P, 800, 2000, 30,30)
twoPointAvg(P, 1024, 2000, 1, 1023)


#twoPointAvg(P2, 256, 1000, 30,30)
#twoPointAvg(P, 800, 2000, 30,30)
#twoPointAvg(P2, 800, 2000, 1, 799)


#twoPointAvg(P, 800, 2000, 10,20)
#twoPointAvg(P, 800, 2000, 240,250)
    
#print("The average of two points is", twoPointAvg(P, 256, 1000, 1,2))
#print ("Variance (Standard Deviation Squared) is " , ((np.std(pointList(P, 256,1000, 30)))**2))
