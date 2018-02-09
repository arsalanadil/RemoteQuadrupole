# -*- coding: utf-8 -*-
"""
Created on Thu May 12 18:24:41 2016

This program estimates the power spectrum of a randomly generated map.

Note that it first generates a random map using a given power spectrum and then (independently) attempts to estimate the original 
power spectrum given those points.

@author: arsalanadil
"""

import numpy as np, matplotlib.pyplot as plt
#import scipy
#from scipy import signal

def P(k):
    #if(k[0][0]==0 and k[0][1] == 0):
        #return 0
    k = np.sqrt(k[0][0]**2 + k[0][1]**2)
    
    return np.exp(-(k/5)**2)
 
    
def kVal(i, j, npoints):
    
    k = np.zeros((1,2))    
    
    if(i<= npoints/2):
        i = i - 1
        
    else:
        i = npoints + 1 - i
        
    if(j<= npoints/2):
        j = j - 1
        
    else:
        j = npoints + 1 - j
        
    k[0][0] = i
    k[0][1] = j
    
    return k

def randomFunction(npoints, P):
    ft = np.zeros((npoints, npoints)).astype(complex)
    #Let's only deal with square patches for now
    
    for x in range(0, npoints):
        for y in range(0, npoints):
            ft[x][y] = (np.random.normal(0,1) + 1j*np.random.normal(0,1))
            ft[x][y] = ft[x][y] * np.sqrt(P(kVal(x,y,npoints)))
            
    complexIFFT = np.fft.ifft2(ft)
    
    realIFFT = np.real(complexIFFT)
    
    #print(realIFFT)    
    
    return realIFFT
"""    
def fourierTransform(npoints, P):
    
   IFFT = randomFunction(npoints, P)#this gets us the randomly generated map
    
   FFT = np.fft.fft2(IFFT)
    
   return FFT
"""
"""    
#print(fourierTransform(128, P))
    
f, Pspec = signal.periodogram(fourierTransform(256,P))
plt.figure()
plt.plot(f, np.sqrt(Pspec))
plt.show()
"""

#print(np.std(fourierTransform(128,P)))

IFFT = randomFunction(256, P)#this gets us the randomly generated map
    
FFT = np.fft.fft2(IFFT)#And this gets us the FFT of the points from that map

def powerSpec(k, FFT):
    #kx = k[0][0]
    #ky= k[0][1]
    npoints = len(FFT)
    
    #num = np.zeros(npoints)
    #denom = np.zeros(npoints)    
    num = 0
    denom = 0
    P=0
    for  i in range(0, npoints):
        for j in range(0, npoints):
            kx = kVal(i,j,npoints)[0][0]
            ky = kVal(i,j,npoints)[0][1]           
            a = np.sqrt(kx**2 + ky**2)
            
            if(k-.5<a<k+.5):              
                num = num + (np.abs([FFT[i][j]]))**2
                denom = denom+1
    
    if(denom != 0):            
       P = num/denom
    
    
    return P
       

def powerSpecPlot(FFT, samplePts):
    
    pVals = np.zeros(samplePts)
    
    for i in range(0, samplePts):
        #j = i+1        
        pVals[i] = powerSpec(i,FFT)

    plt.figure()
    
    
    x = np.arange(1, 250, 1)
    y = np.exp(-(x/5)**2)#this is the actual power spectrum that we've used to generate the map. CHange this for testing
    plt.plot(x,y)
    plt.title("exp(-k/5)^2")
    
    plt.figure(2)    
    plt.scatter(np.arange(0,len(pVals)),pVals, c= 'red',)#this is the power spectrum we PREDICT
    
    ax = plt.gca()
    ax.set_xscale('log')#for exponentials, set x0axis on log scale. for 1/x^n, set y-axis to log
    ax.set_yscale('log')    
    plt.show()

    
powerSpecPlot(FFT, 256)
    