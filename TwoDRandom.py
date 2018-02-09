# -*- coding: utf-8 -*-
"""
Created on Wed May 11 14:05:27 2016

@author: arsalanadil
"""

import numpy as np, matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.graph_objs as go
from mpl_toolkits.mplot3d import Axes3D
from pylab import *


def P(k):
    if(k[0][0]==0 and k[0][1] == 0):
        return 0
    k = np.sqrt(k[0][0]**2 + k[0][1]**2)
    
    return 1/k**2
 
def P2(k):
    if(k[0][0]==0 and k[0][1] == 0):
        return 0
    k = np.sqrt(k[0][0]**2 + k[0][1]**2)
    temp =0.
    if (k>= 20 and k<21):
        temp = 1.
    return temp
 

   
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
    
def pointMatrix(P, npoints, nmaps, whichx, whichy):
    
    temp = np.zeros(nmaps)
    
    for x in range(0, nmaps):
        temp[x] = (npoints**2)*randomFunction(npoints,P)[whichx][whichy]      
    
  
    return temp
    
#print(pointMatrix(P, 16,1000, 1,0))
    
"""   
def plotSurface(P, npoints, nmaps, whichx, whichy):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    a=np.arange(0,npoints)
    x,y=np.meshgrid(a,a)#This and the previous line do the same thing as np.arrange in 2D. This is stupid.
    #print(x)    
    
    ax.plot_surface(x,y,pointMatrix(P,npoints, nmaps, whichx,whichy))
        
"""        
#plotSurface(P, 3, 4, 1, 2)
#plt.imshow(randomFunction(256,P))

#plt.imshow(pointMatrix(P,16,100,11,3))
   
#plt.hist( pointMatrix(P,100,500,12,12),30)
#plt.hist( pointMatrix(P,100,500,1,98),30)
#plt.hist( pointMatrix(P,100,500,50,70),30)
#plt.hist( pointMatrix(P,100,500,1,2),30)
#plt.hist( pointMatrix(P,100,500,97,98),30)

#plt.show()


"""
The following find the two-point corelation function in a given map
"""

def twoPointAvg(P, npoints, nmaps, point1x, point1y, point2x, point2y):
    
    temp2 = np.zeros(nmaps)    
    
    for x in range(0, nmaps):
        temp = (npoints**2)*randomFunction(npoints, P)
        temp2[x] = temp[point1x][point1y]*temp[point2x][point2y]
    
    average = np.mean(temp2)
    
    print("The average of point", point1x, "," , point1y, "and", point2x, "," , point2y, "is", average)
    
   
    return average  
    
twoPointAvg(P, 256, 500, 1,3, 1,3)
     
point = pointMatrix(P, 256,500, 1, 3)
print ("Mean is " , np.mean(point))
print ("Standard Deviation is " , np.std(point))
print ("Standard Deviation Squared is " , (np.std(point))**2)