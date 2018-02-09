# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:15:48 2016

@author: arsalanadil
"""
import numpy as np, matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

"""
Note: if I have 1 only in the first position (zeros everywhere else) of 16x16 matrix 
then the ift gives a matrix where each element = 1/16^2. Hence the normalization by N^2
"""
def plotFunction(npoints, whichx, whichy):
    ft = np.zeros((npoints,npoints))
    ft[whichx][whichy] = 1

    #print("printing ft:")   
    #print (ft)



    ift = (npoints**2)*np.fft.ifft2(ft).astype(complex)

    #print("Now prinitng inverse FFT:")
    #print (ift)
    
    return ift


def plotSurface(npoints, whichx, whichy):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    a=np.arange(0,npoints)
    x,y=np.meshgrid(a,a)#This and the previous line do the same thing as np.arrange in 2D. This is stupid.
    ax.plot_surface(x,y,np.real(plotFunction(npoints,whichx,whichy)),label = "Real part")
    ax.plot_surface(x,y,np.imag(plotFunction(npoints,whichx,whichy)),color = 'r', label = "Complex part", )
    
def plotContour(npoints, whichx, whichy, fignum):
    figure(fignum)
    plt.imshow(np.imag(plotFunction(npoints,whichx, whichy)))
    plt.title("Imaginary Part")
    #plt.colorbar()
    
    
    figure(fignum+1)
    plt.imshow(np.real(plotFunction(npoints,whichx, whichy)))    
    #plt.colorbar()
    plt.title("Real Part")


plotContour(128,5,0,5)

