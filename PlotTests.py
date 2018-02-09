# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:00:23 2017

Routine for creating some common plots. The immediate need is to make reproduce some plots from Hall and Challinor amongst some other basic test cases.
For more details see pg 57 in notebook.  

This routine makes plots for the CROSS CORELATIONS. 

Elm(r) = <alm(0) a22(r)>/[(Cl^1/2) (<|a22(r)|>^1/2)]

@author: arsalanadil
"""

import numpy as np, matplotlib.pyplot as plt
import VarMatrix_Generalized as Remote
import LocalPS as Local
import VarMatrix_Gen_Local as Cross

#covMat = np.load("XTest1.npy")

def powerSpec(k):
    
    if(k==0):
        return 0
    
    else:
        return 1/k**3

def LocalCMB(l):
    local = Local.Integral(l,powerSpec)
    
    return local    
    
def RemoteCMB(r):
    location = np.array([r,0,0])
    temp = Remote.wignerRotation(location,location)[0]


    return temp[0][0]


def CrossCMB(r,l,m):
    location = np.array([r,0,0])
    
    crosscorel = Cross.wignerRotation(location, l, powerSpec)
    print(crosscorel)    
    
    return crosscorel[0][0][l-m]#we generally want El2 (i.e. m=2)
    #this is to find the right index. So, for e.g., for l=4, the matrix is arranged as, "m=-4,-3,-2,-1,..." and we want m=-2 (or 2)
    
def makePlot(r,lmax,m):
    plot = np.zeros(lmax-1)    
    remote = np.abs(RemoteCMB(r))
    
    i = 0
    for l in range(2,lmax+1):
        
        cross = CrossCMB(r,l,m)
        local = 256.2792 * np.abs(LocalCMB(l))#256.2792 is the normalization constant for the LocalCMB
      
        Elm = np.abs(cross/(np.sqrt(local)*np.sqrt(remote)) )
        print("Numerator:", cross, ",local:", local, "remote", remote)        
        #plot[i] = Elm
        plot[i] = Elm**2#squaring for Summed Correlation plot ONLY!
        print("Elm is", Elm)
        i+=1
    
    return plot
    
    
        
"""
z = 0
lmax = 4
m = 2
temp = makePlot(z,lmax,m)    
print(temp)
plt.plot(np.array([2,3,4]), temp)
#sumCor = np.zeros(10)
"""



"""
n = 10    
sumCor = np.zeros(n)
#zVar = np.zeros((n,3))
for z in range(0,n):
    sumCor[z] = np.sum(makePlot(z/10,lmax,m)) 
    #zVar[i] = makePlot(z/10,lmax,m)

print(sumCor)
plt.plot(np.arange(0,10)/10,sumCor)
"""
lmax = 5
m = 2

z = 3
summedLVals = np.zeros(lmax-1)
for l in range(2,lmax+1):
    summedLVals[l-2] = np.sum(makePlot(z,l,m))

print("The Sum over L values is:", summedLVals)
plt.figure(2)
plt.show()    
plt.plot(np.arange(2,lmax+1), summedLVals)
plt.title("Running Sum over L at different RedShifts")