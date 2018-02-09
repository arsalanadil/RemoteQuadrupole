# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:21:34 2016

This is meant to play around, in an effort to understand, the spherical harmonic transform routines
in healpy

@author: arsalanadil
"""

import numpy as np, matplotlib as plt, healpy as hp
np.set_printoptions(threshold = np.nan)#this just ensures we don't get a truncated array when it is printed
#This is just the tutorial. Everyhting works fine
import sympy as sp

"""
NSIDE = 32

m = np.arange(hp.nside2npix(NSIDE))
hp.mollview(m, title="Mollview image RING")
"""
#a = np.zeros(10)

#a[9] = 1

#for x in range(0,256):
#    a[x] = np.random.normal(0,1)

#print(a)

#map = hp.sphtfunc.synfast(a,32)

#print(map)
#hp.orthview(map)


#alm = hp.map2alm(map)
#Cl = hp.alm2cl(alm)

#print(Cl)
"""
Cl = np.zeros(8).astype(complex)
for k in range(1,8):
    Cl[k] = 1/k**3 + 1j*1/k**3

mapp = hp.sphtfunc.synfast(Cl,8)
albar = hp.sphtfunc.map2alm(mapp)
 
al = lambda m: sp.physics.quantum.spin.Rotation.D(2,m,2,np.pi/2,np.pi/2,0).doit()*albar
d = 0
for m in range(-2,3,1):
    d+= al(m)    
    
print(d)
#hp.mollview()
#print(hp.sphtfunc.alm2cl(a))
"""