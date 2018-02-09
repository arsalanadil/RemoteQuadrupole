# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 18:38:44 2016

Somebody wrote a code for calculating C-G coefficients that I found online. I'm testing it out and comparing
it with Mathematica.

UPDATE: The code works wonderfully(and I must admit, better than I had expected). I have made changes to it
such that it corresponds to the our convention of defining the first Ylm term to be conjugated.
Second, I have added a method called J() which is the expansion of the triple integral over three Ylm
coefficients, as defined by Zare in eq 3.115 (pg 102).

Moreover, this code is pretty optimized. One thing we have done is that for values of L=0,2,4 we have generated 
a list of values with the correspong M, m1, and m2 terms.

It then makes sense for me to change the name of this routine from CG_useless.py to ClebGord.py

@author: arsalanadil
"""

from scipy.special import binom
import numpy as np


def clebgor(j1, j2, j, m1, m2, m):
    """
    Parameters:    j1, j2, j: the angular momenta input
                   m1, m2, m: the z components of angular momenta input
    Returns:       The numerical value of the Clebsch-Gordan coeffcient.
    Remarks:       Note that in the sum none of the binomial coeffcients
                   can have negative values.  Thus, zmin is there to make
                   sure that the sums have a cut-off.
    """
    zmin = int(min([j1 - m1, j2 + m2]))
    J = j1 + j2 + j
    return (int(m1 + m2 == m) *
            int(np.abs(j1 - j2) <= j <= (j1 + j2)) *
            int(np.abs(m1) <= j1) *
            int(np.abs(m2) <= j2) *
            int(np.abs(m) <= j) *
            int((j1 + m1) >= 0.0) *
            int((j2 + m2) >= 0.0) *
            int((j + m) >= 0.0) *
            int(J >= 0) *
            np.sqrt(binom(2 * j1, J - 2 * j) *
                    binom(2 * j2, J - 2 * j) /
                    (binom(J + 1, J - 2 * j) *
                    binom(2 * j1, j1 - m1) *
                    binom(2 * j2, j2 - m2) *
                    binom(2 * j, j - m))) *
            np.sum([(-1) ** z * binom(J - 2 * j, z) *
                    binom(J - 2 * j2, j1 - m1 - z) *
                    binom(J - 2 * j1, j2 + m2 - z) for z in range(zmin + 1)]))
                    
def J(j1,j2,j,m1,m2,m):
    if((j1 + j2 + j)%2 != 0):
        return 0
    
    temp = j1#Switch {j1,m1} with {j,m} to follow convention
    temp1 = m1
    j1 = j
    m1 = m
    j = temp
    m = temp1
    
    """
    Now {j1,m1} are the parameters of the first Y2m^* term and {j,m} are conventionally
    defined to correspond to the YLM term (the so called "total angular momentum")
    """    
    
    cg = clebgor(j1, j2, j, m1, m2, m)
    cg0 = clebgor(j1,j2,j,0,0,0)
    val = (np.sqrt(((1 + 2*j1)*(1 + 2*j2))/((1 + 2*j)*4*np.pi)))*cg*cg0
    return val

"""
print("L M  m1  m2   J(2,2,L,m1,m2,M)")    
#array = np.zeros(20)
i = 0
for L in range(0,5,2):
    M = 0
    while(M<=L):
        for m1 in range (0,3):
            for m2 in range(0,3):
                #if(m1-m2 == M):
                    print(L,M, " ",m1, " ", m2, " ",J(2,2,L,m1,m2,M))
                    #array[i] = J(2,2,L,m1,m2,M)
                    i+=1
        M+=1
#print(array)
print("Count =", i)
"""

"""
print("L M   m1   m2     J(2,2,L,m1,m2,M)")    
#array = np.zeros(20)
i = 0
for L in range(0,5,2):
    M = -L
    while(M<=L):
        for m1 in range (-2,3):
            for m2 in range(-2,3):
                if(m1-m2 == M):
                    j = J(2,2,L,m1,m2,M)
                    print(L,M, "  ",m1, "  ", m2, "    ",j)
                    #print(L,M,m1,m2,J(2,2,L,m1,m2,M))
                    #array[i] = J(2,2,L,m1,m2,M)
                    i+=1
        M+=1  
print("Count = ", i)   
"""       