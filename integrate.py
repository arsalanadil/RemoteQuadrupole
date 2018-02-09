# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 12:53:27 2016

The methods doubleIntegral() and integral() are slightly modified scipy routines that
allow us to compute integrals of complex functions. todo: write a class and generalise this
to n-dimensions

to return uncertainty in integral calculation, turn UNCERTAINTY to 'ON' 


@author: arsalanadil
"""
import numpy as np
from scipy.integrate import dblquad,quad

UNCERTAINTY = 'OFF'

def doubleIntegral(func, a, b, c, d, **kwargs):
    def real_func(x,y):
        return np.real(func(x,y))#real and imag parts need to be integrated seperately
    def imag_func(x,y):
        return np.imag(func(x,y))
    
    real_integral = dblquad(real_func, a, b, lambda t: c, lambda t: d, **kwargs)
    imag_integral = dblquad(imag_func, a, b, lambda t: c, lambda t: d, **kwargs)
    
    val = real_integral[0] + 1j*imag_integral[0]
    
    error = (real_integral[1] + 1j*real_integral[1])
    if(UNCERTAINTY == 'ON'):
        return(val, error)
    return val
    
def integral(func, a, b, **kwargs):
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)

    val = real_integral[0] + 1j*imag_integral[0]
    
    error = (real_integral[1] + 1j*real_integral[1])
    if(UNCERTAINTY == 'ON'):
        return(val, error)
    return val