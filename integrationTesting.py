# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 13:00:21 2016

@author: arsalanadil
"""

import Integration as integrate

f = lambda x: x**2

a = integrate.integral(f,0,10)
print(a)