# -*- coding: utf-8 -*-
"""
Created on Tue May 10 16:36:36 2016

@author: arsalanadil
"""

import matplotlib.pyplot as plt, numpy as np
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = Axes3D(fig)
t = np.linspace(0, 5*np.pi, 501)
ax.plot(np.cos(t), np.sin(t),t)