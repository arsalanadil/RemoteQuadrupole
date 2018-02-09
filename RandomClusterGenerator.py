# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 17:16:12 2016
This is just a simple code for generating random galaxy cluster points on the sky.
I just wanted it in a seperate place than everything else; hence an independent class.
It generate n cluster locations, i.e. the redshift  (z_list),
colatidunal angle (theta_list), and polar angle(phi_list);

and saves the resulting nx3 array as a .py file

@author: arsalanadil
"""
import numpy as np

n = 20#number of random clusters we want

zlist = np.random.uniform(0.1,3,n)#sort this list, thats the one we want
z_testList = np.sort(zlist)
phi_testList = np.random.uniform(0, np.pi/2,n)
theta_testList = np.arccos(np.random.uniform(0,1,n))

randClusters = np.vstack((z_testList, theta_testList, phi_testList)).T
#np.save("rand100Clusters", randClusters)
np.save("testClusters", randClusters)