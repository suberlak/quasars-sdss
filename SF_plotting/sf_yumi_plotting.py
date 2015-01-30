# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:24:02 2015

@author: astronomy

A program to read in the output of Yumi's calculations, and plot the same as 
sf_plotting.py  does 

"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import binned_statistic
import sys
import os

inDir = './sf/'
outDir = './sf/' 
if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 

inFiles = os.listdir(inDir)


# sf - those that start with sf ,         tau -  those that start with tau...


# do stats and plot !   

mask_sf = np.zeros_like(inFiles, dtype=bool)
mask_tau = np.zeros_like(inFiles, dtype=bool)

for i in range(len(inFiles)) : 
    if inFiles[i][:5] == 'sfout' : mask_sf[i] = True
    if inFiles[i][:6] == 'tauout' : mask_tau[i] = True
        

sf_files = inFiles[mask_sf]
tau_files = inFiles[mask_tau]

