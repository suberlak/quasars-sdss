# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:39:28 2015

@author: suberlak
"""

import numpy as np
import matplotlib.pyplot as plt 


allqso = np.genfromtxt('QSO_4180560_lines.txt', dtype=str)

delflx = allqso[:,0].astype(float)
tau = allqso[:,1].astype(float)
names = allqso[:,3]

uniq_names = np.unique(names)
color_list = np.zeros_like(tau)
color_idx = np.linspace(0, 1,len(uniq_names))

for i in range(len(uniq_names)):
    mask = np.in1d(names, uniq_names[i])
    color_list[mask] = color_idx[i]
    print 'comparing ', i, '-th obj. id'
    
plt.scatter(np.log10(tau[:1e5]), delflx[:1e5], c=color_list[:1e5])

''' 
Tylko ze jak robie ten wykres to sie robi absolutny miszmasz, i nie wiadomo 
co mam z tego wyniesc... 

'''