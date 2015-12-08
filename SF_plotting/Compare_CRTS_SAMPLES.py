# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 13:29:25 2015

@author: suberlak
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('Sample_4_Stars_blue_314248_lines.txt')

xi = data[:,0]
ei = data[:,1]
n = data[:,2]
rmag = data[:,3]

data1 = np.genfromtxt('Sample_starB_tau_0-1.7_204242_lines_1270_objects.txt')

xi1 = data1[:,0]
ei1= data1[:,1]
tau = data1[:,2]
n1 = data1[:,3]
rmag1 = data1[:,4]