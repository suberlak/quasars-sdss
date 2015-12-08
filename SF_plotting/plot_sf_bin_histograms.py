# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 14:04:43 2015

@author: suberlak
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['ytick.labelsize'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['axes.labelsize'] = 30
rcParams['axes.linewidth'] = 2
rcParams['font.size'] = 20

N = np.genfromtxt('QSO_bin_size_18.5-19.txt')


x = np.arange(1,201)

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(111)
ax.plot(x,N)
ax.set_title('QSO 18.5-19 mag')
plt.tight_layout()
plt.savefig('QSO_hist_18.5-19.png')