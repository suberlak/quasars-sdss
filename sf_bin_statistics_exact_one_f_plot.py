# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 14:54:18 2015

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


ch_sig, full_sig = np.genfromtxt('QSO_17-19_err_0.3_no_corr_bin_sigma_mode_col1-chunks_col2-full.txt')
bin_max = 200

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(111)
#fig.subplots_adjust(hspace=0)
#ax = axs.ravel()

x = np.arange(1,bin_max+1)
ax.set_ylabel(r'SF ($=\sigma$)')
ax.set_xlabel('Bin number')
ax.plot(x,ch_sig, label='chunks', lw=2)
ax.plot(x,full_sig, label='full', lw=2)
ax.legend(framealpha=0.7, loc = 'upper left')
title = 'QSO_mags_17-19_err_0.3_no_corr_bins_1-'+str(bin_max)+'.png'
plt.savefig(title)