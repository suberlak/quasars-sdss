# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:53:08 2015

@author: astronomy

Program to plot Structure Function given the flux difference vs tau table. 

It does the following : 

* Read in the data from master file (allowing to select whether star or 
  qso for column length)
* Allow some selection criteria, eg. use only those objects that have <mag> 
  below certain threshold, or a given error value,  or redshift... 
* Plot the points  delta_mag vs  time_difference,  overplotting  the mean of 
  delta_mag  calculated in each bin, and RMS calculated in two ways 
* Plot  the Structure Function,  taking the RMS  from each bin, and plotting 
  it versus the time difference  


"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.stats import binned_statistic

# Read in the LC files 
inDir =  './sf_TRY/'
outDir = './sf_TRY/'
choice = 'stars' 
if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 

if choice == 'stars' : 
    inFile = 'SF_CRTS_stars_master.txt'
else:
    inFile = 'SF_CRTS_quasars_master.txt'
    
raw_data = np.loadtxt(inDir+inFile, dtype='str')

delflx  = raw_data[:,0].astype(float)
tau     = raw_data[:,1].astype(float)
avg_mag = raw_data[:,2].astype(float)
avg_err = raw_data[:,3].astype(float)
ID      = raw_data[:,4].astype(float)
# color (not in CRTS...)









####
#### ALLOW SOME SELECTION  CRITERIA HERE.... 
####








###################
#### PLOTTING  ####
###################

def plotting(tau,delflx) :   
    plt.clf()
    plt.scatter(tau, delflx)
    plt.xlabel('Time difference [days]')
    plt.ylabel(r'$\Delta$ m')
    

# Define functions for bin statistics 

rms_robust = lambda x : np.percentile(x,75) - np.percentile(x,25)
rms_std = lambda x : np.std(x)
nbins = 100 # ensure uniform sampling in all statistics (same bins...)

# Calculate bin statistics 
bin_means = binned_statistic(tau, delflx, statistic = 'mean', bins=nbins)[0]
bin_rms_std = binned_statistic(tau, delflx, statistic = rms_std, bins=nbins)[0]
bin_rms_robust = binned_statistic(tau, delflx, statistic = rms_robust, bins=nbins)[0]

# Pull out some tau to plot means... 
bin_tau = binned_statistic(tau, tau, statistic='mean', bins=nbins)[0]

# Quickly plot the three lines : and Zeljko meant to plot those, and +/-
plotting(tau,delflx)

plt.plot(bin_tau, bin_means, color='white', label='Mean', lw = 2)
plt.plot(bin_tau, bin_means + bin_rms_std, color='r',lw = 2, label='Mean+/-RMS_std')
plt.plot(bin_tau, bin_means - bin_rms_std, color='r',lw = 2)
plt.plot(bin_tau, bin_means - bin_rms_robust, color='yellow',lw = 2 ,label='Mean+/-RMS_robust')
plt.plot(bin_tau, bin_means + bin_rms_robust, color='yellow',lw = 2)
plt.xlabel('Time difference [days]')
plt.ylabel(r'$\Delta$ m')
plt.legend()
plt.savefig('SF_tau-vs-del_mag.png')
plt.show()
