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

if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 

'''
Must make a choice between stars and quasars , and determine sample name 
(deafult  : s0)
'''

choice = 'qso' 
sample = 's0'

if choice == 'stars' : 
    inFile = 'SF_CRTS_stars_master.txt'
if choice == 'qso' :
    inFile = 'SF_CRTS_quasars_sample.txt'
    
raw_data = np.loadtxt(inDir+inFile, dtype='str')

delflx  = raw_data[:,0].astype(float)
tau     = raw_data[:,1].astype(float)
avg_mag = raw_data[:,2].astype(float)
avg_err = raw_data[:,3].astype(float)

if choice =='stars' : 
    ID      = raw_data[:,4].astype(float)  # as floats  
if choice == 'qso' :
    ID      = raw_data[:,4]  # as strings  
# color (not in CRTS...)




####
#### ALLOW SOME SELECTION  CRITERIA HERE.... 
####





###################
#### PLOTTING  ####
###################



def plotting(tau,delflx) :   
    plt.clf()
    
    ### A simple scatter plotr  
    #plt.scatter(tau, delflx)     
    
    ### plotting in a MY HISTOGRAM  way - too cluttered...     
    H, xedges,yedges = np.histogram2d(tau,delflx,bins=50)  
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    plt.pcolormesh(xedges, yedges, Hmasked)   # as  a color map 
   
   ## from KDE  : EXTREMELY SLOW !  
     # from scipy.stats import gaussian_kde
#    x = tau 
#    y = delflx
#    xy = np.vstack([x,y])  
#    z = gaussian_kde(xy)(xy)
#    idx = z.argsort()
#    x, y, z = x[idx], y[idx], z[idx]
#    fig, ax = plt.subplots()
#    ax.scatter(x, y, c=z, s=5, edgecolor='')
    
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


# ######## #
# Quickly plot the three lines : and Zeljko meant to plot those, and +/-
# ######## # 

plotting(tau,delflx)
plt.plot(bin_tau, bin_means, color='Gold', label='Mean', lw = 2)
plt.plot(bin_tau, bin_means + bin_rms_std, color='r',lw = 2, label='Mean+/-RMS_std')
plt.plot(bin_tau, bin_means - bin_rms_std, color='r',lw = 2)
plt.plot(bin_tau, bin_means - bin_rms_robust, color='Magenta',lw = 2 ,label='Mean+/-RMS_robust')
plt.plot(bin_tau, bin_means + bin_rms_robust, color='Magenta',lw = 2)
plt.title(r'Flux difference vs tau, nbins ='+str(nbins)+', '+choice) 
plt.xlabel('Time difference [days]')
plt.ylabel(r'$\Delta$ m')
plt.legend()
title1 = 'SF_'+choice+'_tau-vs-del_mag_'+sample+'.png'
plt.savefig(title1)
plt.show()

# ############
# Plot  the rms_std vs tau - i.e. the Structure Function 
# ############

# nbins above, check also how many stars are in the sample :
N_objects = len(np.unique(ID))

plt.clf()
plt.scatter(np.log10(bin_tau), bin_rms_std)
plt.xlabel('$log_{10}$ [days]')
plt.ylabel('SF = standard deviation')
plt.title('SF '+choice+', '+str(N_objects)+' objects')
title2 = 'SF_'+choice+'_tau-vs-rms_std_'+sample+'.png'
plt.savefig(title2)
plt.show()
# problem - it doesn't look how I expected it to look... 

# ############
# Plot the rms_robust vs tau  
# ############

plt.clf()
plt.scatter(np.log10(bin_tau), bin_rms_robust)
plt.xlabel('$log_{10}$ [days]')
plt.ylabel('SF = (75% - 25%) ')
plt.title('SF '+choice+', '+str(N_objects)+' objects')
title2 = 'SF_'+choice+'_tau-vs-rms_robust_'+sample+'.png'
plt.savefig(title2)
plt.show()