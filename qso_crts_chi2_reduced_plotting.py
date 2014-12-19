# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 23:46:03 2014

@author: astronomy

quick program to plot the day-averaged  rms and mean of the reduced chi-squared 
from  qso_crts_counting.py

I calculated how in each day there was a scatter of measurements around some mean,
and found the reduced chi-squared :  
chi2red = np.sum(weights*(np.power((mags[condition]-avgmag),2.0))) / (N-1)

and calculated its means and rms  :

chi2mean = np.mean(chi2arr)
chi2rms  = np.sqrt(np.mean(np.power(chi2mean - chi2arr,2.0)))

"""

import numpy as np 
import matplotlib.pyplot as plt

dir_in = 'QSO_CRTS_processed_N2/'
dir_out = 'QSO_CRTS_processed_N2/'

filein = dir_in + 'chi_stats_out_files.txt'

data = np.loadtxt(filein,dtype=str)

qso_chosen = data[:,0]
good_days = data[:,1]
total_days = data[:,2]
chi2mean_arr = data[:,3].astype(float)
chi2rms_arr = data[:,4].astype(float)
plt.clf()
x = np.arange(0,len(chi2mean_arr))
#plt.xlim((0.95,1.05))
cond = np.ones_like(chi2mean_arr, dtype=bool)
i=0
#restrict range
for chi in chi2mean_arr:
    
    if chi < 0.2 : cond[i] = False
    if chi > 300 : cond[i]= False     
    i += 1.0
    
values =  chi2mean_arr[cond]    
plt.hist(values  , bins=45, range=(min(values),10))
plt.title(r'Histogram of the  $\left<\chi ^{2} _{DOF}\right>$')
plt.xlabel(r'$\left< \chi ^{2} _{DOF}\right>$'  )
plt.ylabel('Number of objects')
plt.savefig('qso_CRTS_chi2_mean_hist.png')

plt.hist(chi2rms_arr , bins=45, range=(0,20))
plt.title(r'Histogram of the rms $\chi ^{2} _{DOF}$')
plt.xlabel(r'$ \chi ^{2} _{DOF}$'  )
plt.ylabel('Number of objects')
plt.savefig('qso_CRTS_chi2_rms_hist.png')
