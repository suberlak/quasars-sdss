# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 22:12:01 2014

@author: Chris 

Plotting stats for CRTS stars,  to see if there is anything weird in those 
that have sigma < 0 ...

I leave it open to also plot these same stats for QSO, 
because it's so much quicker, to test the program ...
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf



#args = sys.argv
#err = int(args[1])
ch = 0
dir_in_out = ['QSO_CRTS_analysis/', 'stars_CRTS_analysis/']
files = ['javelin_CRTS_chain_results_err_w.txt', 'javelin_CRTS_stars_err_w_chain_results.txt']
chain_results = dir_in_out[ch]+  files[ch] 
data = np.loadtxt(chain_results,dtype='str' )


############################
# READ IN THE CHAIN VALUES#
############################

fname = np.empty(0,dtype=float)
sigma_m = np.empty(0,dtype=float)
tau_m =  np.empty(0,dtype=float)

print 'Reading in all the values ... '

for i in range(len(data[:,0])):
   try:
       fname = np.append(fname, data[i,0])
       sigma_m = np.append(sigma_m, float(data[i,2]))
       tau_m = np.append(tau_m, float(data[i,5]))
      
   except ValueError:     
       pass

if len(sigma_m) != len(tau_m) : 
    m = min(len(tau_m), len(sigma_m))
    sigma_m = sigma_m[0:m]
    tau_m = tau_m[0:m]
    fname = fname[0:m]

assert len(sigma_m) == len(tau_m)
 
print '\nOut of ', len(data[:,0]), ' rows we were able to read in ', len(tau_m)



##############################
# LIMIT TO ONLY GOOD POINTS  #
##############################
good_LC = np.loadtxt(dir_in_out[ch] + 'good_err_LC.txt', dtype='str')
good_LC_cut = np.empty(0, dtype=str)
for i in range(len(good_LC)):
    good_LC_cut = np.append(good_LC_cut, good_LC[i][4:-4])

fname_cut = np.empty(0, dtype=str)
for i in range(len(fname)):
    fname_cut = np.append(fname_cut, fname[i][:-10])

matched_names = set(good_LC_cut).intersection(fname_cut)
matched_names = np.array(list(matched_names))

indices = np.zeros_like(fname, dtype=bool)
for i in range(len(fname)):
    indices[i] = np.where(fname_cut == matched_names[0])[0][0]

    
#good_LC_mask = sel_points_qso(dir_in_out[ch], fname)


sigma_m = sigma_m[indices]
tau_m = tau_m[indices]

# recalculate to get sigma hat 

sigma_hat = sigma_m * np.sqrt(tau_m / (2.0 * 365.0))

############################
#  READ IN STATS FOR LC'S  # 
############################

stats = np.loadtxt(dir_in_out[ch]+'LC_stats.txt', dtype='str')
lc_names = stats[:,0]

lc_names_cut = np.empty(0,dtype=str)
for i in range(len(lc_names)):
    lc_names_cut = np.append(lc_names_cut, lc_names[i][4:-4])
    
inds = np.zeros_like()

mjd_span = stats[:,1]
mag_rms = stats[:,2]
mag_mean = stats[:,3]
err_mean = stats[:,4]
N_lines = stats[:,5]

##################
# PLOTTING STATS #
##################

sigma_neg = [sigma_hat< 0.0]
sigma_pos = [sigma_hat>0]  