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

def sel_points(dir_in_out, fname):
    good_LC = np.loadtxt(dir_in_out + 'good_err_LC.txt', dtype='str')
    good_LC_mask = np.zeros_like(fname, dtype='bool')
    for i in range(len(fname)):
        obj_compared = fname[i][4:]
        print '\nComparison in progress...', str((float(i) / float(len(fname)) )*100.0)[:5], '%'
        for name in good_LC :
            if  obj_compared == name[4:-8] :
                good_LC_mask[i] = True
        
    print 'Out of ', len(fname), 'objects, we use ',  good_LC_mask.sum()
    return good_LC_mask
    
good_LC_mask = sel_points(dir_in_out[ch], fname)