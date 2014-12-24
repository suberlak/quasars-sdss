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




#args = sys.argv
#err = int(args[1])
ch = 1
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



############################
# SELECTING POINTS TO USE  #
############################

def sel_points_stars(dir_in_out, fname):
    good_LC = np.loadtxt(dir_in_out + 'good_err_LC.txt', dtype='str')
    good_LC_cut = np.empty(0, dtype=str)

    for i in range(len(good_LC)):
        good_LC_cut = np.append(good_LC_cut, good_LC[i][4:-8])
        
    good_LC_mask = np.zeros_like(fname, dtype='bool')
    for i in range(len(fname)):
        print '\nComparison in progress...', str((float(i) / float(len(fname)) )*100.0)[:5], '%'
        good_LC_mask[i] =  fname[i][4:] in  good_LC_cut 
        
        
    print 'Out of ', len(fname), 'objects, we use ',  good_LC_mask.sum()
    return good_LC_mask
    
def sel_points_qso(dir_in_out, fname):
    good_LC = np.loadtxt(dir_in_out + 'good_err_LC.txt', dtype='str')
    good_LC_cut = np.empty(0, dtype=str)

    for i in range(len(good_LC)):
        good_LC_cut = np.append(good_LC_cut, good_LC[i][4:-4])
        
    good_LC_mask = np.zeros_like(fname, dtype='bool')
    for i in range(len(fname)):
        print '\nComparison in progress...', str((float(i) / float(len(fname)) )*100.0)[:5], '%'
        good_LC_mask[i] =  fname[i][:-10] in  good_LC_cut 
        
        
    print 'Out of ', len(fname), 'objects, we use ',  good_LC_mask.sum()
    return good_LC_mask
    
if ch == 0:
    good_LC_mask = sel_points_qso(dir_in_out[ch], fname)
if ch == 1 :
    good_LC_mask = sel_points_stars(dir_in_out[ch], fname)
    
fname = fname[good_LC_mask]
sigma_m = sigma_m[good_LC_mask]
tau_m = tau_m[good_LC_mask]

sigma_hat = sigma_m * np.sqrt(tau_m / (2.0 * 365.0))

############################
#  READ IN STATS FOR LC'S  # 
############################

stats = np.loadtxt(dir_in_out[ch]+'LC_stats.txt', dtype='str')
lc_names = stats[:,0]

def check_stats_qso(fname,lc_names):
    ind = np.zeros(len(fname), dtype=int)    
    for i in range(len(fname)):
        print '\nChecking', i, ' of ', len(fname)
        for j in range(len(lc_names)): 
            if lc_names[j][4:-4] == fname[i][:-10] : 
                ind[i] = j
    return ind 
    
def check_stats_stars(fname,lc_names):
    ind = np.zeros(len(fname), dtype=int)    
    for i in range(len(fname)):
        print '\nChecking', i, ' of ', len(fname)
        for j in range(len(lc_names)): 
            if lc_names[j][4:-8] == fname[i][4:] : 
                ind[i] = j
    return ind 
    
if ch == 0 : 
   ind = check_stats_qso(fname,lc_names)
if ch == 1 :
   ind = check_stats_stars(fname, lc_names)    
    
mjd_span = stats[:,1]
mag_rms = stats[:,2]
mag_mean = stats[:,3]
err_mean = stats[:,4]
N_lines = stats[:,5]

# choose only those matching from stats - now their order matches the fname order

lc_names = lc_names[ind]
mjd_span = mjd_span[ind]
mag_rms  = mag_rms[ind]
mag_mean = mag_mean[ind]
err_mean = err_mean[ind]
N_lines = N_lines[ind]

print 'These should match: ', fname[0], lc_names[0]

DAT = np.column_stack((fname,lc_names, sigma_m, tau_m, sigma_hat, mjd_span, mag_rms, mag_mean, err_mean, N_lines))

np.savetxt(dir_in_out[ch]+'javelin_sigma_tau_plus_stats_matched_good_err.txt',DAT,delimiter=" ", fmt="%s")


#newDAT=DAT[DAT[:,2].argsort()]   # sort according to sigma_m

