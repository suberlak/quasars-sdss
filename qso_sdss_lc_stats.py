# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 17:39:48 2014

@author: suberlak


Retrieve statistics for each lightcurve:

avg_mag ( to see whether the observed quasar is bright or dim  )
avg_err (whether the error isn't too big )
N_obs = number of rows = lightcurve length 
t_obs = timespan of observations

Made from a combination  of  qso_lc_stats.py  which takes statistics for the CRTS objects,
and qso_sdss_javelin_input_maker.py , which reads in all the 
sdss lightcurves , and splits them according to filter bandpasses 

"""

import numpy as np

dir_choice=['QSO_try2/', 'QSO_S82/', 'QSO_SDSS_JAV/', 'QSO_SDSS_analysis/']

dir_input=dir_choice[1]
dir_output=dir_choice[3]
names=np.loadtxt(dir_input+'out_good.list',dtype=str)

ra_list = np.zeros_like(names, dtype=float)
dec_list = np.zeros_like(names, dtype=float)

qso_stats = np.zeros((len(names),18))
# name , ra , dec , and then for each  u,g,r,i,z :   lc_length , avg_mag , avg_mag_err 

mjds = [0,0,0,0,0]
mags = [0,0,0,0,0]
mag_errs = [0,0,0,0,0]
mjd_cols = [0,3,6,9,12]
mag_cols = [1,4,7,10,13]
err_cols = [2,5,8,11,14]
prefix = ['u','g','r','i','z']
    
for j in range(len(names)):   # )
 
    address=dir_input+names[j]
    data=np.loadtxt(address)
    #print '\nQuasar ', names[j]
   
    
    qso_stats[j,0] = names[j]
    qso_stats[j,1] = data[0,15]  # ra 
    qso_stats[j,2] = data[0,16]  # dec 
    
    for l in range(len(mag_cols)):
        # checking whether a given mag does not have a bad mesaurement 
         
        col = mag_cols[l]
        
        #print '\ndata column checked: ', col, 'band', prefix[l]
        cond_good_data=np.ones_like(data[:,0],dtype=bool)
        #print 'Original data length ', len(data[:,col])
        for k in range(len(data)):
            if( data[k,col] < 1.0  or data[k,col] > 40):
         #       print data[k,col], 'is < 1.0 or > 40'
                cond_good_data[k] = False
        #print 'Filtered data length', len(data[cond_good_data, col])
         
        mags[l] = data[cond_good_data,mag_cols[l]] 
        mag_errs[l] = data[cond_good_data,err_cols[l]]
 
        avg_mag = np.mean(mags[l])
        avg_mag_err = np.mean(mag_errs[l])
        lc_band_length = len(mags[l])
        #print 'length, <mag>, <mag_err> for this band:', lc_band_length, avg_mag, avg_mag_err
        p = l+1
        qso_stats[j,3*p+2] = avg_mag_err
        qso_stats[j,3*p+1] = avg_mag
        qso_stats[j,3*p] = lc_band_length
        
    
        
        
    print 'Done with quasar', j+1, 'out of', len(names), '\n'

prefix = ['u','g','r','i','z']
for i in range(1,6):
    number_short = len(np.where(qso_stats[:,3*i] < 10)[0])
    print 'We find that for',prefix[i-1], 'band, there are', number_short, 'lightcurves which are shorter than 10 rows'

filename = dir_output + 'SDSS_lc_stats_multiband.txt'
np.savetxt(filename,qso_stats,fmt="%s")
print 'Saved list of all names, ra, and dec, and u,g,r,i,z band statistics: ',\
'lightcurve length per bandpass, mean band magnitude, mean band magnitude error ', filename   