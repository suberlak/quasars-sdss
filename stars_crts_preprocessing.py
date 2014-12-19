# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 18:43:56 2014

@author: astronomy

A program to read in CRTS QSO data, and pre-process the lightcurves 
for lightcurve fitting :average-out the measurements from one day, and 
take as an error the rms of the magnitude in a day.  
If such rms scatter is less than 0.02 mag,  add 0.01 mag in quadrature
to avoid mistakes caused by unreasonably small error. 
Discard lightcurves that after merging are shorter than 10 (this means there
are less than 10 observing days : same criterion as Chelsea's code )
"""

import numpy as np
import os 


dir_in='stars_CRTS/0/'
dir_err_rms = 'stars_CRTS_processed_err_rms/'
dir_err_w = 'stars_CRTS_processed_err_w/'
 
qso_names = np.loadtxt(dir_in+'file.list',dtype=str)


#
#  IGNORE EMPTY  LIGHTCURVES AND THOSE TOO SHORT ...
#

gi = np.ones_like(qso_names, dtype=bool) # good indices 

for i in range(len(qso_names)):
    if os.stat(dir_in+qso_names[i])[6]==0:  gi[i] = False
    address=dir_in+qso_names[i]
    data=np.loadtxt(address)
    if data.size < 20.0 :  gi[i] = False

num_empty=np.sum(np.logical_not(gi))
num_notempty=np.sum(gi)

print 'Out of', len(qso_names), 'files, we have', num_notempty, 'files that were long', \
'and therefore, ', num_empty, 'almost empty files'

print 'Files that were empty:', qso_names[np.logical_not(gi)]

qso_chosen = qso_names[gi]

print '\nPerforming calculations on files with more than one measurement...'


#
# CALCULATE MEAN DAILY MAG, AND MEAN MJD, AND RMS... 
# IF RMS < 0.02, ADD 0.01 IN QUADRATURE 
#

processed_files = np.empty(0,dtype=str)

counter=0


for obj in qso_chosen:
    address=dir_in+obj
    data=np.loadtxt(address)

    averages=np.zeros(shape=(len(data),3))

    mjd = data[:,0]
    mags = data[:,1]
    errs = data[:,2]
    days = data[:,0]
    days = [int(day) for day in days]
    days = np.unique(days)          # pulling out only the unique values 
 
    # storage arrays  for each qso (all N's)
    mjd_arr = np.zeros_like(days).astype(float)
    avg_mags = np.zeros_like(days).astype(float)
    avg_err_weights = np.zeros_like(days).astype(float)
    avg_err_rms = np.zeros_like(days).astype(float)
    chi2arr = np.zeros_like(days).astype(float)
    Nobs = np.zeros_like(days).astype(float)
    
    print 'obj= ',counter, 'For Quasar', obj

    # loop through days calculating mean, etc. 
    for i in range(len(days)):
        day = days[i]
        int_mjd = np.require(mjd,int)       # forcing mjd array -> integers
        condition = (int_mjd == day)        # finding where int(mjd) = day
        N = float(len(mags[condition]))            # number of obs in that night 
        Nobs[i] = N
        avgmag = np.average(mags[condition],weights=errs[condition]) # works for single measurement too ! 
        
        # weights error         
        weights=1.0 / ( errs[condition] * errs[condition]) 
        avg_mags[i] = avgmag
        error_weights = 1.0 / np.sqrt(np.sum(weights))
        
        # rms error 
        diff= avgmag - mags[condition]
        sq = diff ** 2.0
        rms = np.sqrt(sq.mean())
        error_rms = rms 
        
        # increase error if too small 
        if error_rms < 0.02 : 
            error_rms = np.sqrt(error_rms**2.0 + 0.01**2.0)
                  
        if error_weights < 0.02 : 
            error_weights = np.sqrt(error_weights**2.0 + 0.01**2.0)  
            
        avg_err_weights[i] = error_weights
        avg_err_rms[i] = error_rms
        mjd_arr[i] = np.mean(mjd[condition])
        #chi2 = np.sum(weights*(np.power((mags[condition]-avgmag),2.0))) 
        #chi2arr[i] = chi2
        # print 'i = ', i, 'On day MJD', day, 'N obs=', N, 'avgmag=', avgmag, \
        # 'avg_err=',error, 'chi2=',chi2
    
    # save output of averaging of each file to a separate file 
    # only if  more than 10 obs per day 
    if( len(mjd_arr) >  10 ) :
        name_out1=dir_err_rms+'out_'+obj[:18]+'.txt'
        name_out2 = dir_err_w + 'out_'+obj[:18]+'.txt'
        print 'Saved', name_out1, name_out2
        assert len(mjd_arr) == len(avg_mags) == len(avg_err_weights) == \
        len(avg_err_rms) == len(Nobs)
        
        np.savetxt(name_out1, np.column_stack((mjd_arr,avg_mags, avg_err_rms)),fmt='%11.4f')
        np.savetxt(name_out2, np.column_stack((mjd_arr,avg_mags, avg_err_weights)),fmt='%11.4f')
        processed_files = np.append(processed_files, name_out1[len(dir_err_rms):])
    
    counter += 1    
    
    print '  '

# Save lists of files 
np.savetxt(dir_err_w+'out.list',processed_files,delimiter = ' ', newline='\n', fmt='%s' ) 
np.savetxt(dir_err_rms+'out.list',processed_files,delimiter = ' ', newline='\n', fmt='%s' )