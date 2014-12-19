# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 22:40:26 2014

@author: astronomy


same as qso_crts_preprocessing ,  but this one ignores days with N=1 , 
and calculates reduced chi-squared, and reports what percentage of days 
was thus ignored  

"""

import numpy as np
import os 


dir_in='QSO_CRTS/'
dir_out = 'QSO_CRTS_processed_N2/'
 
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
# IGNORE THOSE DAYS WITH N=1, SO THAT REDUCED  CHI-2  CAN BE CALCULATED 
#


processed_files = np.empty(0,dtype=str) # store names of processed files  : 
                                        # makes automatically a filelist :) 
percent_N = np.empty(0,dtype=float)
good_days = np.empty(0,dtype=float)
total_days = np.empty(0,dtype=float)
chi2mean_arr = np.empty(0,dtype=float)
chi2rms_arr = np.empty(0,dtype=float)
counter=0
cond = np.zeros_like(qso_chosen,dtype=bool)

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
 
    # check how long should be the storage array
    n_count = 0.0 
    which_day = np.zeros_like(days,dtype=bool)
    
    
    # figure out which days have more than one measurement (N>1)    
    for i in range(len(days)):
        day = days[i]
        int_mjd = np.require(mjd,int)       # forcing mjd array -> integers
        condition = (int_mjd == day)        # finding where int(mjd) = day
        N = float(len(mags[condition]))            # number of obs in that night 
        if N > 1 : 
            n_count += 1.0
            which_day[i] = True
    assert which_day.sum() == n_count
    
    percent= (n_count / len(days)  ) * 100.0
   
    print 'Out of ', len(days), 'days for that QSO, only ', n_count , ' had N>1.'
    print 'This is ', str(percent)[:5], '% of measurement for that QSO' 
    
    long_days = days[which_day]
     
    # storage arrays  for each qso 
    mjd_arr = np.zeros(n_count).astype(float)
    avg_mags = np.zeros(n_count).astype(float)
    avg_err_weights = np.zeros(n_count).astype(float)
    avg_err_rms = np.zeros(n_count).astype(float)
    chi2arr = np.zeros(n_count).astype(float)
    Nobs = np.zeros(n_count).astype(float)
    
    print 'obj= ',counter, 'For Quasar', obj

    # loop through days calculating mean, etc. 
    for i in range(int(n_count)):
        day = long_days[i]   # only use days with N>1
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
        rms = np.sqrt(np.mean(np.power( avgmag - mags[condition], 2.0)))
        error_rms = rms 
        
        # increase error if too small 
        if error_rms < 0.02 : 
            error_rms = np.sqrt(error_rms**2.0 + 0.01**2.0)
                  
        if error_weights < 0.02 : 
            error_weights = np.sqrt(error_weights**2.0 + 0.01**2.0)  
            
        avg_err_weights[i] = error_weights
        avg_err_rms[i] = error_rms
        mjd_arr[i] = np.mean(mjd[condition])
        chi2 = np.sum(weights*(np.power((mags[condition]-avgmag),2.0))) / (N-1)
        chi2arr[i] = chi2
        # print 'i = ', i, 'On day MJD', day, 'N obs=', N, 'avgmag=', avgmag, \
        # 'avg_err=',error, 'chi2=',chi2
      
    
    # save output of averaging of each file to a separate file 
    
    #assert len(mjd_arr) == len(avg_mags) == len(avg_err_weights) == \
    #len(avg_err_rms) == len(Nobs) == len(chi2arr)
    
    
    # Only take those lightcurves longer than 10 entries 
    if len(mjd_arr) >  10 :
    # Calculate mean and rms of the reduced chi2   
        chi2mean = np.mean(chi2arr)
        chi2rms  = np.sqrt(np.mean(np.power(chi2mean - chi2arr,2.0)))
        chi2mean_arr = np.append(chi2mean_arr, chi2mean)
        chi2rms_arr = np.append(chi2rms_arr, chi2rms)  
         
        name_out=dir_out+'out_'+obj[:18]+'.txt'
        print 'Saved', name_out
        np.savetxt(name_out, np.column_stack((mjd_arr,avg_mags,avg_err_weights, avg_err_rms,Nobs,\
    chi2arr)),fmt='%11.4f')
        name_out=dir_out+'out_'+obj[:18]+'.txt'
        
        # Global  stats 
        processed_files = np.append(processed_files, name_out[len(dir_out):])
        percent_N = np.append(percent_N, percent)
        good_days = np.append(good_days , n_count)
        total_days = np.append(total_days, len(days))
        cond[counter] = True

    counter += 1    
    
    print '  ' 
    
total_good_percent = (good_days.sum() / total_days.sum()) * 100.0

print 'In total we have ', str(total_good_percent)[:5], \
'% days which had more than one measurement for all ',len(qso_chosen[cond]),\
' quasars considered'

np.savetxt(dir_out+'chi_stats_out_files.txt', np.column_stack((qso_chosen[cond],\
 good_days, total_days, chi2mean_arr, chi2rms_arr)), fmt='%s')
 
np.savetxt(dir_out+'out.list',processed_files,delimiter = ' ', newline='\n', fmt='%s' )