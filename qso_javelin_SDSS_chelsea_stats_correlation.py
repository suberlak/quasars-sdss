# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 18:04:33 2014

@author: suberlak

A program to take the lightcurve stats from  qso_lc_stats.py ,  
and the output of comparison of javelin and Chelsea's result for each quasar. 
It pulls from the stats table all the relevant stats for the subset of 
quasars on which javelin was run, to plot the log(tau_difference) vs log(sigma_difference),
and colour - code the third dimension according to the value of the investigated
stat parameter for a given quasar,  in a quest of looking for a correlation 
that may show to what types of quasar lightcurves javelin is failing, and if
it is  failing systematically. 


"""
from math import e, pi, isinf
import numpy as np 
import matplotlib.pyplot as plt

dir_in = 's82drw/'
band=['u','g','r','i','z']
SDSS_javelin_files = 'javelin_SDSS_chelsea_comparison_'+band[1]+'_band.txt'
CRTS_javelin_files = 'javelin_chelsea_comparison_ALL.txt'
javelin_chelsea_comparison =dir_in +SDSS_javelin_fits

javch =  np.loadtxt(javelin_chelsea_comparison, dtype='str')

print 'We also have javelin-matched Chelsea results for ',len(javch), ' quasars.' 

qso_name = javch[:,0]
ra_jav = javch[:,1].astype(np.float)
ra_ch= javch[:,2].astype(np.float)
dec_jav =javch[:,3].astype(np.float)
dec_ch = javch[:,4].astype(np.float)
tau_jav = javch[:,5].astype(np.float)
tau_ch = javch[:,6].astype(np.float)
sigma_jav =javch[:,7].astype(np.float)
sigma_ch = javch[:,8].astype(np.float)
sig_rat = javch[:,9].astype(np.float)


log_tau_ratio = np.log10(tau_jav / tau_ch)
log_sigma_ratio = np.log10(sigma_jav/sigma_ch)

# good_indices= np.where(map(isinf,log_tau_ratio)==False)


# get rid of infinities : 
#for i in range(len(log_sigma_ratio)):
#    if map(isinf,log_sigma_ratio)[i] == True : 
#        log_sigma_ratio[i] = -1e20
#
#for i in range(len(log_tau_ratio)):
#    if map(isinf,log_tau_ratio)[i] == True : 
#        log_tau_ratio[i] = 1e20
        
st = np.asarray(map(isinf,log_sigma_ratio),dtype=bool) 
lt = np.asarray(map(isinf,log_tau_ratio),dtype=bool) 
res = st+lt

ist = np.where(st == True)
ilt= np.where(lt ==True)

ires = np.where(res == True) # this np.boolean array should contain indices 
     # for which both st and lt  are True, i.e. either log_sigma_ratio  or log_tau_ratio  is infinite 

gi = -res   # good_indices 

non_inf = len(np.where(gi == True)[0])

print 'Out of ', len(ra_jav),' rows, we have ', non_inf, ' of those that do not',\
' have any infinities, and only those are saved in the txt output '

print 
DAT = np.column_stack((ra_jav[gi], ra_ch[gi], dec_jav[gi],dec_ch[gi], tau_jav[gi],\
tau_ch[gi], log_tau_ratio[gi], sigma_jav[gi], sigma_ch[gi], sig_rat[gi], \
log_sigma_ratio[gi], timespan_obs[gi], nobs_object[gi], lc_length[gi], \
avg_N_day[gi], avg_mag_ttl[gi], avg_err_ttl[gi], avg_mjd_diff[gi], mean_time_bet_obs[gi] )) 

# sort the DAT column accoring to col(7) : log_tau_ratio : 
 
newDAT=DAT[DAT[:,6].argsort()]

output = 'qso_jav_chelsea_correlation_tabulated_ALL.txt'
np.savetxt(output,newDAT,fmt="%s")


