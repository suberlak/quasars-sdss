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

import numpy as np 
import matplotlib.pyplot as plt


javelin_chelsea_comparison = 's82drw/javelin_chelsea_comparison3.txt'
lc_stats_file = 'qso_name_timespan_nobs_2.npy'

javch =  np.loadtxt(javelin_chelsea_comparison, dtype='str')
lc_stats = np.load(lc_stats_file)

print 'We have stats data for', len(lc_stats[0,:]), ' quasars.'
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

print 'We are now retrieving stats data for each matched quasar...'

# initiate stats arrays

timespan_obs = np.zeros_like(ra_jav)
nobs_object = np.zeros_like(ra_jav)
lc_length = np.zeros_like(ra_jav)
avg_N_day = np.zeros_like(ra_jav)
avg_mag_ttl = np.zeros_like(ra_jav)
avg_err_ttl = np.zeros_like(ra_jav)
avg_mjd_diff = np.zeros_like(ra_jav)
mean_time_bet_obs = np.zeros_like(ra_jav)

for i in range(len(qso_name)):
    qso = qso_name[i]
    #print qso
    for j in range(len(lc_stats[0,:])):
        #print j
        if (lc_stats[0,j][0:18] == qso) :
     #       print lc_stats[0,j], qso
            timespan_obs[i]=lc_stats[1,j]
            nobs_object[i]=lc_stats[2,j]
            lc_length[i]=lc_stats[3,j]
            avg_N_day[i]=lc_stats[4,j]
            avg_mag_ttl[i] = lc_stats[5,j]
            avg_err_ttl[i]=lc_stats[6,j]
            avg_mjd_diff[i]=lc_stats[7,j]
            mean_time_bet_obs[i]=lc_stats[8,j]

stat_names=['timespan_obs','nobs_object','lc_length', 'avg_N_day',\
 'avg_mag_ttl', 'avg_err_ttl', 'avg_mjd_diff','mean_time_betw_obs']
 
stat_vals = [timespan_obs,nobs_object,lc_length,avg_N_day,avg_mag_ttl,avg_err_ttl,\
avg_mjd_diff,mean_time_bet_obs]



for i in range(len(stat_names)):
    plt.clf()
    plt_title  = 'log(Tau ratio) vs log(Sigma ratio) vs' + stat_names[i]
    plt.title(plt_title)
    #plt.xlim((-10,20))
    #plt.ylim((-10,10))
    plt.xlabel('log(tau_jav / tau_ch)')
    plt.ylabel('log(sig_jav / sig_ch)')            
    x=np.log10(tau_jav/tau_ch)
    y=np.log10(sig_rat)
    z=stat_vals[i]
    plt.scatter(x,y,c=z,cmap=plt.cm.coolwarm)
    # plt.show()
    fname='qso_jav_ch_corr_'+str(i)+'_'+stat_names[i]+'.png'
    plt.savefig(fname)
    i=i+1