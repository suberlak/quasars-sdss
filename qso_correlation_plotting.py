# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 15:55:24 2014

@author: suberlak

Plotting histograms of correlated values from 
qso_javelin_chelsea_stats_correlation.py

"""
import numpy as np
import matplotlib.pyplot as plt

results  = 'qso_jav_chelsea_correlation_tabulated.txt'
output =  np.loadtxt(results, dtype='str')


ra_jav = output[:,0].astype(np.float)
ra_ch = output[:,1].astype(np.float)
dec_jav= output[:,2].astype(np.float)
dec_ch= output[:,3].astype(np.float)
tau_jav= output[:,4].astype(np.float)
tau_ch= output[:,5].astype(np.float)
log_tau_ratio= output[:,6].astype(np.float)
sigma_jav= output[:,7].astype(np.float)
sigma_ch= output[:,8].astype(np.float)
sig_rat= output[:,9].astype(np.float)
log_sigma_ratio= output[:,10].astype(np.float)
timespan_obs= output[:,11].astype(np.float)
nobs_object= output[:,12].astype(np.float)
lc_length= output[:,13].astype(np.float)
avg_N_day= output[:,14].astype(np.float)
avg_mag_ttl= output[:,15].astype(np.float)
avg_err_ttl= output[:,16].astype(np.float)
avg_mjd_diff= output[:,17].astype(np.float)
mean_time_bet_obs= output[:,18].astype(np.float) 


divide = np.max(np.where(log_tau_ratio < -2))


#  AVG ERR TTL 

# figuring out what are the percentiles of the avg_err_ttl distribution

pct1sig = len(avg_err_ttl) * np.array([0.25,0.50,0.75])  
medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
vsort = np.sort(avg_err_ttl)  # sorts the array along either sigma or tau dimension 
hpd = vsort[medlowhig] # picks out values at the positions for the 
print 'avg_err_ttl   25, 50, 75 percentiles ', hpd 

hist_min = hpd[1] - 2*(hpd[2]-hpd[0])
hist_max = hpd[1] +2*(hpd[2]-hpd[0])  
print hist_max, hist_min

plt.clf()
plt.hist(avg_err_ttl[0:divide],bins=40)
plt.xlim((hist_min,hist_max))
plt.title('Average error distribution, bad tau ratio (<-2)')
plt.xlabel('Average total error ')
plt.ylabel('Frequency')
fname1 = 'qso_corr-bad_tau-avg_ttl_err-hist.png'
plt.savefig(fname1)

plt.clf()
plt.hist(avg_err_ttl[divide:],bins=40)
plt.xlim((hist_min,hist_max))
plt.title('Average error distribution, good tau ratio (>-2)')
plt.xlabel('Average total error ')
plt.ylabel('Frequency')
fname2 = 'qso_corr-good_tau-avg_ttl_err-hist.png'
plt.savefig(fname2)

print 'plotting avg_err_ttl for bad tau vs log_sigma_ration as a  colour-coded histogram'

plt.clf()
fig1 = plt.figure()
x=avg_err_ttl[0:divide]
y=log_sigma_ratio[0:divide]
plt.plot(x,y,'.r')
plt.xlim((hist_min,hist_max))
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
nbins =100
H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlim((hist_min,hist_max))
plt.title('Sample with log (tau ratio) <-2')
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname3='qso_corr-bad_tau-avg_ttl_err-log_sig_ratio_hist.png'
plt.savefig(fname3)


print 'plotting avg_err_ttl  for bad tau vs mean_time_betw_obs as a colour-coded histogram'

plt.clf()
fig1 = plt.figure()
x=avg_err_ttl[0:divide]
y=mean_time_bet_obs[0:divide]
plt.plot(x,y,'.r')
plt.xlim((hist_min,hist_max))
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Mean time between obs  ')
nbins =100
H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlim((hist_min,hist_max))
plt.title('Sample with log (tau ratio) <-2')
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Mean time between obs  ')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname3='qso_corr-bad_tau-avg_ttl_err-mean_time_bet_obs_hist.png'
plt.savefig(fname3)



print 'plotting avg_err_ttl for good tau vs log_sigma_ration as a  colour-coded histogram'

plt.clf()
fig1 = plt.figure()
x=avg_err_ttl[divide:]
y=log_sigma_ratio[divide:]
plt.plot(x,y,'.r')
plt.xlim((hist_min,hist_max))
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
nbins =100
H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlim((hist_min,hist_max))
plt.title('Sample with log (tau ratio) >-2')
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname3='qso_corr-good_tau-avg_ttl_err-log_sig_ratio_hist.png'
plt.savefig(fname3)


print 'plotting avg_err_ttl  for good tau vs mean_time_betw_obs as a colour-coded histogram'

plt.clf()
fig1 = plt.figure()
x=avg_err_ttl[divide:]
y=mean_time_bet_obs[divide:]
plt.plot(x,y,'.r')
plt.xlim((hist_min,hist_max))
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Mean time between obs  ')
nbins =100
H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlim((hist_min,hist_max))
plt.title('Sample with log (tau ratio) >-2')
plt.xlabel('Average of quasar error on magnitude')
plt.ylabel('Mean time between obs  ')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname3='qso_corr-good_tau-avg_ttl_err-mean_time_bet_obs_hist.png'
plt.savefig(fname3)



#  MEAN TIME BET OBS 

print 'Plotting histograms for mean time bet obs '

pct1sig = len(mean_time_bet_obs) * np.array([0.25,0.50,0.75])  
medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
vsort = np.sort(mean_time_bet_obs)  # sorts the array along either sigma or tau dimension 
hpd = vsort[medlowhig] # picks out values at the positions for the 
print 'mean_time_bet_obs 25, 50, 75 percentiles ', hpd 

hist_min = hpd[1] - 2*(hpd[2]-hpd[0])
hist_max = hpd[1] +2*(hpd[2]-hpd[0])  
print hist_max, hist_min

plt.clf()
plt.hist(mean_time_bet_obs[0:divide],bins=40)
plt.xlim((hist_min,hist_max))
plt.title('mean_time_bet_obs distribution, bad tau ratio (<-2)')
plt.xlabel('mean_time_bet_obs ')
plt.ylabel('Frequency')
fname1 = 'qso_corr-bad_tau-mean_time_bet_obs-hist.png'
plt.savefig(fname1)

plt.clf()
plt.hist(mean_time_bet_obs[divide:],bins=40)
plt.xlim((hist_min,hist_max))
plt.title('mean_time_bet_obs, good tau ratio (>-2)')
plt.xlabel('mean_time_bet_obs ')
plt.ylabel('Frequency')
fname2 = 'qso_corr-good_tau-mean_time_bet_obs-hist.png'
plt.savefig(fname2)


