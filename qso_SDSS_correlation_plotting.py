# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 15:55:24 2014

@author: suberlak

Plotting histograms of javelin SDSS vs Chelsea SDSS from  
chelsea_results_SDSS_load.py

"""
import numpy as np
import matplotlib.pyplot as plt
from math import  isinf
band = 'u'

# set_prior  = FALSE  in Javelin  
#results  = 's82drw/javelin_SDSS_chelsea_comparison_'+band+'_band_TEST1_with_z.txt'  # TEST is with NO PRIOR

# set_prior = TRUE in Javelin
results='s82drw/javelin_SDSS_chelsea_comparison_u_band_with_z.txt'  # standard setting is WITH PRIOR
#results  = 's82drw/javelin_SDSS_chelsea_comparison_'+band+'_band_MEAN_SUB_with_z.txt' # no prior, mean-subtracted fitting, with a column for z 
fname3='sdss_'+band+'_band-log_sigma_ratio_vs_log_tau_ratio_NO_PRIOR_MEAN_SUB.png'
fname4='sdss_'+band+'_band-log_sigma_ratio_vs_log_tau_ratio_NO_PRIOR_MEAN_SUB_redshift.png'
fname5='sdss_'+band+'_band-log_sigma_vs_log_tau_CHELSEA.png'
#fname6='sdss_'+band+'_band-log_sigma_vs_log_tau_NO_PRIOR_MEAN_SUB.png'

output =  np.loadtxt(results, dtype='str')


ra_jav = output[:,1].astype(np.float)
ra_ch = output[:,2].astype(np.float)
dec_jav= output[:,3].astype(np.float)
dec_ch= output[:,4].astype(np.float)
tau_jav= output[:,5].astype(np.float)
tau_ch= output[:,6].astype(np.float)
sigma_jav= output[:,7].astype(np.float)
sigma_ch= output[:,8].astype(np.float)
sig_rat= output[:,9].astype(np.float)
redshift = output[:,10].astype(np.float)


# HISTOGRAM OF  TAU AND SIGMA RATIOS   : no redshift correction

 
# 
#print '\n Plotting coloured hist for log_tau_ratio  vs log_sigma_ratio' 
#plt.clf()
#fig1 = plt.figure()
#x = np.log10(tau_jav / tau_ch)
#y = np.log10(sigma_jav / sigma_ch)
#
#xinf = np.asarray(map(isinf,x),dtype=bool)
#yinf = np.asarray(map(isinf,y),dtype=bool)
#ttlinf = xinf + yinf
## ttlwh = np.where(ttlinf == True)  list of good indices
#gi = -ttlinf  # good_indices 
#non_inf = len(np.where(gi == True)[0])
#print 'Out of ', len(ra_jav),' rows, we have ', non_inf, ' of those that do not',\
#' have any infinities, and only those are used for plotting '
#
#plt.plot(x[gi],y[gi],'.r')
#nbins =150
#
#H, xedges,yedges = np.histogram2d(x[gi],y[gi],bins=nbins)
#H = np.rot90(H)
#H = np.flipud(H)
#Hmasked = np.ma.masked_where(H==0,H)
#fig2 = plt.figure()
#plt.pcolormesh(xedges, yedges, Hmasked)
#plt.xlim((-2,2))
#plt.ylim((-1.5,1.5))
#
#title1 =  'SDSS quasars : Javelin vs Chelsea, '+band+' band for '+str(non_inf)+ ' objects'
#plt.title(title1)
#plt.axhline(0)
#plt.axvline(0)
#plt.xlabel('log(tau_jav / tau_ch)')
#plt.ylabel('log(sigma_jav / sigma_ch)')
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')
#plt.savefig(fname3)
#     

# HISTOGRAM OF  TAU AND SIGMA RATIOS   : with redshift correction       

#print '\n Plotting coloured hist for log_tau_ratio  vs log_sigma_ratio with redshift correction' 
#plt.clf()
#fig1 = plt.figure()
#x = np.log10(tau_jav / ((1+redshift)*tau_ch))
#y = np.log10(sigma_jav / sigma_ch)
#
#xinf = np.asarray(map(isinf,x),dtype=bool)
#yinf = np.asarray(map(isinf,y),dtype=bool)
#ttlinf = xinf + yinf
## ttlwh = np.where(ttlinf == True)  list of good indices
#gi = -ttlinf  # good_indices 
#non_inf = len(np.where(gi == True)[0])
#print 'Out of ', len(ra_jav),' rows, we have ', non_inf, ' of those that do not',\
#' have any infinities, and only those are used for plotting '
#
#plt.plot(x[gi],y[gi],'.r')
#nbins =200
#
#H, xedges,yedges = np.histogram2d(x[gi],y[gi],bins=nbins)
#H = np.rot90(H)
#H = np.flipud(H)
#Hmasked = np.ma.masked_where(H==0,H)
#fig2 = plt.figure()
#plt.pcolormesh(xedges, yedges, Hmasked)
#plt.xlim((-2,2))
#plt.ylim((-1.5,1.5))
#title = 'SDSS quasars : Javelin vs Chelsea, '+band+' band for '+str(non_inf)+ ' objects \n redshift corrected '
#plt.title(title)
#plt.axhline(0)
#plt.axvline(0)
#plt.xlabel('log(tau_jav / tau_ch)')
#plt.ylabel('log(sigma_jav / sigma_ch)')
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')
#plt.savefig(fname4)
##     
#
#print ' Files saved are',fname3, ' and ', fname4
# 
 
# HISTOGRAM  OF   LOG SIGMA  VS LOG TAU,  LIKE FIG 13 FROM MCLEOD+2011 


print '\n Plotting coloured hist for log_tau  vs log_sigma  for javelin fitting ' 
plt.clf()
fig1 = plt.figure()
x = np.log10(sigma_ch)
y = np.log10(tau_ch)

xinf = np.asarray(map(isinf,x),dtype=bool)
yinf = np.asarray(map(isinf,y),dtype=bool)
ttlinf = xinf + yinf
# ttlwh = np.where(ttlinf == True)  list of good indices
gi = -ttlinf  # good_indices 
non_inf = len(np.where(gi == True)[0])
print 'Out of ', len(ra_jav),' rows, we have ', non_inf, ' of those that do not',\
' have any infinities, and only those are used for plotting '

plt.plot(x[gi],y[gi],'.r')
nbins =100

H, xedges,yedges = np.histogram2d(x[gi],y[gi],bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlim((-1.5,0.1))
plt.ylim((-1,4.3))
title = 'SDSS quasars : Chelsea results,  '+band+' band for '+str(non_inf)+ ' objects '
plt.title(title)
#plt.axhline(0)
#plt.axvline(0)
plt.ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15)
plt.xlabel(r'$\log_{10}{ \, \sigma_{ch}}$',fontsize=15)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
plt.savefig(fname5)
print 'File saved is ', fname5 
 
 
### SINGLE HISTOGRAMS 
#
#for i in range(5,6):
#    print '\n Plotting histograms for',names[i]
#    pct1sig = len(stats[i][0:divide]) * np.array([0.25,0.50,0.75])  
#    medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
#    vsort = np.sort(stats[i][0:divide])  # sorts the array along either sigma or tau dimension 
#    hpdb = vsort[medlowhig] # picks out values at the positions for the 
#    print names[i], ': 25, 50, 75 percentiles for bad tau ', hpdb 
#    
#    hist_minb = hpdb[1] - 2*(hpdb[2]-hpdb[0])
#    hist_maxb = hpdb[1] +2*(hpdb[2]-hpdb[0])  
#    print  'hist_min, hist_max', hist_minb, hist_maxb
#    
#    
#    pct1sig = len(stats[i][divide:]) * np.array([0.25,0.50,0.75])  
#    medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
#    vsort = np.sort(stats[i][divide:])  # sorts the array along either sigma or tau dimension 
#    hpdg = vsort[medlowhig] # picks out values at the positions for the 
#    print names[i], ': 25, 50, 75 percentiles for good tau', hpdg 
#    
#    hist_ming = hpdg[1] - 2*(hpdg[2]-hpdg[0])
#    hist_maxg = hpdg[1] +2*(hpdg[2]-hpdg[0])  
#    print  'hist_min, hist_max', hist_ming, hist_maxg
#    
#    
#    plt.clf()
#    plt.hist(stats[i][0:divide],bins=50)
#    plt.xlim((hist_minb,hist_maxb))
#    plot_title = names[i] + 'distribution, bad tau ratio (<-2.5)'
#    plt.title(plot_title)
#    plt.xlabel(names[i])
#    plt.ylabel('Frequency')
#    fname1 = 'corr_h1d-'+names[i]+'-bad_tau.png'
#    plt.savefig(fname1)
#    
#    plt.clf()
#    plt.hist(stats[i][divide:],bins=50)
#    plt.xlim((hist_ming,hist_maxg))
#    plot_title = names[i] + 'distribution, good tau ratio (>-2.5)'
#    plt.title(plot_title)
#    plt.xlabel(names[i])
#    plt.ylabel('Frequency')
#    fname2 = 'corr_h1d-'+names[i]+'-good_tau.png'
#    plt.savefig(fname2)

# COLOURED HISTOGRAMS     
    
#for i in range(len(stats)):
#    
#    pct1sig = len(stats[i][0:divide]) * np.array([0.25,0.50,0.75])  
#    medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
#    vsort = np.sort(stats[i][0:divide])  # sorts the array along either sigma or tau dimension 
#    hpdb = vsort[medlowhig] # picks out values at the positions for the 
#    print '\n', names[i], ': 25, 50, 75 percentiles for bad tau ', hpdb 
#    
#    hist_minb = hpdb[1] - 2*(hpdb[2]-hpdb[0])
#    hist_maxb = hpdb[1] +2*(hpdb[2]-hpdb[0])  
#    print  'hist_min, hist_max', hist_minb, hist_maxb
#    
#    
#    pct1sig = len(stats[i][divide:]) * np.array([0.25,0.50,0.75])  
#    medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
#    vsort = np.sort(stats[i][divide:])  # sorts the array along either sigma or tau dimension 
#    hpdg = vsort[medlowhig] # picks out values at the positions for the 
#    print '\n',names[i], ': 25, 50, 75 percentiles for good tau', hpdg 
#    
#    hist_ming = hpdg[1] - 2*(hpdg[2]-hpdg[0])
#    hist_maxg = hpdg[1] +2*(hpdg[2]-hpdg[0])  
#    print  'hist_min, hist_max', hist_ming, hist_maxg    
#    
#    for j in range(len(stats)):
#        if (i != j) : 
#            print '\n Plotting coloured hist for',names[i], ' vs ', names[j], 'for bad tau' 
#            plt.clf()
#            fig1 = plt.figure()
#            x=stats[i][0:divide]
#            y=stats[j][0:divide]
#            plt.plot(x,y,'.r')
#            plt.xlim((hist_minb,hist_maxb))
#            #plt.xlabel(names[i])
#            #plt.ylabel(names[j])
#            nbins =100
#            H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
#            H = np.rot90(H)
#            H = np.flipud(H)
#            Hmasked = np.ma.masked_where(H==0,H)
#            fig2 = plt.figure()
#            plt.pcolormesh(xedges, yedges, Hmasked)
#            plt.xlim((hist_minb,hist_maxb))
#            plt.title('Sample with log (tau ratio) <-2')
#            plt.xlabel(names[i])
#            plt.ylabel(names[j])
#            cbar = plt.colorbar()
#            cbar.ax.set_ylabel('Counts')
#            fname3='corr_h2d-'+names[i]+'_vs_'+names[j]+'-bad_tau.png'
#            plt.savefig(fname3)
#            
#            print '\n Plotting coloured hist for',names[i], ' vs ', names[j], 'for good tau'
#            plt.clf()
#            fig1 = plt.figure()
#            x=stats[i][divide:]
#            y=stats[j][divide:]
#            plt.plot(x,y,'.r')
#            plt.xlim((hist_ming,hist_maxg))
#            #plt.xlabel(names[i])
#            #plt.ylabel(names[j])
#            nbins =100
#            H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
#            H = np.rot90(H)
#            H = np.flipud(H)
#            Hmasked = np.ma.masked_where(H==0,H)
#            fig2 = plt.figure()
#            plt.pcolormesh(xedges, yedges, Hmasked)
#            plt.xlim((hist_ming,hist_maxg))
#            plt.title('Sample with log (tau ratio) >-2')
#            plt.xlabel(names[i])
#            plt.ylabel(names[j])
#            cbar = plt.colorbar()
#            cbar.ax.set_ylabel('Counts')
#            fname3='corr_h2d-'+names[i]+'_vs_'+names[j]+'-good_tau.png'
#            plt.savefig(fname3)
            
            
            
#  AVG ERR TTL 

# figuring out what are the percentiles of the avg_err_ttl distribution
#print '\n Plotting histograms for avg err ttl '
#
#pct1sig = len(avg_err_ttl[0:divide]) * np.array([0.25,0.50,0.75])  
#medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
#vsort = np.sort(avg_err_ttl[0:divide])  # sorts the array along either sigma or tau dimension 
#hpdb = vsort[medlowhig] # picks out values at the positions for the 
#print 'avg_err_ttl for bad tau   25, 50, 75 percentiles ', hpdb 
#
#hist_minb = hpdb[1] - 2*(hpdb[2]-hpdb[0])
#hist_maxb = hpdb[1] +2*(hpdb[2]-hpdb[0])  
#print hist_maxb, hist_minb
#
#pct1sig = len(avg_err_ttl[divide:]) * np.array([0.25,0.50,0.75])  
#medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
#vsort = np.sort(avg_err_ttl[divide:])  # sorts the array along either sigma or tau dimension 
#hpdg = vsort[medlowhig] # picks out values at the positions for the 
#print 'avg_err_ttl for good tau   25, 50, 75 percentiles ', hpdg 
#
#hist_ming = hpdg[1] - 2*(hpdg[2]-hpdg[0])
#hist_maxg = hpdg[1] +2*(hpdg[2]-hpdg[0])  
#print hist_maxg, hist_ming
#
#
#
#plt.clf()
#plt.hist(avg_err_ttl[0:divide],bins=40)
#plt.xlim((hist_minb,hist_maxb))
#plt.title('Average error distribution, bad tau ratio (<-2)')
#plt.xlabel('Average total error ')
#plt.ylabel('Frequency')
#fname1 = 'qso_corr-avg_ttl_err-bad_tau-hist.png'
#plt.savefig(fname1)
#
#
#plt.clf()
#plt.hist(avg_err_ttl[divide:],bins=40)
#plt.xlim((hist_ming,hist_maxg))
#plt.title('Average error distribution, good tau ratio (>-2)')
#plt.xlabel('Average total error ')
#plt.ylabel('Frequency')
#fname2 = 'qso_corr-avg_ttl_err-good_tau-hist.png'
#plt.savefig(fname2)

#print '\n plotting avg_err_ttl for bad tau vs log_sigma_ration as a  colour-coded histogram'
#
#plt.clf()
#fig1 = plt.figure()
#x=avg_err_ttl[0:divide]
#y=log_sigma_ratio[0:divide]
#plt.plot(x,y,'.r')
#plt.xlim((hist_minb,hist_maxb))
#plt.xlabel('Average of quasar error on magnitude')
#plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
#nbins =100
#H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
#H = np.rot90(H)
#H = np.flipud(H)
#Hmasked = np.ma.masked_where(H==0,H)
#fig2 = plt.figure()
#plt.pcolormesh(xedges, yedges, Hmasked)
#plt.xlim((hist_minb,hist_maxb))
#plt.title('Sample with log (tau ratio) <-2')
#plt.xlabel('Average of quasar error on magnitude')
#plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')
#fname3='qso_corr-avg_ttl_err-log_sig_ratio-bad_tau-hist.png'
#plt.savefig(fname3)

#print '\n plotting avg_err_ttl for good tau vs log_sigma_ration as a  colour-coded histogram'
#
#plt.clf()
#fig1 = plt.figure()
#x=avg_err_ttl[divide:]
#y=log_sigma_ratio[divide:]
#plt.plot(x,y,'.r')
#plt.xlim((hist_ming,hist_maxg))
#plt.xlabel('Average of quasar error on magnitude')
#plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
#nbins =100
#H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
#H = np.rot90(H)
#H = np.flipud(H)
#Hmasked = np.ma.masked_where(H==0,H)
#fig2 = plt.figure()
#plt.pcolormesh(xedges, yedges, Hmasked)
#plt.xlim((hist_ming,hist_maxg))
#plt.title('Sample with log (tau ratio) >-2')
#plt.xlabel('Average of quasar error on magnitude')
#plt.ylabel('Log (sigma_jav / sigma_chelsea) ')
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')
#fname3='qso_corr-avg_ttl_err-log_sig_ratio-good_tau-hist.png'
#plt.savefig(fname3)