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


"""
File with results of Chelsea's analysis
"""

dir_in = 's82drw/'
dir_out = 'QSO_SDSS_analysis/'
band=['u','g','r','i','z']
band_pos = 0
prior = 'no'
if prior == 'yes' : 
    SDSS_jav_Ch_table = 'javelin_SDSS_chelsea_comparison_'+band[band_pos]+'_band_with_z.txt'
else :
    SDSS_jav_Ch_table  = 'javelin_SDSS_chelsea_comparison_'+band[band_pos]+'_band_TEST1_with_z.txt'

javelin_chelsea_comparison =dir_in +SDSS_jav_Ch_table
javch =  np.loadtxt(javelin_chelsea_comparison, dtype='str')

# Define plot  range 
xmin = -0.8
xmax = 1
ymin = -1
ymax = 2

"""
File with the SDSS quasars statistics: one file for all bands 
"""

filename = 'QSO_SDSS_analysis/'+ 'SDSS_lc_stats_multiband.txt'
qso_stats =  np.loadtxt(filename) 

print 'We have ugriz stats for SDSS data from file ', filename, 'of length', len(qso_stats)
print 'We also have a comparison of Chelsea and javelin  SDSS fits for ',len(javch), ' quasars.' 

print 'Retrieving statistics for quasars from ', band[band_pos], 'band'


"""
Data from comparison of Javelin and Chelsea's fits 
"""
qso_name = javch[:,0]
ra_jav = javch[:,1].astype(np.float)
ra_ch= javch[:,2].astype(np.float)
dec_jav =javch[:,3].astype(np.float)
dec_ch = javch[:,4].astype(np.float)
tau_jav = javch[:,5].astype(np.float)
tau_ch = javch[:,6].astype(np.float)
sigma_jav =javch[:,7].astype(np.float)
sigma_ch_raw = javch[:,8].astype(np.float)
sig_rat = javch[:,9].astype(np.float)
chelsea_quasars_redshift = javch[:,10].astype(np.float)

#tau_ch = tau_ch * (1 + chelsea_quasars_redshift)
sigma_hat = sigma_ch_raw
sigma_ch = sigma_hat * np.sqrt(tau_ch / (2.0*365.0))

log_tau_ratio = np.log10(tau_jav / tau_ch)
log_sigma_ratio = np.log10(sigma_jav/sigma_ch)

#print 'We are now retrieving stats data for each matched quasar...'
#
## initiate stats arrays
#lc_length = np.zeros_like(ra_jav)
#avg_mag_ttl = np.zeros_like(ra_jav)
#avg_err_ttl = np.zeros_like(ra_jav)
#
#
#print 'There are ', len(qso_name), ' quasars to match with stat data ' 
#
#p = band_pos+1
#
#for i in range(len(qso_name)):
#    index = np.where(qso_stats[:,0] == qso_name[i].astype(float))[0][0]
#    lc_length[i] = qso_stats[index,3*p]
#    avg_mag_ttl[i] = qso_stats[index,3*p+1]
#    avg_err_ttl[i] = qso_stats[index,3*p+2]
     
            

def load_x_y(x_arr, y_arr , x_limits, y_limits):
 
    print '\n Loading x and y ... ' 
   
    x = x_arr 
    y = y_arr
    
    # Removing all the infinities 
    xinf = np.asarray(map(isinf,x),dtype=bool)
    yinf = np.asarray(map(isinf,y),dtype=bool)
    ttlinf = xinf + yinf
    gi = -ttlinf  # good_indices 
    
    assert len(gi) == len(x)
    
    # Imposing selection criteria 
    ysmall = np.where(y < y_limits[0])
    ylarge = np.where(y > y_limits[1])
    xsmall = np.where(x < x_limits[0])        
    xlarge = np.where(x > x_limits[1]) 
    gi[xsmall] = False
    gi[ysmall] = False
    gi[xlarge] = False
    gi[ylarge] = False
    non_inf = len(np.where(gi == True)[0])
        
    percent = (float(non_inf) / float(len(x)))  * 100.0
   
    print 'Out of ', len(x),' rows, we have ', non_inf, ' of those that match', \
    'the criteria of ',  x_limits[0],' < x <', x_limits[1],' and ', y_limits[0],\
    ' < y < ',y_limits[1], 'and only those are used for plotting ...  '
    
    return x[gi], y[gi], non_inf, percent
    
xlim = [xmin,xmax]
ylim = [ymin, ymax]

x,y,num, percent =  load_x_y(log_sigma_ratio,log_tau_ratio, xlim, ylim)

print num, 'Quasars meet the selection criteria, i.e.',percent,'%' 

# Finding the longest lc within those selected quasars:
#index = np.where(lc_length == max(lc_length[cond]))
#print '\nFrom those quasars, the longest is', qso_name[index][0], 'with', lc_length[index][0], 'rows.'
#print 'Its <mag> =', avg_mag_ttl[index][0], 'and <avg_err_ttl>=', avg_err_ttl[index][0]
#
## Finding the brightest lc within the selected ones
#index = np.where(avg_mag_ttl == min(avg_mag_ttl[cond]))
#print '\nFrom those quasars, the brightest is', qso_name[index][0], 'with avg_mag',avg_mag_ttl[index][0] 
#print 'Its length is', lc_length[index][0] , 'and <avg_err_ttl>=', avg_err_ttl[index][0]


'''

Some plotting  : logs of ratios ...

'''

plt.plot(x,y,'.r')
nbins =70
H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure(figsize=[10,8])
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlim((xmin,xmax))
plt.ylim((ymin,ymax))
plt.axvline(0, color='r', lw=2)
plt.axhline(0,color='r',lw=2)
title = 'SDSS '+band[band_pos]+' band Jav vs Ch,  ' + str(percent)[:5]+'% objects of '\
        +str(len(ra_jav))+ ', prior='+prior  
plt.title(title)
plt.ylabel(r'$\log_{10}{ \, \left(  \tau_{jav} / \tau_{ch} \right)}$',fontsize=15)
plt.xlabel(r'$\log_{10}{ \, \left(  \sigma_{jav} / \sigma_{ch} \right)}$',fontsize=15)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname3='test_sdss_jav_ch_band_'+band[band_pos]+'_logs_ratios_prior-'+prior+'.png'
plt.savefig(fname3)
print 'Saving',  fname3

# Plot lc_length vs avg_mag... 
#plt.clf()
#title = 'SDSS quasars '+band[band_pos]+' selected '+ str(len(avg_mag_ttl[new]))
#plt.title(title)
#plt.xlabel('<magnitude>')
#plt.ylabel('Lightcurve length [days]')
#plt.scatter(avg_mag_ttl[new],lc_length[new],c=avg_err_ttl[new],cmap=plt.cm.coolwarm)
#fname1='test_sdss_'+band[band_pos]+'_band-SEL_avg_mag_vs_lc_length1.png'
#plt.savefig(fname1)


# plot lc_length  vs   avg_mag  as histogram 
#
#x=avg_mag_ttl[new]
#y=lc_length[new]
#plt.plot(x,y,'.r')
#nbins=100
#H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
#H = np.rot90(H)
#H = np.flipud(H)
#Hmasked = np.ma.masked_where(H==0,H)
#fig3 = plt.figure()
#plt.pcolormesh(xedges, yedges, Hmasked)
##plt.xlim((-3,1))
##plt.ylim((-2,1))
#title = 'SDSS quasars '+band[band_pos]+' selected '+ str(len(avg_mag_ttl[new]))
#plt.title(title)
#plt.xlabel('<magnitude>')
#plt.ylabel('Lightcurve length [days]')
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('Counts')
#fname3='test_sdss_'+band[band_pos]+'_band-sel_avg_mag_vs_lc_length.png'
#plt.savefig(fname3)
#
#'''
#Save all data 
#'''


#DAT = np.column_stack((qso_name,ra_jav,ra_ch,tau_jav,tau_ch,sigma_jav, sigma_ch,\
# lc_length, avg_mag_ttl, avg_err_ttl)) 
#
## sort the DAT column accoring to col(7) : log_tau_ratio : 
# 
##newDAT=DAT[DAT[:,6].argsort()]
#
#file_all = dir_out + 'qso_SDSS_jav_chelsea_and_stats_'+band[band_pos]+'_band_ALL.txt'
#np.savetxt(file_all,DAT,fmt="%s")
#print '\nSaving all the Chelsea vs Javelin SDSS comparison enriched with statistical',\
#' data (even those where log_tau_ratio would lead to infinity), for ',len(ra_jav) ,'quasars',\
#' for ', band[band_pos] ,' band to file ' , file_all 
#
#
#'''
#Save only selected data
#'''

#sel_DAT = np.column_stack((qso_name[new],ra_jav[new],ra_ch[new],tau_jav[new],tau_ch[new],\
#sigma_jav[new], sigma_ch[new], lc_length[new], avg_mag_ttl[new], avg_err_ttl[new]))
#
#file_sel = dir_out +  'qso_SDSS_jav_chelsea_compared_stats_'+band[band_pos]+'_band_SEL.txt'
#np.savetxt(file_sel,sel_DAT,fmt="%s")
#print '\nSaving all the Chelsea vs Javelin SDSS comparison enriched with statistical',\
#' data (even those where log_tau_ratio would lead to infinity), for ',len(qso_name[new]) ,\
#'quasars for ', band[band_pos] ,' band to file ' , file_sel 
#
#
#mag = avg_mag_ttl[new]
#avg_err_ttl[new]
#lc = lc_length[new]
#ind = np.where(mag < 18)
#np.where(lc[ind] > 80)
#
#lc[ind][41]
#
#indices = [6,41]
#for i in indices:
#    print qso_name[new][ind][i],lc[ind][i],mag[ind][i],err[ind][i],tau_jav[new][ind][i],\
#    tau_ch[new][ind][i],sigma_jav[new][ind][i],sigma_ch[new][ind][i]
