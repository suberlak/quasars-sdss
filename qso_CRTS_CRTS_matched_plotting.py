# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 13:58:42 2014

@author: astronomy

modified  qso_CRTS_SDSS_matched_plotting...

match the Chelsea CRTS results to   JAVELIN CRTS Results that were already 
matched to Chelsea SDSS S82 results. 

meant to plot the log(sigma_hat) vs log(tau)  , as well as log(sigma_hat) 
vs log(sigma_hat)

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf

######################
#     lOAD DATA      # 
######################

dir_in = 'QSO_CRTS_analysis/'  
dir_out = 'QSO_CRTS_analysis/'


results_file =  'javelin_CRTS_err_w_Chelsea_s82drw_r_compare.txt'
output = dir_in+ results_file

data = np.loadtxt(output,dtype='str' )

qso_name = data[:,0]       
ra_crts = data[:,1].astype(float)  # ra and dec in degrees 
ra_sdss = data[:,2].astype(float)  
dec_crts = data[:,3].astype(float) 
dec_sdss = data[:,4].astype(float) 
tau_med_jav_crts = data[:,5].astype(float) 
tau_chelsea_sdss = data[:,6].astype(float)  
sigma_med_jav_crts = data[:,7].astype(float)  
sigma_hat_jav_crts = sigma_med_jav_crts * np.sqrt(tau_med_jav_crts / (2.0*365))
sigma_chelsea_sdss = data[:,8].astype(float)   
sigma_hat_chelsea_sdss = sigma_chelsea_sdss * np.sqrt(tau_chelsea_sdss /(2.0*365) )
xmin = 0.001   # sigma limits 
xmax = 5
ymin = 1   # tau limits 
ymax = 70000

xlim = [xmin, xmax]
ylim = [ymin, ymax]


##########################
# SELECT POINTS TO USE   #
##########################
# Those sent to Chelsea were already from good_err_LC.txt list, so now only need 

good_LC = np.loadtxt(dir_in + 'good_err_LC.txt', dtype='str')
good_LC_cut = np.empty(0, dtype=str)

for i in range(len(good_LC)):
    good_LC_cut = np.append(good_LC_cut, good_LC[i][4:-4])

good_LC_mask = np.zeros_like(qso_name, dtype='bool')
for i in range(len(qso_name)):
    print '\nComparison in progress...', str((float(i) / float(len(qso_name)) )*100.0)[:5], '%'
    good_LC_mask[i] =  qso_name[i] in good_LC_cut     
  
        
print 'Out of ', len(qso_name), 'objects, we use ',  good_LC_mask.sum()

###########################
# DEFINE NEEDED FUNCTION  #
###########################

def load_x_y(x_arr, y_arr, x_limits, y_limits):

    print '\n Loading x and y ... ' 
   
    x = x_arr 
    y = y_arr
    
    # sieve out suspiciously bad values , based only on x and y  
    xinf = np.asarray(map(isinf,x),dtype=bool)
    yinf = np.asarray(map(isinf,y),dtype=bool)
    ttlinf = xinf + yinf
    # ttlwh = np.where(ttlinf == True)  list of good indices
    gi = -ttlinf  # good_indices 
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


def histogram2D(x_arr, y_arr, number, percent, xlim, ylim, title, dir_out):
    # args could include javelin results_file , from which you can 
    # take the info about the prior  
    font = 20
    
    x = np.log10(x_arr)
    y = np.log10(y_arr)
    
    nbins =50
    plt.clf()
    fig1 = plt.figure()
        
    # Define the canvas to work on   and the  grid  
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)    
            
    # First histogram  : Chelsea results 
    H, xedges,yedges = np.histogram2d(x,y,bins=nbins)  
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    ax1 = fig1.add_subplot(gs[:,:90])   
    pcObject1 = ax1.pcolormesh(xedges, yedges, Hmasked)
    ax1.tick_params(axis='x', labelsize=font)
    ax1.tick_params(axis='y', labelsize=font)

    xmin = np.log10(xlim[0])
    xmax = np.log10(xlim[1])
    ymin =  np.log10(ylim[0])
    ymax =  np.log10(ylim[1])
    
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    
    
    x_label_ch = r'$\log_{10}{ \, \left(  \hat\sigma_{ch} \right)}$'
    y_label_ch = r'$\log_{10}{ \, \left( \tau_{ch} \right)}$'
    x_label_jav = r'$\log_{10}{ \, \left(  \hat\sigma_{jav} \right)}$'
    y_label_jav = r'$\log_{10}{ \, \left(  \tau_{jav} \right)}$'
     
    
    if title == 'ch' : 
        plt.ylabel(y_label_ch,fontsize=font+5)
        plt.xlabel(x_label_ch,fontsize=font+5)
        title_hist = 'S82 SDSS Chelsea results, '+ str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'Chelsea_s82_SDSS_matched_sigma_hat_tau.png' 
    if title == 'jav' :
        plt.ylabel(y_label_jav,fontsize=font+5)
        plt.xlabel(x_label_jav,fontsize=font+5)
   
        title_hist = 'CRTS Javelin results, '+ str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_Javelin_matched_sigma_hat_tau.png'
    if title == 'ss' : 
        plt.xlabel(x_label_ch,fontsize=font+5)
        plt.ylabel(x_label_jav,fontsize=font+5)
        title_hist = 'S82 SDSS Chelsea vs CRTS JAVELIN, '+str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_SDSS_matched_sigma_hat_sigma_hat.png'
    if title == 'tt' :
        plt.xlabel(y_label_ch, fontsize=font+5)
        plt.ylabel(y_label_jav,fontsize=font+5)
        title_hist = 'S82 SDSS Chelsea vs CRTS JAVELIN, '+str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_SDSS_matched_tau_tau.png'  
        
    plt.title(title_hist, fontsize = font)
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:,95:])
    axC.tick_params(axis='y', labelsize=font)
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
    cbar.ax.set_ylabel('Counts', fontsize=font)
    
    plt.savefig(fname)
    print 'File saved is ', fname

# Make log(sigma_hat)  vs log(tau)  histogram for Chelsea 
#x_arr, y_arr, number, percent = load_x_y(sigma_hat_chelsea_sdss[good_LC_mask], tau_chelsea_sdss[good_LC_mask], xlim, ylim)
#histogram2D(x_arr, y_arr, number, percent, xlim, ylim, 'ch', dir_out)
##
### Make log(sigma_hat)  vs log(tau)  histogram for Javelin CRTS  
#x_arr, y_arr, number, percent = load_x_y(sigma_hat_jav_crts, tau_med_jav_crts, xlim, ylim)
#histogram2D(x_arr, y_arr, number, percent, xlim, ylim, 'jav', dir_out)
##
### Make log(sigma_hat) vs log(sigma_hat) histogram  
##
x_arr, y_arr, number, percent = load_x_y(sigma_hat_chelsea_sdss[good_LC_mask], sigma_med_jav_crts[good_LC_mask], xlim, xlim)
#histogram2D(x_arr, y_arr, number, percent, xlim, xlim, 'ss', dir_out)
##
### Make log(tau) vs log(tau) histogram 
##
#x_arr, y_arr, number, percent = load_x_y(tau_chelsea_sdss[good_LC_mask], tau_med_jav_crts[good_LC_mask], ylim, ylim)
#histogram2D(x_arr, y_arr, number, percent, ylim, ylim, 'tt', dir_out)

def median_and_rms(array):
    median = np.median(array)
    rms = np.percentile(array, 75) - np.percentile(array,25)
    print median, rms
    
print '\nMedian and RMs from IQR (75-25%) range for Chelsea sigma hat is ', median_and_rms(x_arr)

print '\nMedian and RMs from IQR (75-25%) range for Javelin CRTS sigma hat is ', median_and_rms(y_arr)