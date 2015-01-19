# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 13:45:05 2014

@author: astronomy

modified  qso_CRTS_CRTS_matched_plotting...

match the Chelsea CRTS results to   JAVELIN CRTS Results that were already 
matched to Chelsea SDSS S82 results. 

Plots three figures for the poster:

1) Chelsea for  log(tau_CRTS) vs  log(tau_SDSS)
2) Chelsea for log(sigma_hat_SDSS) vs log(sigma_hat_CRTS)
3) Chelsea for log(SF_inf_SDSS) vs log(SF_inf_CRTS)

Using ONLY the data that fulfils the CRTS selection criteria (Pnoise-Plike) , 
and err_flag ... 

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf

##################################################################
#           lOAD :     JAV CRTS - Chelsea SDSS matched data      # 
##################################################################

dir_in = 'QSO_CRTS_analysis/'  
dir_out = 'QSO_CRTS_analysis/'

matched_data =  'javelin_CRTS_err_w_Chelsea_s82drw_r_compare.txt'
output = dir_in+ matched_data


data = np.loadtxt(output,dtype='str' )

qso_name = data[:,0]       
ra_crts = data[:,1].astype(float)  # ra and dec in degrees 
ra_sdss = data[:,2].astype(float)  
dec_crts = data[:,3].astype(float) 
dec_sdss = data[:,4].astype(float) 
tau_med_jav_crts = data[:,5].astype(float) 
tau_ch_sdss = data[:,6].astype(float)  
sigma_med_jav_crts = data[:,7].astype(float)  
sigma_hat_jav_crts = sigma_med_jav_crts * np.sqrt( (2.0*365) / tau_med_jav_crts  )

sigma_ch_sdss = data[:,8].astype(float)   
sigma_hat_ch_sdss = sigma_ch_sdss * np.sqrt((2.0*365) / tau_ch_sdss  )
sfinf_ch_sdss = sigma_hat_ch_sdss * np.sqrt(tau_ch_sdss / 365.0)
#############################
#   LOAD Chelsea CRTS fits  #
#############################
# REMEMBER THAT HER CODE RETURNS SIGMA HAT !!! 
chelsea_crts = 'Chelsea_CRTS_fits.dat'
out = dir_in + chelsea_crts

data1 = np.loadtxt(out, dtype='str')
qso_name1 = data1[:,0]
log_10_tau_ch_crts = data1[:,2].astype(float)
log_10_sigma_hat_ch_crts = data1[:,3].astype(float)

tau_ch_crts = np.power(10,log_10_tau_ch_crts)
sigma_hat_ch_crts = np.power(10,log_10_sigma_hat_ch_crts)

sigma_ch_crts = sigma_hat_ch_crts * np.sqrt(tau_ch_crts /(2.0*365) )
sfinf_ch_crts = sigma_hat_ch_crts * np.sqrt(tau_ch_crts / 365.0)

edge_flag = data1[:,8].astype(float)
Plike = data1[:,9].astype(float)
Pnoise = data1[:,10].astype(float)


# Those sent to Chelsea were already error-selected and length - selected 
# (longer than 10 lines) , so at least edge should be !=0 for all...

#############################
# IMPOSE SELECTION CRITERIA #
#############################

cond1 = (Plike - Pnoise)> 2.0

cond2 = edge_flag == 0.0

cond = cond1 & cond2
qso_name1 = qso_name1[cond]

print 'Our conditions reduce the Chelsea CRTS fits data from', len(cond), ' to ', len(qso_name1) 

#################################################
#  MATCH CHELSEA CRTS TO JAV CRTS-CHELSEA SDSS  #
#################################################

good_ch_CRTS_mask = np.zeros_like(qso_name1, dtype=bool)
good_jav_ch_mask = np.zeros_like(qso_name, dtype = bool )

for i in range(len(qso_name1)): 
    good_ch_CRTS_mask[i] = qso_name1[i][4:-4] in qso_name
    
for i in range(len(qso_name)):
    good_jav_ch_mask[i] = 'out_'+qso_name[i]+'.txt' in  qso_name1
    
print 'For Chelsea CRTS out of ', len(qso_name1), ' we have ',\
 good_ch_CRTS_mask.sum() , 'of those LC that were matched with SDSS'

print 'For Jav Chelsea CRTS-SDSS matched, out of', len(qso_name), 'we have',\
 good_jav_ch_mask.sum() , 'of those LCs that are also present in Chelseas CRTS fits' 


###########################
# DEFINE NEEDED FUNCTION  #
###########################

def load_x_y(x_arr, y_arr, x_limits, y_limits):
    print '\n ---------------------------------------'
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
    
    # remove NaNs
    y_nan = np.isnan(x)
    x_nan = np.isnan(y)
    gi[y_nan] = False
    gi[x_nan] = False 
    
    print 'We have ', y_nan.sum(), 'NaNs in y'
    print 'and  ', x_nan.sum(), 'NaNs in x'
    
    
    good_len = len(np.where(gi == True)[0])
    
    percent = (float(good_len) / float(len(x)))  * 100.0
   
    print 'Out of ', len(x),' rows, we have ', good_len, ' of those that match', \
    'the criteria of ',  x_limits[0],' < x <', x_limits[1],' and ', y_limits[0],\
    ' < y < ',y_limits[1], 'and only those are used for plotting ...  '
    
    return x[gi], y[gi], good_len, percent


def histogram2D(x_arr, y_arr, number, percent, xlim, ylim, title, dir_out):
    # args could include javelin results_file , from which you can 
    # take the info about the prior  
    font = 25
    
    x = np.log10(x_arr)
    y = np.log10(y_arr)
    
    nbins =30
    plt.clf()
    fig1 = plt.figure()
        
    # Define the canvas to work on   and the  grid  
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)    
            
    # First histogram  : Chelsea results 
    # return x, y
    
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
    
    x = np.linspace(xmin,xmax)
    y = x
    plt.plot(x,y, '--',color='r',lw=8)     
    
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    
    title_hist = 'CRTS data vs SDSS data' 
    # Three titles that I need ... 
    
    inc = 8
   
    if title == 'ss' : 
        plt.xlabel(r'$\log_{10}{ \, \left(  \hat\sigma_{CRTS} \right)}$',fontsize=font+inc)
        plt.ylabel(r'$\log_{10}{ \, \left(  \hat\sigma_{SDSS} \right)}$',fontsize=font+inc)
        fname = dir_out + 'poster_CRTS_ch_SDSS_ch_match_sigma_hat.png'
   
    if title == 'sfinf' :
        plt.xlabel(r'$\log_{10}{ \, \left(  SF_{\infty, CRTS} \right)}$', fontsize=font+inc)
        plt.ylabel(r'$\log_{10}{ \, \left(  SF_{\infty, SDSS} \right)}$',fontsize=font+inc)
        fname = dir_out + 'poster_CRTS_ch_SDSS_ch_match_sf_inf.png'  
 
    if title == 'tt' :
        plt.xlabel(r'$\log_{10}{ \, \left( \tau_{CRTS} \right)}$', fontsize=font+inc)
        plt.ylabel(r'$\log_{10}{ \, \left( \tau_{SDSS} \right)}$',fontsize=font+inc)
        fname = dir_out + 'poster_CRTS_ch_SDSS_ch_match_tau_tau.png'     
        
    plt.title(title_hist, fontsize = font)
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:,95:])
    axC.tick_params(axis='y', labelsize=font)
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
    cbar.ax.set_ylabel('Counts', fontsize=font)
    
    plt.savefig(fname)
    print 'File saved is ', fname


def median_and_rms(array, name):
    median = np.median(array)
    rms = np.percentile(array, 75) - np.percentile(array,25)
    print '\nMedian and RMs from IQR (75-25%) range for', name, median, rms
    

##
### Make log(sigma_hat) vs log(sigma_hat) histogram   Chelsea CRTS : Chelsea SDSS
##
xlim = [np.power(10,-2.0),np.power(10,1.4)]
ylim= [np.power(10,-2.0), np.power(10,0.5)]
x_arr, y_arr, number, percent = load_x_y(sigma_hat_ch_crts[good_ch_CRTS_mask], sigma_hat_ch_sdss[good_jav_ch_mask], xlim, ylim)
histogram2D(x_arr, y_arr, number, percent, xlim, ylim, 'ss', dir_out)
median_and_rms(x_arr, 'sigma hat Chelsea CRTS')
median_and_rms(y_arr, 'sigma hat Chelsea SDSS' )

###
#### Make log(tau) vs log(tau) histogram Chelsea CRTS : Chelsea SDSS
###
xlim=[np.power(10,-0.5), np.power(10,4)]
ylim=[np.power(10,1.0), np.power(10,4.5)]
x_arr, y_arr, number, percent = load_x_y(tau_ch_crts[good_ch_CRTS_mask], tau_ch_sdss[good_jav_ch_mask], xlim, ylim)
histogram2D(x_arr, y_arr, number, percent, xlim, ylim, 'tt', dir_out)
median_and_rms(x_arr, 'tau Chelsea CRTS')
median_and_rms(y_arr, 'tau Chelsea SDSS')

###
#### Make log(sf_inf) vs log(sf_inf) histogram  Chelsea CRTS : Chelsea SDSS
###
xlim = [np.power(10,-1.5),5]
ylim= [np.power(10,-1.5), 5]
x_arr, y_arr, number, percent = load_x_y(sfinf_ch_crts[good_ch_CRTS_mask], sfinf_ch_sdss[good_jav_ch_mask], xlim,ylim)
histogram2D(x_arr, y_arr, number, percent,xlim,ylim, 'sfinf', dir_out)
median_and_rms(x_arr, 'sfinf Chelsea CRTS')
median_and_rms(y_arr, 'sfinf Chelsea SDSS')
