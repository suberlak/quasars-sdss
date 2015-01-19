# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 13:45:05 2014

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
sigma_hat_jav_crts = sigma_med_jav_crts * np.sqrt(tau_med_jav_crts / (2.0*365))
sigma_chelsea_sdss = data[:,8].astype(float)   
sigma_hat_chelsea_sdss = sigma_chelsea_sdss * np.sqrt(tau_ch_sdss /(2.0*365) )

#############################
#   LOAD Chelsea CRTS fits  #
#############################

chelsea_crts = 'Chelsea_CRTS_fits.dat'
out = dir_in + chelsea_crts

data1 = np.loadtxt(out, dtype='str')
qso_name1 = data1[:,0]
log_10_tau_chelsea_crts = data1[:,2].astype(float)
log_10_sigma_chelsea_crts = data1[:,3].astype(float)

tau_ch_crts = np.power(10,log_10_tau_chelsea_crts)
sigma_ch_crts = np.power(10,log_10_sigma_chelsea_crts)

sigma_hat_chelsea_crts = sigma_ch_crts * np.sqrt(tau_ch_crts /(2.0*365) )

edge_flag = data1[:,8].astype(float)
Plike = data1[:,9].astype(float)
Pnoise = data1[:,10].astype(float)


# Those sent to Chelsea were already error-selected and length - selected 
# (longer than 10 lines) , so at least edge should be !=0 for all...

######################
#     SET LIMITS     # 
######################

xmin = 0.001   # sigma limits 
xmax = 5
ymin = 1   # tau limits 
ymax = 70000

xlim = [xmin, xmax]
ylim = [ymin, ymax]


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
    
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    
    
    x_label_ch = r'$\log_{10}{ \, \left(  \hat\sigma_{ch} \right)}$'
    y_label_ch = r'$\log_{10}{ \, \left( \tau_{ch} \right)}$'
    x_label_jav = r'$\log_{10}{ \, \left(  \hat\sigma_{jav} \right)}$'
    y_label_jav = r'$\log_{10}{ \, \left(  \tau_{jav} \right)}$'
     
    
    if title == 'ch' : 
        plt.ylabel(y_label_ch,fontsize=font+5)
        plt.xlabel(x_label_ch,fontsize=font+5)
        title_hist = 'S82 CRTS Chelsea results, '+ str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_ch_matched_sigma_hat_tau.png' 
    if title == 'jav' :
        plt.ylabel(y_label_jav,fontsize=font+5)
        plt.xlabel(x_label_jav,fontsize=font+5)
        title_hist = 'CRTS Javelin results, '+ str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_Javelin_matched_sigma_hat_tau.png'
    if title == 'ss1' : 
        plt.xlabel(x_label_ch,fontsize=font+5)
        plt.ylabel(x_label_ch,fontsize=font+5)
        title_hist = 'S82 CRTS Chelsea vs SDSS Chelsea , '+str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_ch_SDSS_ch_matched_sigma_hat_sigma_hat.png'
    if title == 'ss2' : 
        plt.xlabel(x_label_ch,fontsize=font+5)
        plt.ylabel(x_label_jav,fontsize=font+5)
        title_hist = 'S82 CRTS Chelsea vs CRTS Javelin  , '+str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_ch_CRTS_jav_matched_sigma_hat_sigma_hat.png'   
    if title == 'tt1' :
        plt.xlabel(y_label_ch, fontsize=font+5)
        plt.ylabel(y_label_jav,fontsize=font+5)
        title_hist = 'S82 CRTS Chelsea vs SDSS Chelsea, '+str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_ch_SDSS_ch_matched_tau_tau.png'  
    if title == 'tt2' :
        plt.xlabel(y_label_ch, fontsize=font+5)
        plt.ylabel(y_label_jav,fontsize=font+5)
        title_hist = 'S82 CRTS Chelsea vs CRTS Javelin , '+str(number) + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_ch_CRTS_jav_matched_tau_tau.png'     
        
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
    
# Make log(sigma_hat)  vs log(tau)  histogram for Chelsea  CRTS
x_arr, y_arr, number, percent = load_x_y(sigma_hat_chelsea_crts[good_ch_CRTS_mask], tau_ch_crts[good_ch_CRTS_mask], xlim, ylim)
histogram2D(x_arr, y_arr, number, percent, xlim, ylim, 'ch', dir_out)
median_and_rms(x_arr, 'sigma hat Chelsea CRTS')
#median_and_rms(y_arr, 'tau Chelsea CRTS' )
##
##
### Make log(sigma_hat) vs log(sigma_hat) histogram   Chelsea CRTS : Chelsea SDSS
##
x_arr, y_arr, number, percent = load_x_y(sigma_hat_chelsea_crts[good_ch_CRTS_mask], sigma_hat_chelsea_sdss[good_jav_ch_mask], xlim, xlim)
histogram2D(x_arr, y_arr, number, percent, xlim, xlim, 'ss1', dir_out)
median_and_rms(x_arr, 'sigma hat Chelsea CRTS')
median_and_rms(y_arr, 'sigma hat Chelsea SDSS' )

### Make log(sigma_hat) vs log(sigma_hat) histogram   Chelsea CRTS : Jav CRTS
##
x_arr, y_arr, number, percent = load_x_y(sigma_hat_chelsea_crts[good_ch_CRTS_mask], sigma_hat_jav_crts[good_jav_ch_mask], xlim, xlim)
histogram2D(x_arr, y_arr, number, percent, xlim, xlim, 'ss2', dir_out)
#median_and_rms(x_arr, 'sigma hat Chelsea CRTS')
median_and_rms(y_arr, 'sigma hat Javelin CRTS' )

###
#### Make log(tau) vs log(tau) histogram 
###
x_arr, y_arr, number, percent = load_x_y(tau_ch_crts[good_ch_CRTS_mask], tau_ch_sdss[good_jav_ch_mask], ylim, ylim)
histogram2D(x_arr, y_arr, number, percent, ylim, ylim, 'tt1', dir_out)
#median_and_rms(x_arr, 'tau Chelsea CRTS')
#median_and_rms(y_arr, 'tau Chelsea SDSS')

###
#### Make log(tau) vs log(tau) histogram 
###
x_arr, y_arr, number, percent = load_x_y(tau_ch_crts[good_ch_CRTS_mask], tau_med_jav_crts[good_jav_ch_mask], ylim, ylim)
histogram2D(x_arr, y_arr, number, percent, ylim, ylim, 'tt2', dir_out)
median_and_rms(x_arr, 'tau Chelsea CRTS')
median_and_rms(y_arr, 'tau Javelin CRTS')
