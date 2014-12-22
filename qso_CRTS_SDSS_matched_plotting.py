# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 13:58:42 2014

@author: astronomy

modified  qso_DRW_plotting...

meant to plot the log(sigma) vs log(tau)  for Chelsea S82 results 
matched to Javelin fits of CRTS data,  from  chelsea_results_load_crts.py

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf

######################
#     lOAD DATA      # 
######################

dir1 = 's82drw/'
dir_out = 'QSO_CRTS_analysis/'
err_choice = 1 #  (or 1)

out_names = ['javelin_CRTS_err_rms_Chelsea_s82drw_r_compare.txt',  'javelin_CRTS_err_w_Chelsea_s82drw_r_compare.txt']
output = dir1+ out_names[err_choice]

data = np.loadtxt(output,dtype='str' )

qso_name = data[:,0]       
ra_crts = data[:,1].astype(float)  # ra and dec in degrees 
ra_sdss = data[:,2].astype(float)  
dec_crts = data[:,3].astype(float) 
dec_sdss = data[:,4].astype(float) 
tau_med_jav_crts = data[:,5].astype(float) 
tau_chelsea_sdss = data[:,6].astype(float)  
sigma_med_jav_crts = data[:,7].astype(float)  
sigma_chelsea_sdss = data[:,8].astype(float)   

xmin = 0.001   # sigma limits 
xmax = 5
ymin = 1   # tau limits 
ymax = 70000

xlim = [xmin, xmax]
ylim = [ymin, ymax]

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


def histogram(x_arr, y_arr, number, percent, xlim, ylim, title, dir_out):
    # args could include javelin results_file , from which you can 
    # take the info about the prior  

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
    
    xmin = np.log10(xlim[0])
    xmax = np.log10(xlim[1])
    ymin =  np.log10(ylim[0])
    ymax =  np.log10(ylim[1])
    
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    
    
    x_label_ch = r'$\log_{10}{ \, \left(  \sigma_{ch} \right)}$'
    y_label_ch = r'$\log_{10}{ \, \left( \tau_{ch} \right)}$'
    x_label_jav = r'$\log_{10}{ \, \left(  \sigma_{jav} \right)}$'
    y_label_jav = r'$\log_{10}{ \, \left(  \tau_{jav} \right)}$'
     
    
    if title == 'ch' : 
        plt.ylabel(y_label_ch,fontsize=15)
        plt.xlabel(x_label_ch,fontsize=15)
        title_hist = 'S82 SDSS Chelsea results, '+ str(number) + ', i.e.  ' + str(percent)+ '% points'
        fname = dir_out + 'Chelsea_s82_SDSS_matched.png' 
    else:
        plt.ylabel(y_label_jav,fontsize=15)
        plt.xlabel(x_label_jav,fontsize=15)
        if err_choice == 0 : 
            err = 'err_rms'
        else: 
            err = 'err_w'
        
        title_hist = 'CRTS Javelin results,'+err+', '+ str(number) + ', i.e.  ' + str(percent)+ '% points'
        fname = dir_out + 'CRTS_Javelin_matched_'+err+'.png'
        
    
    plt.title(title_hist)
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:,95:])
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
    cbar.ax.set_ylabel('Counts')
    
    plt.savefig(fname)
    print 'File saved is ', fname

# Make log(sigma)  vs log(tau)  histogram for Chelsea 
x_arr, y_arr, number, percent = load_x_y(sigma_chelsea_sdss, tau_chelsea_sdss, xlim, ylim)
histogram(x_arr, y_arr, number, percent, xlim, ylim, 'ch', dir_out)

# Make log(sigma)  vs log(tau)  histogram for Javelin CRTS  
x_arr, y_arr, number, percent = load_x_y(sigma_med_jav_crts, tau_med_jav_crts, xlim, ylim)
histogram(x_arr, y_arr, number, percent, xlim, ylim, 'jav', dir_out)
