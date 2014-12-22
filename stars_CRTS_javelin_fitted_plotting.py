# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 20:32:05 2014

@author: astronomy

Plot 2D histogram for  the CRTS stars fitted with  Javelin , output of 
condor_wrapper_CRTS_stars.py and 
javelin_chain_retrieve_stars.py


"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf
import sys 

######################
#     lOAD DATA      # 
######################

dir_in = ['stars_CRTS_err_rms_chains/', 'stars_CRTS_err_w_chains/']
dir_out = 'stars_CRTS_analysis/'

args = sys.argv
err = float(args[1])

if err ==0 : 
    err_txt = 'err_rms'
else: 
    err_txt = 'err_w'
    
chain_results = dir_in[err]+ 'javelin_CRTS_stars_'+err_txt+'_chain_results.txt'  

data = np.loadtxt(chain_results,dtype='str' )

# fname =  data[:,0]
#sigma_l
sigma_m = data[:,2]
# sigma_h
#tau_l 
tau_m = data[:,5]
# tau_h

######  functions below based on qso_CRTS_SDSS_matched_plotting.py

#########################
#      SET LIMITS       #
#########################

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


def histogram(x_arr, y_arr, number, percent, xlim, ylim, title, dir_out,err_choice):
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
        
        title_hist = 'CRTS stars Javelin results,'+err+', '+ str(number)[:5] + ', i.e.  ' + str(percent)[:5]+ '% points'
        fname = dir_out + 'CRTS_stars_Javelin_'+err+'.png'
        
    
    plt.title(title_hist)
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:,95:])
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
    cbar.ax.set_ylabel('Counts')
    
    plt.savefig(fname)
    print 'File saved is ', fname

# Make log(sigma)  vs log(tau)  histogram for Chelsea 
x_arr, y_arr, number, percent = load_x_y(sigma_m, tau_m, xlim, ylim)
histogram(x_arr, y_arr, number, percent, xlim, ylim, 'jav', dir_out,err)
