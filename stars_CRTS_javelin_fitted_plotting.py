# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 20:32:05 2014

@author: Chris 

Plot 2D histogram for  the CRTS stars fitted with  Javelin , output of 
condor_wrapper_CRTS_stars.py and 
javelin_chain_retrieve_stars.py


"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf
#import sys 

######################
#     lOAD DATA      # 
######################

dir_in_out = 'stars_CRTS_analysis/'

#args = sys.argv
#err = int(args[1])


chain_results = dir_in_out+ 'javelin_CRTS_stars_err_w_chain_results.txt'  

data = np.loadtxt(chain_results,dtype='str' )

fname = np.empty(0,dtype=str)
sigma_m = np.empty(0,dtype=float)
tau_m =  np.empty(0,dtype=float)


print 'Reading in all the values ... '

for i in range(len(data[:,0])):
   try:
       fname = np.append(fname, data[i,0])
       sigma_m = np.append(sigma_m, float(data[i,2]))
       tau_m = np.append(tau_m, float(data[i,5]))
      
   except ValueError:     
       pass

if len(sigma_m) != len(tau_m) : 
    m = min(len(tau_m), len(sigma_m))
    sigma_m = sigma_m[0:m]
    tau_m = tau_m[0:m]
    fname = fname[0:m]

assert len(sigma_m) == len(tau_m)
 
print '\nOut of ', len(data[:,0]), ' rows we were able to read in ', len(tau_m)

# RECALCULATE SIGMA_HAT FROM SIGMA OF JAVELIN

sigma_hat = sigma_m * np.sqrt(tau_m / (2.0 * 365.0))

############################
# SELECTING POINTS TO USE  #
############################

def sel_points(dir_in_out, fname):
    good_LC = np.loadtxt(dir_in_out + 'good_err_LC.txt', dtype='str')
    good_LC_cut = np.empty(0, dtype=str)

    for i in range(len(good_LC)):
        good_LC_cut = np.append(good_LC_cut, good_LC[i][4:-8])
        
    good_LC_mask = np.zeros_like(fname, dtype='bool')
    for i in range(len(fname)):
        print '\nComparison in progress...', str((float(i) / float(len(fname)) )*100.0)[:5], '%'
        good_LC_mask[i] =  fname[i][4:] in  good_LC_cut 
        
        
    print 'Out of ', len(fname), 'objects, we use ',  good_LC_mask.sum()
    return good_LC_mask
    
good_LC_mask = sel_points(dir_in_out, fname)

###############################
# QUICK 1D HISTOGRAM OF SIGMA #
###############################

def histogram1D_sigma(sigma_m, dir_in_out ):
    font = 20 
    lsm = np.log10(sigma_m)
    
    xmin = np.percentile(lsm,5)
    xmax = np.percentile(lsm,95)
    
    plt.clf()
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)
    ax1 = fig1.add_subplot(gs[:,:90])   
    
    ax1.hist(lsm,100,range=[xmin,xmax],facecolor='gray',align='mid',alpha=0.5)
    
    plt.title('CRTS Stars, err_w ',fontsize=font)
    plt.xlabel(r'$\log_{10}{ \, \left(  \sigma_{jav} \right)}$',fontsize=font)
    plt.ylabel('Number of counts',fontsize=font )
    
    ax1.tick_params(axis='x', labelsize=font)
    ax1.tick_params(axis='y', labelsize=font)
    fname = dir_in_out+'CRTS_stars_sigma_histogram.png'
    plt.savefig(fname)
    print 'We saved ', fname 
    
histogram1D_sigma(sigma_m, dir_in_out)

######  functions below based on qso_CRTS_SDSS_matched_plotting.py

#########################
#      SET LIMITS       #
#########################

xmin = np.log10(0.1)   # sigma limits 
xmax = np.log10(5.0)
ymin = np.log10(1.0)   # tau limits 
ymax = np.log10(10000.0)

def set_limits(x_values,y_values, xmin, xmax, ymin, ymax):
    # x_values and y_values are sigma and tau 
    # that are not constrained by anything
    # this function checks whether the constraint 
    # below is taking away more than 95% of points... 
    x0,x1, y0, y1 = xmin, xmax, ymin, ymax 
    x = x_values
    y = y_values
    
    if xmax < np.percentile(x,94):
        xmax = np.percentile(x,93)
        print 'We had to change the x_max from', x1, ' to ', xmax
    if ymax < np.percentile(y,90):
        ymax = np.percentile(y, 93)
        print 'We had to change the y_max from', y1, ' to ', ymax
    if xmin > np.percentile(x, 7):
        xmin = np.percentile(x,5)
        print 'We had to change the x_min from', x0, ' to ', xmin
    if ymin > np.percentile(y, 7):
        ymin = np.percentile(y,5)
        print 'We had to change the x_min from', y0, ' to ', ymin
    
    x_lim = [xmin, xmax]
    y_lim = [ymin,ymax]
    
    return  x_lim, y_lim 
    
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
    ' < y < ',y_limits[1], ' and have no infinities,  and only those are used for plotting ...  '
    
    return x[gi], y[gi], non_inf, percent


def histogram2D(x_arr, y_arr, number, percent, xlim, ylim, dir_out):
    # args could include javelin results_file , from which you can 
    # take the info about the prior  
    font=20
    x,y = x_arr,y_arr
    nbins =50

    # Define the canvas to work on   and the  grid  
    plt.clf()
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)    
            
    # First histogram  : Chelsea results 
    H, xedges,yedges = np.histogram2d(x,y,bins=nbins)  
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    ax1 = fig1.add_subplot(gs[:,:90])   
    pcObject1 = ax1.pcolormesh(xedges, yedges, Hmasked)
    
    xmin = xlim[0]
    xmax = xlim[1]
    ymin =  ylim[0]
    ymax =  ylim[1]
    
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    
    x_label_jav = r'$\log_{10}{ \, \left(  \hat\sigma_{jav} \right)}$'
    y_label_jav = r'$\log_{10}{ \, \left(  \tau_{jav} \right)}$'
    
    plt.ylabel(y_label_jav,fontsize=font)
    plt.xlabel(x_label_jav,fontsize=font)
    title_hist = 'CRTS stars Javelin results, '+ str(number)[:5] + ', i.e.  ' + str(percent)[:5]+ '% points'
    
    plt.title(title_hist, fontsize=font)
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:,95:])
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
    cbar.ax.set_ylabel('Counts', fontsize=font)
    axC.tick_params(axis='y', labelsize=font)
    
    # save file
    fname = dir_out + 'CRTS_stars_Javelin_sigma_hat_tau.png'
    plt.savefig(fname)
    print 'File saved is ', fname


sigma_hat = sigma_hat[good_LC_mask]
tau_m = tau_m[good_LC_mask]

xlim, ylim =  set_limits(np.log10(sigma_hat), np.log10(tau_m), xmin, xmax, ymin, ymax)
x_arr, y_arr, number, percent = load_x_y(np.log10(sigma_hat), np.log10(tau_m), xlim, ylim)
histogram2D(x_arr, y_arr, number, percent, xlim, ylim, dir_in_out)
