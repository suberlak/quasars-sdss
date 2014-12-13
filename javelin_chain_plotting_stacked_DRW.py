# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 14:48:12 2014

@author: suberlak

Plotting stacked chain posterior distributions to investigate their gaussianity


Chains are named  :
 ch_DRWtest_LC547_err2.dat_chain.dat 

"""

import numpy as np 
import matplotlib.pyplot as plt 


dir_in = ['qso_drw_analysis/','qso_drw_medium_analysis/']

dir_out= ['qso_drw_analysis/','qso_drw_medium_analysis/']


dir_an = dir_in[1]  # CHOOSE WHICH ANALYSIS DIR WE NEED 
length= 'medium'
# ARRAYS WITH ITERATED PARAMETERS 
error = [1,2]  # need to add 3 if also considering error 3
with_prior = ['no']  # need to add 'yes' if we have both cases 

# DEFINE A FUNCTION TO PRINT STACKED CHAINS FOR A
# GIVEN  ERROR AND PRIOR SETTING 

def printing(with_prior, chain_stack):
    global dir_an

    if with_prior == 'yes' : 
        filename = dir_an + chain_stack+'_with_prior.npy'
        files=np.load(filename)
        post = '_with_prior'
    else : 
        filename = dir_an + chain_stack+'_no_prior.npy'
        files=np.load(filename)
        post = '_no_prior'    
    
    err  = chain_stack[26]
    
    
    #########################
    #     LOADING CHAINS    #
    #########################
    
    log10_sigma = files[0]
    log10_tau = files[1]
    
    # REMOVING THE OUTLIERS 
    gi = np.ones(len(log10_sigma), dtype=bool)
    xsmall = np.where(log10_sigma < -1.5)
    xbig   = np.where(log10_sigma > 1.0)
    ysmall = np.where(log10_tau < 1 )
    ybig   = np.where(log10_tau > 6 )
    gi[xsmall] = False
    gi[xbig]   = False
    gi[ysmall] = False
    gi[ybig]   = False
    
    # ARRAYS OF FILTERED DATA
    x = log10_sigma[gi]
    y = log10_tau[gi]
    
    print 'Out of ', len(log10_sigma), ' we used ', len(x) , ' rows which satisfy',\
           ' selection criteria ' 
    
    # PLOTTING THE FILTERED DATA 
    plt.plot(x,y,'.r')
    nbins =100
    
    H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
    
    # finding the maximum of the 2D distribution    
    # a,b = np.where(H == H.max())    
    #x_max = xedges[a[0]]
    #y_max = yedges[b[0]]
    
    y_true = np.log10(100)
    x_true = np.log10(0.2)
    
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    plt.figure()
    plt.pcolormesh(xedges, yedges, Hmasked)
       
    # PRINTING PART  
    title = 'Stacked chains '+length+' prior=' + with_prior +', err '+ err
    plt.title(title)
    plt.axhline(y_true,color='r',lw=2)
    plt.axvline(x_true,color='r',lw=2)
    plt.xlabel(r'$\log_{10}{\,\sigma}$',fontsize=16)
    plt.ylabel(r'$\log_{10}{\,\tau}$',fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    fname=dir_an+length+'_stacked_chains'+post+'_err_'+err+'.png'
    plt.savefig(fname)
    print 'Saving ', fname

# CALL THE PRINTING FUNCTION, AND LOOP 
# OVER PRIOR SETTINGS AND THREE VALUES OF ERROR 
for i in with_prior:
    for j in error:
        chain_stack = 'log10_sigma_log_10_tau_err' + str(j) + '_stacked_chains'
        print 'with prior =', i, ', err=', j, chain_stack       
        printing(i, chain_stack)
