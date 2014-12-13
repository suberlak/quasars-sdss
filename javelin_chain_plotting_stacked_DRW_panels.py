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
from matplotlib.gridspec import GridSpec

dir_in = ['qso_drw_analysis/','qso_drw_medium_analysis/']

dir_out= ['qso_drw_analysis/','qso_drw_medium_analysis/']


dir_an = dir_in[0]  # CHOOSE WHICH ANALYSIS DIR WE NEED 
length= 'medium'
# ARRAYS WITH ITERATED PARAMETERS 
error = [1,2]  # need to add 3 if also considering error 3
with_prior = ['no', 'yes']  # need to add 'yes' if we have both cases 

# DEFINE A FUNCTION TO PRINT STACKED CHAINS FOR A
# GIVEN  ERROR AND PRIOR SETTING 

def load_x_y(prior, err):
    global dir_an
     
    chain_stack = 'log10_sigma_log_10_tau_err' + str(err) + '_stacked_chains'     
     
    if prior == 'yes' : 
        filename = dir_an + chain_stack+'_with_prior.npy'
        files=np.load(filename)
        post = '_with_prior'
    else : 
        filename = dir_an + chain_stack+'_no_prior.npy'
        files=np.load(filename)
        post = '_no_prior'    
    
    err  = chain_stack[26]   # take it out of name 
    
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
    return x, y, post , err 
           
           
           
           

def printing(err_in):
    
    # Define number of bins at the beginning, especially  if it is shared between the histograms... 
    nbins =100
    # the true values of sigma and tau     
    y_true = np.log10(100)
    x_true = np.log10(0.2)
    x_try  = np.log10(0.2 * np.sqrt(2.0)) # if perhaps I took Javelin's output as sigma, and in 
      # reality it is SF_inf .... 
    
    
    # Define the canvas to work on   and the  grid  
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)    

    # FIRST HISTOGRAM : NO PRIOR :   making the histogram values 
    x,y,post,err = load_x_y(with_prior[0], err_in)    ## CHANGE ERR 
    H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    # First histogram : make axis, and plot all that is needed  
    ax1 = fig1.add_subplot(gs[15:,:48])       
    pcObject = ax1.pcolormesh(xedges, yedges, Hmasked)
    

    #plt.xlim((x_min,x_max))
    #plt.ylim((y_min,y_max))
    title = 'Stacked chains '+length+', prior=' + with_prior[0] +', err '+ err
    plt.title(title)
    plt.axhline(y_true,color='r',lw=2)
    plt.axvline(x_true,color='r',lw=2)
    plt.axvline(x_try,color='b',lw=2)
    plt.ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15)
    plt.xlabel(r'$\log_{10}{ \, \sigma_{ch}}$',fontsize=15)
    
    # SECOND HISTOGRAM : WITH PRIOR : making the histogram values 
    x,y,post,err = load_x_y(with_prior[1], err_in)
    H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    # Second histogram : make axis 
    ax2 = fig1.add_subplot(gs[15:,53:])
    #ax2.set_yticklabels('',visible=False)
    ax2.set_ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15 )
    ax2.set_xlabel(r'$\log_{10}{ \, \sigma_{ch}}$',fontsize=15)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_coords(1.15,0.5)
    pcObject1 = ax2.pcolormesh(xedges, yedges, Hmasked)
    #plt.xlim((x_min,x_max))
    #plt.ylim((y_min,y_max))
    title = 'Stacked chains '+length+', prior=' + with_prior[1] +', err '+ err
    plt.title(title)
    plt.axhline(y_true,color='r',lw=2)
    plt.axvline(x_true,color='r',lw=2)
    plt.axvline(x_try,color='b',lw=2)
    # plt.ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15)
    
    
    
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:5,:])
    cbar = fig1.colorbar(pcObject,ax=ax1, cax=axC, orientation='horizontal')
    # cbar.ax.set_ylabel('Counts')
    fname=dir_an+length+'_stacked_chains_panels_err_'+err+'.png'
    plt.savefig(fname)
    print 'File saved is ', fname 
    
    
    
# CALL MY FUNCTIONS, EFFORTLESSLY LOOPING OVER ERROR VALUES     
    
for err_in in error:     
    printing(err_in)

