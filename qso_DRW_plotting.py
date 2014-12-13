# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 17:48:01 2014

@author: suberlak

plotting   log(tau)  vs  log(sigma)  for the 
simulated DRW  sample 

"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf

# set_prior = TRUE in Javelin
#results='qso_drw_analysis/javelin_drw_test_chain_results_with_prior_all.txt'

# set_prior  = FALSE  in Javelin 
#results='qso_drw_analysis/javelin_drw_test_chain_results_no_prior_603.txt'
results='qso_drw_analysis/javelin_drw_test_chain_results_no_prior_all.txt'

fname='qso_drw_analysis/drw_no_prior_err'
fname1 ='_log_sigma_vs_log_tau.png'


output =  np.loadtxt(results, dtype='str')

name =output[:,0].astype(str)
sigma_max = output[:,1].astype(np.float)
tau_max = output[:,2].astype(np.float)
sigma_l= output[:,3].astype(np.float)
sigma_m= output[:,4].astype(np.float)
sigma_h= output[:,5].astype(np.float)
tau_l = output[:,6].astype(np.float)
tau_m = output[:,7].astype(np.float)
tau_h = output[:,8].astype(np.float)


# HISTOGRAM  OF   LOG SIGMA  VS LOG TAU,  LIKE FIG 13 FROM MCLEOD+2011 
ind=[0,0,0]
for i in range(len(name)):
    a = name[i][-1]
    for j in range(1,4):
        if a == str(j) :
            ind[j-1] = np.append(ind[j-1],i)

err1 = ind[0][1:]
err2 = ind[1][1:]
err3 = ind[2][1:]

upind = [err1,err2,err3]

for k in range(1,4):
    print '\n Plotting coloured hist for log_tau  vs log_sigma  for javelin fitting ' 
    plt.clf()
    fig1 = plt.figure()
    x = np.log10(sigma_max[upind[k-1]])
    y = np.log10(tau_max[upind[k-1]])
    
    x1 = np.log10(sigma_m[upind[k-1]])
    y1 = np.log10(tau_m[upind[k-1]])
    
    
    # sieve out suspiciously bad values , based only on x and y 
    
    if k < 1 :
        xinf = np.asarray(map(isinf,x),dtype=bool)
        yinf = np.asarray(map(isinf,y),dtype=bool)
        ttlinf = xinf + yinf
        # ttlwh = np.where(ttlinf == True)  list of good indices
        gi = -ttlinf  # good_indices 
       
        non_inf = len(np.where(gi == True)[0])
    else :
        # separate treatment of the high error test  
        xinf = np.asarray(map(isinf,x),dtype=bool)
        yinf = np.asarray(map(isinf,y),dtype=bool)
        ttlinf = xinf + yinf
        # ttlwh = np.where(ttlinf == True)  list of good indices
        gi = -ttlinf  # good_indices 
        ysmall = np.where(y < 0)
        xsmall = np.where(x < -1.4)
        gi[xsmall] = False
        gi[ysmall] = False
        non_inf = len(np.where(gi == True)[0])
        
    
    print 'Out of ', len(x),' rows, we have ', non_inf, ' of those that do not',\
    ' have any infinities, and only those are used for plotting '
    
    
    # Define number of bins at the beginning, especially  if it is shared between the histograms... 
    nbins =60
    
    
    # Define the canvas to work on   and the  grid  
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)    
    
    # Define plot size
    x_min =-1.2
    x_max = 0.8
    y_min = 0
    y_max=4
    
    # First histogram :  making the histogram values 
    H, xedges,yedges = np.histogram2d(x[gi],y[gi],bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    # First histogram : make axis, and plot all that is needed  
    ax1 = fig1.add_subplot(gs[15:,:48])  #     
    
    pcObject = ax1.pcolormesh(xedges, yedges, Hmasked)
    

    plt.xlim((x_min,x_max))
    plt.ylim((y_min,y_max))
    title = 'DRW err'+str(k)+', '+str(non_inf)+ ' obj, max values '
    plt.title(title)
    plt.axhline(np.log10(100))
    plt.axvline(np.log10(0.2))
    plt.ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15)
    plt.xlabel(r'$\log_{10}{ \, \sigma_{ch}}$',fontsize=15)
    
    # Second histogram : making the histogram values 
    H, xedges,yedges = np.histogram2d(x1[gi],y1[gi],bins=nbins)
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
    plt.xlim((x_min,x_max))
    plt.ylim((y_min,y_max))
    title = 'DRW err'+str(k)+', '+str(non_inf)+ ' obj, median values '
    plt.title(title)
    plt.axhline(np.log10(100))
    plt.axvline(np.log10(0.2))
    # plt.ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15)
    
    
    
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:5,:])
    cbar = fig1.colorbar(pcObject,ax=ax1, cax=axC, orientation='horizontal')
    # cbar.ax.set_ylabel('Counts')
    fname2 = fname+str(k)+fname1
    plt.savefig(fname2)
    print 'File saved is ', fname2 
     
 