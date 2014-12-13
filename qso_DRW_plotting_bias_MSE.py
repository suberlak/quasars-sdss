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
#results='qso_drw_analysis/javelin_drw_test_chain_results_no_prior_all.txt'
results = 'qso_drw_S_M_Chelsea_results/fitslong.dat'
fname='qso_drw_S_M_Chelsea_results/fitslong_plot_err_'
end = '.png'


output_ch =  np.loadtxt(results, dtype='str')

# CHELSEA RESULTS 
name_ch =output_ch[:,0].astype(str)
log_10_tau = output_ch[:,2].astype(np.float)       #  NOTE : already log_10
log_10_sigma_hat = output_ch[:,3].astype(np.float)
tau = np.power(10,log_10_tau)
sigma_hat = np.power(10,log_10_sigma_hat)
sigma = sigma_hat * np.sqrt(tau / (2.0*365.0))
log_10_sigma  = np.log10(sigma)




##  Chelsea : err_pos = -7
##  Javelin : err_pos = -1 
def err_rows_extract(name_list, err_pos, n_errors):
    if n_errors == 3 : 
        
        ind=[0,0,0]
        name = name_list
        #print len(name)
        for i in range(0,len(name)):
            err = name[i][err_pos]  # take err id value (1,2, or 3) from filename 
            
            for j in range(1,4):
                if err == str(j) :
                    ind[j-1] = np.append(ind[j-1],i)
        
        err1 = ind[0][1:]
        err2 = ind[1][1:]
        err3 = ind[2][1:]
        
        upind = [err1,err2,err3]
    if n_errors == 2:  
        ind=[0,0]
        name = name_list
        # print len(name)
        for i in range(0,len(name)):
            err = name[i][err_pos]  # take err id value (1,2, or 3) from filename 
            
            for j in range(1,3):
                if err == str(j) :
                    ind[j-1] = np.append(ind[j-1],i)
        
        err1 = ind[0][1:]
        err2 = ind[1][1:]
        
        
        upind = [err1,err2]
    
    return upind

upind = err_rows_extract(name_ch,-7,3)


def load_x_y(x_arr, y_arr, err_indices, ka , x_limits, y_limits):
    indices = err_indices 
    print '\n Loading x and y ... ' 
   
    x = x_arr[indices]  
    y = y_arr[indices]
    
    # sieve out suspiciously bad values , based only on x and y 
    
    if ka < 0 :  
        xinf = np.asarray(map(isinf,x),dtype=bool)
        yinf = np.asarray(map(isinf,y),dtype=bool)
        ttlinf = xinf + yinf
        # ttlwh = np.where(ttlinf == True)  list of good indices
        gi = -ttlinf  # good_indices 
       
        non_inf = len(np.where(gi == True)[0])
    else :   #(ALWAYS)
        # separate treatment of the high error test  
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

for k in range(1,4):
    print '\n Plotting coloured hist for log_tau  vs log_sigma  for javelin fitting ' 
    plt.clf()
    fig1 = plt.figure()
    
     # Define plot size
    x_min =-1.1
    x_max = 0.4
    y_min = 1.2
    y_max=4
    
    x_lim = [x_min, x_max]
    y_lim = [y_min,y_max]
    
    # We only load x and y that are within the limits of my histogram, to have the same 
    # pixel size for both distributions... (regardless of how much I'm removing by my boundaries)
    
    x_ch, y_ch,  num_jav, percent = load_x_y(log_10_sigma,log_10_tau,upind[k-1], k, x_lim, y_lim)
    
    
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
    H, xedges,yedges = np.histogram2d(x_ch,y_ch,bins=nbins)
    a,b = np.where(H == H.max())    
    x_max = xedges[a[0]]
    y_max = yedges[b[0]]
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    # First histogram : make axis, and plot all that is needed  
    ax1 = fig1.add_subplot(gs[15:,:])  #     
    
    pcObject = ax1.pcolormesh(xedges, yedges, Hmasked)
    

    plt.xlim((x_min,x_max))
    plt.ylim((y_min,y_max))
    title = 'DRW err'+str(k)+', '+str(percent)+ ' obj, max values '
    plt.title(title)
    plt.axhline(np.log10(100))
    plt.axvline(np.log10(0.2))
    plt.ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15)
    plt.xlabel(r'$\log_{10}{ \, \sigma_{ch}}$',fontsize=15)
    
  
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:5,:])
    cbar = fig1.colorbar(pcObject,ax=ax1, cax=axC, orientation='horizontal')
    # cbar.ax.set_ylabel('Counts')
    fname2 = fname+str(k)+end
    plt.savefig(fname2)
    print 'File saved is ', fname2 
     
 