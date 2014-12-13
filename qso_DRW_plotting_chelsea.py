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

prior = 'no'  # or 'no' 
length = 'medium'  # or 'medium'
n_errors= 2  # also need to change number of error when calling function  err_rows_extract  below 
if prior== 'yes':
    if length == 'short' :   # results with Javelin Prior , short LC length
        results_jav = 'qso_drw_analysis/javelin_drw_test_chain_results_with_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsshort.dat'
    else:                    # results with Javelin Prior , medium  LC length
        results_jav = 'qso_drw_medium_analysis/javelin_drw_test_chain_results_with_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsmed.dat'
else : 
    if length == 'short' :   # results without Javelin Prior , short LC length
        results_jav = 'qso_drw_analysis/javelin_drw_test_chain_results_no_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsshort.dat'
    else:                    # results without Javelin Prior , medium  LC length  
        results_jav = 'qso_drw_analysis/javelin_drw_test_chain_results_no_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsmed.dat'
        


fig_name=['qso_drw_S_M_Chelsea_results/drw_'+length+'_prior-'+prior+'_err','_log_sigma_vs_log_tau.png']

print 'Using Chelsea results from', results_ch,' and Javelin results from ', results_jav

output_ch =  np.loadtxt(results_ch, dtype='str')
output_jav = np.loadtxt(results_jav , dtype='str')


# CHELSEA RESULTS 
name_ch =output_ch[:,0].astype(str)
log_10_tau = output_ch[:,2].astype(np.float)       #  NOTE : already log_10
log_10_sigma_hat = output_ch[:,3].astype(np.float)
tau = np.power(10,log_10_tau)
sigma_hat = np.power(10,log_10_sigma_hat)
sigma = sigma_hat * np.sqrt(tau / (2.0*365.0))
log_10_sigma  = np.log10(sigma)

# JAVELIN RESULTS 
name_jav =output_jav[:,0].astype(str)
sigma_max_jav = output_jav[:,1].astype(np.float)  # extracted as max value of the chain 2D histogram 
tau_max_jav = output_jav[:,2].astype(np.float)
sigma_l_jav= output_jav[:,3].astype(np.float)
sigma_m_jav_raw = output_jav[:,4].astype(np.float)    # extracted as the median value of the chain 
sigma_m_jav = sigma_m_jav_raw 

sigma_h_jav= output_jav[:,5].astype(np.float)
tau_l_jav = output_jav[:,6].astype(np.float)
tau_m_jav = output_jav[:,7].astype(np.float)
tau_h_jav = output_jav[:,8].astype(np.float)



# HISTOGRAM  OF   LOG SIGMA  VS LOG TAU,  LIKE FIG 13 FROM MCLEOD+2011 
# SINCE ALL THE  ROWS INCLUDE ALL ERR VALUES, NEED TO FIRST SELECT THOSE 
# ROWS THAT HAVE ERR1 , ERR2, ERR3 SEPARATELY 



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

ch_indices = err_rows_extract(name_ch,-7,3)
jav_indices = err_rows_extract(name_jav,-1,2)


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
        
    
    print 'Out of ', len(x),' rows, we have ', non_inf, ' of those that do not',\
    ' have any infinities, and only those are used for plotting '
    
    return x[gi], y[gi], non_inf


for k in range(1,n_errors+1):  # looping over err1, err2 , err3, selecting appropriate rows 
    #global prior 

    print '\nFor err', k
    
     # Define plot size
    x_min =-1.1
    x_max = 0.4
    y_min = 1.2
    y_max=4
    
    x_lim = [x_min, x_max]
    y_lim = [y_min,y_max]
    
    # We only load x and y that are within the limits of my histogram, to have the same 
    # pixel size for both distributions... (regardless of how much I'm removing by my boundaries)
    
    x_ch, y_ch,  num_jav = load_x_y(log_10_sigma,log_10_tau,ch_indices[k-1], k, x_lim, y_lim)
    x_jav, y_jav, num_ch = load_x_y(np.log10(sigma_m_jav),np.log10(tau_m_jav), jav_indices[k-1], k, x_limits=x_lim, y_limits=y_lim )
    print '\n Plotting coloured hist for log_tau  vs log_sigma  for Chelsea fitting '
    
    print ' For Chelsea plot we use ', len(x_ch), 'sigma values and ',  len(y_ch), ' , tau values'
    print ' We Javelin plot we use ', len(x_jav), 'sigma values and ',  len(y_jav), ' , tau values'
    
    # Define number of bins at the beginning, especially  if it is shared between the histograms... 
    nbins =40
    plt.clf()
    fig1 = plt.figure()
    
    # Define the canvas to work on   and the  grid  
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)    
    
   
    
    # First histogram  : Chelsea results 
    H, xedges,yedges = np.histogram2d(x_ch,y_ch,bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    # First histogram : make axis, and plot all that is needed  
    ax1 = fig1.add_subplot(gs[15:,:48])  #     
    
    pcObject1 = ax1.pcolormesh(xedges, yedges, Hmasked)


    plt.xlim((x_min,x_max))
    plt.ylim((y_min,y_max))
    title = 'DRW err'+str(k)+', '+str(num_ch)+ ' obj, Chelsea '
    plt.title(title)
    plt.axhline(np.log10(100), lw=2)
    plt.axvline(np.log10(0.2),lw=2)
    plt.ylabel(r'$\log_{10}{ \, \tau_{ch}}$',fontsize=15)
    plt.xlabel(r'$\log_{10}{ \, \sigma_{ch}}$',fontsize=15)
    
    
    # Second histogram : Javelin results 
    H, xedges,yedges = np.histogram2d(x_jav,y_jav,bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    # Second histogram : make axis, and plot all that is needed  
    ax2 = fig1.add_subplot(gs[15:,53:]) #     
    pcObject2 = ax2.pcolormesh(xedges, yedges, Hmasked)
    
    plt.xlim((x_min,x_max))
    plt.ylim((y_min,y_max))
    title = 'DRW err'+str(k)+', '+str(num_jav)+ ' obj Javelin, prior='+ prior
    plt.title(title)
    ax2.set_ylabel(r'$\log_{10}{ \, \tau_{JAV, median}}$',fontsize=15 )
    ax2.set_xlabel(r'$\log_{10}{ \, \sigma_{JAV, median}}$',fontsize=15)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_coords(1.15,0.5)
    plt.axhline(np.log10(100), lw=2)
    plt.axvline(np.log10(0.2),lw=2)
    
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:5,:])
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='horizontal')
    # cbar.ax.set_ylabel('Counts')
    fname2 = fig_name[0]+str(k)+fig_name[1]
    plt.savefig(fname2)
    print 'File saved is ', fname2 
     
 