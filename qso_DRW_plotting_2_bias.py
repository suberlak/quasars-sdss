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

prior = 'else'  # or 'no' 
length = 'long'  # or 'medium' 'long'
n_errors= 2  # also need to change number of error when calling function  err_rows_extract  below 

# Define plot size
x_min = -0.15
x_max = 2
y_min = -0.3
y_max = 2
    
if prior== 'yes':
    if length == 'short' :   # results with Javelin Prior , short LC length
        results_jav = 'qso_drw_analysis/javelin_drw_test_chain_results_with_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsshort.dat'
        n_errors = 3
    if length == 'medium' :                    # results with Javelin Prior , medium  LC length
        results_jav = 'qso_drw_medium_analysis/javelin_drw_test_chain_results_with_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsmed.dat'
        n_errors = 2
if prior =='no':
    if length == 'short' :   # results without Javelin Prior , short LC length
        results_jav = 'qso_drw_analysis/javelin_drw_test_chain_results_no_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsshort.dat'
        n_errors = 3
    if length == 'medium' :                   # results without Javelin Prior , medium  LC length  
        results_jav = 'qso_drw_medium_analysis/javelin_drw_test_chain_results_no_prior_all.txt'
        results_ch='qso_drw_S_M_Chelsea_results/fitsmed.dat'
        n_errors = 2
     
if prior == 'else' : 
    results = 'qso_drw_S_M_Chelsea_results/fitslong.dat'
    src = 'Chelsea'
    n_errors = 3
      


fig_name=['qso_drw_S_M_Chelsea_results/drw_'+length+'_prior-'+prior+'_err','_log_sigma_vs_log_tau.png']

print 'Using results from', results, 'n_errors=,', n_errors

output_ch =  np.loadtxt(results_ch, dtype='str')
output_jav = np.loadtxt(results_jav , dtype='str')


# CHELSEA RESULTS 
name_ch_raw =output_ch[:,0].astype(str)
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

# transform chelsea naming into javelin names
name_ch = name_ch_raw
for i in range(len(name_ch)):
    name_ch[i] = name_ch_raw[i][:-6]


def name_matching():
    matched_sigma_ch = np.empty(0, dtype=float)
    matched_tau_ch = np.empty(0, dtype=float)
    matched_name_ch = np.empty(0, dtype=str)
    a = min(len(name_ch),len(name_jav)) < len(name_ch)
    b = len(name_ch) == len(name_jav)
    
    if a or b :
        # case where either javelin data shorter than Chelsea, or they have same 
        # length 
        print '\n Matching Chelsea to Javelin...'
    
        for i in range(len(name_jav)):
            name = name_jav[i]
            index = np.where(name_ch == name)[0][0]
            matched_sigma_ch = np.append(matched_sigma_ch,log_10_sigma[index] )
            matched_tau_ch = np.append(matched_tau_ch,log_10_tau[index] )
            matched_name_ch = np.append(matched_name_ch, name_ch[index])
    if not a :  # case where chelsea data shorter than Javelin
        print 'Something is wrong : Chelsea shorter than Javelin...'
    return matched_sigma_ch, matched_tau_ch, matched_name_ch

 


def err_rows_extract(name_list, err_pos, n_errors): 
    '''
    Since Javelin and Chelsea data is not ordered , name lists are 
    different for each set 
    
    Form a three-column list of indices for each error: 
    col0 are indices of rows with err1
    col1 are indices of rows with err2
    col2 are indices of rows with err3
    
    Handle separately two cases, whether we have data with 
    err1,2,3 or only err1,2
    '''
    if n_errors == 3 :  
        ind=[0,0,0]
        name = name_list
        #upprint len(name)
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
            err = name[i][err_pos]  # take err id value (1 or 2) from filename      
            for j in range(1,3):
                if err == str(j) :
                    ind[j-1] = np.append(ind[j-1],i)  
        err1 = ind[0][1:]
        err2 = ind[1][1:]
        upind = [err1,err2]
    
    return upind


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

log_sigma_ch , log_tau_ch , matched_name_ch = name_matching()


for k in range(1,n_errors+1):   #n_errors+1
    # looping over err1, err2 , err3, selecting appropriate rows 

    print '\nFor err', k
     
    x_lim = [x_min, x_max]
    y_lim = [y_min,y_max]
    
    # We only load x and y that are within the limits of my histogram, to have the same 
    # pixel size for both distributions... (regardless of how much I'm removing by my boundaries)
    
    jav_indices = err_rows_extract(name_jav,-1,n_errors)
    assert len(jav_indices) == n_errors
    
    # Check that Chelsea matched rows are exactly the ones that correspond to javelin
    truth = matched_name_ch[jav_indices[0]] == name_jav[jav_indices[0]]
    assert len(np.where(truth == False)[0]) == 0 
    
    # And since I have matched Chelsea to Javelin results,  now the rows with corresponding 
    # err1,2,3, are the same 
    
    sigma =  np.power(10,log_sigma_ch)
    tau   =  np.power(10,log_tau_ch) 
    
    
    # Remove all values that are not within bounds specified 
    # since values plotted are  logs of ratios,  I'm feeding here
    # already logs of ratios 
    
    x , y, num, percent = load_x_y(sigma ,tau,jav_indices[k-1], k, x_lim, y_lim)
    
    print '\n Plotting coloured hist for log_tau  vs log_sigma  for Chelsea fitting... '
    assert     len(log_sigma_ch) == len(log_tau_ch) ==len(sigma_m_jav) ==len(tau_m_jav)
        
    # Define number of bins at the beginning, especially  if it is shared between the histograms... 
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
    
    # First histogram : make axis, and plot all that is needed  
    ax1 = fig1.add_subplot(gs[:,:90])  #     
    
    pcObject1 = ax1.pcolormesh(xedges, yedges, Hmasked)


    #plt.xlim((x_min,x_max))
    #plt.ylim((y_min,y_max))
    title = 'DRW '+ length +', err'+str(k)+', prior='+prior+', cross-matched, '+str(percent)[:5]+'% points'  
    plt.title(title)
    plt.axhline(0, color = 'r', lw=2)
    plt.axvline(0,color='r',lw=2)
    plt.ylabel(r'$\log_{10}{ \, \left(  \tau_{jav} / \tau_{ch} \right)}$',fontsize=15)
    plt.xlabel(r'$\log_{10}{ \, \left(  \sigma_{jav} / \sigma_{ch} \right)}$',fontsize=15)
    
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:,95:])
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
    # cbar.ax.set_ylabel('Counts')
    fname2 = fig_name[0]+str(k)+fig_name[1]
    plt.savefig(fname2)
    print 'File saved is ', fname2 
     
 