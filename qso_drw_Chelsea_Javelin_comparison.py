# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 17:17:48 2014

@author: astronomy

A short program to read in the results of Chelsea's and javelin fittings 

For Chelsea there are three (Short, Medium, Long), and for Javelin there are four 
(Short W/prior ,  Short No Prior , Medium W/Prior,  Medium No Prior)

It then reads the lines corresponding to appropriate errors, and calculates 
the stats (median and rms of the log (tau / tau_true).

It also plots 2D histograms of the calculated quantities, finds the maximum of the histogram,
and reports the bias of the maximum of the histogram wrt real value.

It makes tests for other possible ways of interpretation for sigma or tau from
javelin (at the beginning the definition can be changed).

"""

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf

definition = 1.0

results_ch = ['drw_S_M_L_Chelsea_results/fitsshort.dat',\
'drw_S_M_L_Chelsea_results/fitsmed.dat',\
'drw_S_M_L_Chelsea_results/fitslong.dat'] 

length_ch = ['S', 'M' , 'L ']


results_jav = ['qso_drw_analysis/javelin_drw_test_chain_results_no_prior_all.txt',\
'qso_drw_analysis/javelin_drw_test_chain_results_with_prior_all.txt',\
'qso_drw_medium_analysis/javelin_drw_test_chain_results_no_prior_all.txt',\
'qso_drw_medium_analysis/javelin_drw_test_chain_results_with_prior_all.txt']
length_jav = ['S','S','M','M']
ke = [3,3,2,2]  # how many error types there is per length type... 

dir_out = 'drw_S_M_L_Chelsea_results/'


# CHELSEA RESULTS LOAD 
def load_chelsea(results_file):
    output_ch =  np.loadtxt(results_file, dtype='str')
    name_ch =output_ch[:,0].astype(str)
    log_10_tau = output_ch[:,2].astype(np.float)       #  NOTE : already log_10
    log_10_sigma_hat = output_ch[:,3].astype(np.float)
    tau_ch = np.power(10,log_10_tau)
    sigma_hat = np.power(10,log_10_sigma_hat)
    sigma_ch = sigma_hat * np.sqrt(tau_ch / (2.0*365.0))
    K = tau_ch * np.sqrt(sigma_ch * np.sqrt(2.0*365.0))
    return name_ch, sigma_ch, tau_ch, sigma_hat, K

# JAVELIN RESULTS LOAD 
def load_javelin(results_file):
    output_jav = np.loadtxt(results_file , dtype='str')
    name_jav = output_jav[:,0].astype(str)
    #sigma_max_jav = output_jav[:,1].astype(np.float)  # extracted as max value of the chain 2D histogram 
    #tau_max_jav = output_jav[:,2].astype(np.float)
    #sigma_l_jav= output_jav[:,3].astype(np.float)
    tau_m_jav = output_jav[:,7].astype(np.float) # [days] extracted as the median value of the chain 
    sigma_m_jav = output_jav[:,4].astype(np.float)  
    if definition == 1.0 :     
        sigma_hat_jav = sigma_m_jav * np.sqrt(2.0*365.0 / tau_m_jav ) # using Chelsea's definition
        K = tau_m_jav * np.sqrt(sigma_m_jav * np.sqrt(2.0*365.0))
    if definition == 2.0 : 
        sigma_hat_jav = sigma_m_jav * np.sqrt(2.0 / tau_m_jav ) # test  definition
        K = tau_m_jav * np.sqrt(sigma_m_jav * np.sqrt(2.0))
    
    #sigma_h_jav= output_jav[:,5].astype(np.float)
    #tau_l_jav = output_jav[:,6].astype(np.float)
    #tau_h_jav = output_jav[:,8].astype(np.float)
    return name_jav, sigma_m_jav,  tau_m_jav, sigma_hat_jav, K
    
# EXTRACT ROWS WITH GIVEN ERROR VALUES
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
            #print err
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

# LOAD ONLY THOSE  X AND Y VALUES THAT MAKE SENSE ... 
def load_xyz(x_arr, y_arr, z_arr, err_indices, x_limits, y_limits,z_limits):
    indices = err_indices 
    print '\n Loading x and y and z ... ' 
   
    x = x_arr[indices]  
    y = y_arr[indices]
    z = z_arr[indices]
    # sieve out suspiciously bad values , based only on x and y 
    xinf = np.asarray(map(isinf,x),dtype=bool)
    yinf = np.asarray(map(isinf,y),dtype=bool)
    zinf = np.asarray(map(isinf,z),dtype=bool)
    ttlinf = xinf + yinf + zinf
    # ttlwh = np.where(ttlinf == True)  list of good indices
    gi = -ttlinf  # good_indices 
    zsmall = np.where(z < z_limits[0])        
    zlarge = np.where(z > z_limits[1])    
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
    
    return x[gi], y[gi],z[gi],  non_inf, percent

# STATS FUNCTIONS 
#def rmse(predictions, targets):
#    return np.sqrt(((predictions - targets) ** 2).mean())
    
def give_stats(fit_values,true_value):
    ''' Takes an array of fit values  , 
        and a scalar true value for fitted 
        parameter 
    '''   
    median = np.median(np.log10(fit_values)) - np.log10(true_value)
    rms = (np.percentile(np.log10(fit_values), 75) - np.percentile(np.log10(fit_values), 25)) * 0.7413
    return  median, rms

# PLOTTING AND CALCULATING THE BIAS 
def histogram(x_arr, y_arr, number, percent, xlim, ylim, title, k, size, *args):
    # args could include javelin results_file , from which you can 
    # take the info about the prior  
    
    length_label = size 
        
    
    if title == 'jav':
        prior = args[0]
        print 'An argument for prior taken is ', args[0]
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
    a,b = np.where(H == H.max())    
    x_max = xedges[a[0]]
    y_max = yedges[b[0]]
    
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    
    x_true = np.log10(0.2)
    y_true = 2.0
    print '\n The true values are (',  x_true, ',', y_true, ')'
    print 'The maximum of the log(sigma_ch) vs log(tau_ch) is (', x_max,', ', y_max, ')'
    del_x =  x_max - x_true
    del_y =  y_max - y_true 
    mse = np.sqrt((del_x**2.0) + (del_y **2.0))
    print 'The bias is thus  (', del_x ,',', del_y, ')'
    print 'The MSE bias is ' , mse
    

    ax1 = fig1.add_subplot(gs[:,:90])   
    pcObject1 = ax1.pcolormesh(xedges, yedges, Hmasked)
    
    xmin = np.log10(x_lim[0])
    xmax = np.log10(x_lim[1])
    ymin =  np.log10(y_lim[0])
    ymax =  np.log10(y_lim[1])
    
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    
    plt.axhline(2, color = 'r', lw=2)
    plt.axvline(np.log10(0.2),color='r',lw=2)
    plt.axhline(y_max, color='b' , lw=2 )
    plt.axvline(x_max,color='b' , lw=2  )
    x_label_ch = r'$\log_{10}{ \, \left(  \sigma_{ch} \right)}$'
    y_label_ch = r'$\log_{10}{ \, \left( \tau_{ch} \right)}$'
    x_label_jav = r'$\log_{10}{ \, \left(  \sigma_{jav} \right)}$'
    y_label_jav = r'$\log_{10}{ \, \left(  \tau_{jav} \right)}$'
     
    
    if title == 'ch' : 
        plt.ylabel(y_label_ch,fontsize=15)
        plt.xlabel(x_label_ch,fontsize=15)
        title_hist = 'DRW Chelsea , err '+str(k+1) +', '+ str(number) + ', i.e.  ' + str(percent)+ '% points'
        fname = dir_out + 'drw_log_sigma_log_tau_'+title+'_'+length_label+'_err_'+str(k+1)+'.png'
    else:
        plt.ylabel(y_label_jav,fontsize=15)
        plt.xlabel(x_label_jav,fontsize=15)
        title_hist = 'DRW Javelin , err '+str(k+1) +', '+ str(number) + ', i.e.  ' + str(percent)+ '% points'
        fname = dir_out + 'drw_log_sigma_log_tau_'+title+'_'+length_label+'_err_'+str(k+1)+'_'+prior+'_prior.png' 
        
    
    plt.title(title_hist)
    # Add the colorbar  
    axC = fig1.add_subplot(gs[:,95:])
    cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
    cbar.ax.set_ylabel('Counts')
    
    plt.savefig(fname)
    print 'File saved is ', fname

# CHECK IF THERE IS PRIOR OR NOT.. 
def check_prior(results_file):
    if results_file[-16:-14] == 'no' :
        print 'We are using javelin results with no prior...'
        prior = 'no'
    else:
        print 'We are using javelin results with prior ...' 
        prior = 'yes'
    return prior

# SET CONSTRAINTS ON THE DATA... 
def set_limits(x_values,y_values,z_values, index_list):
    # x_values and y_values are sigma and tau 
    # that are not constrained by anything
    # this function checks whether the constraint 
    # below is taking away more than 95% of points...
    x_min = 0.01
    x_max = 10. 
    y_min = 1. 
    y_max = 1000.
    z_min = 0.01
    z_max = 100
    
#    x = x_values[index_list]
#    y = y_values[index_list]
#    
#    if x_max < np.percentile(x,95):
#        x_max = np.percentile(x,97)
#        print 'We had to change the x_max from', 10.0, ' to ', x_max
#    if y_max < np.percentile(y,95):
#        y_max = np.percentile(y, 97)
#        print 'We had to change the y_max from', 1000.0, ' to ', y_max
#    if x_min > np.percentile(x, 5):
#        x_min = np.percentile(x,3)
#        print 'We had to change the x_min from', 0.01, ' to ', x_min
#    if y_min > np.percentile(y, 5):
#        y_min = np.percentile(y,3)
#        print 'We had to change the x_min from', 1.0, ' to ', y_min
    
    x_lim = [x_min, x_max]
    y_lim = [y_min,y_max]
    z_lim = [z_min,z_max]
    return  x_lim, y_lim , z_lim



# MAIN LOOP - CALCULATING THE STATS, AND MAKING A HISTOGRAM, AS WELL AS 
# CALCULATING THE BIAS. 
#
#tau_median      = np.empty(0,dtype=float)
#sigma_median    = np.empty(0,dtype=float)
#tau_rms         = np.empty(0,dtype=float)
#sigma_rms       = np.empty(0,dtype=float)
#error           = np.empty(0,dtype=float)
#prior_flag      = np.empty(0,dtype=str)
#distr_log_sigma_max = np.empty(0,dtype=float)
#distr_log_tau_max   = np.empty(0,dtype=float)
#bias_sigma      = np.empty(0,dtype=float)
#bias_tau        = np.empty(0,dtype=float)
#MSE_bias        = np.empty(0,dtype=float)
#outfilename     = np.empty(0,dtype=float)

res_ch = []

err_pos = [-7,-7,-5]  # because Chelsea changes her naming schemes..
tau_true = 100.0
sigma_true = 0.2 


sigma_hat_true = sigma_true * np.sqrt(2.0 * 365.0 / tau_true) # Chelsea's definition
K_true = tau_true * np.sqrt(sigma_true*np.sqrt(365.0 * 2.0))


#  CHELSEA RESULTS  
for j in range(0,len(results_ch)): # 
    title = 'ch'
    print 'CHELSEA RESULTS'
    print '\nUsing ', results_ch[j]
    print 'Its length is ', length_ch[j]
    
    name_ch,  sigma_ch, tau_ch, sigma_hat, K_ch = load_chelsea(results_ch[j])
    ch_ind_list = err_rows_extract(name_ch, err_pos[j], 3)    
    
    for k in range(3):
        print '\n Error : k=', k
        # x_lim, y_lim , z_lim = set_limits(sigma_ch, tau_ch, sigma_hat, ch_ind_list[k]) # now doesn't do anything
        # s_ch, t_ch, s_ch_hat, num_ch, percent_ch = load_xyz(sigma_ch, tau_ch, sigma_hat, ch_ind_list[k], x_lim, y_lim, z_lim)
        
        # for now skip setting limits, and loading xyz according to its limits 
        # just load everything, and hope that the stats will not be sensitive 
        # to outliers    (histogram would be,  but I am not plotting the 
        # histogram now at all !)          
        
        s_ch = sigma_ch[ch_ind_list[k]]
        t_ch = tau_ch[ch_ind_list[k]]
        s_ch_hat = sigma_hat[ch_ind_list[k]]
        K = K_ch[ch_ind_list[k]]
        
#        delta_b = 2.0* (np.percentile(t_ch, 75) - np.percentile(t_ch, 25)) / np.power(len(t_ch),1.0/3.0 )
#        print 'delta b based on tau' , delta_b
#        delta_b_2 = 2.0* (np.percentile(s_ch, 75) - np.percentile(s_ch, 25)) / np.power(len(s_ch),1.0/3.0 )
#        print 'delta b based on sigma', delta_b_2
#        print 'current bin width tau', (t_ch.max() - t_ch.min()) / 50
#        print 'current bin width sigma', (s_ch.max() - s_ch.min()) / 50
        
        t_med_ch, t_rms_ch = give_stats(t_ch,tau_true)    
        print '\n Chelsea: log(tau_fit) -  log(tau_true)'
        print ' median: ', t_med_ch,'rms:', t_rms_ch 

        s_med_ch, s_rms_ch = give_stats(s_ch,sigma_true)    
        print ' \n log(sigma_fit) - log(sigma_true)'
        print ' median: ', s_med_ch,'rms:', s_rms_ch 
        
        s_hat_med, s_hat_rms = give_stats(s_ch_hat, sigma_hat_true)
        print ' \n log(sigma_hat_fit) - log(sigma_hat_true)'
        print  'median: ' , s_hat_med, 'rms: ', s_hat_rms
        
        K_med, K_rms = give_stats(K,K_true)        
        print ' \n log(K_fit) - log(K_true)'
        print  'median: ' , K_med, 'rms: ', K_rms
        
        #tau_median = np.append(tau_median, t_med_ch)
        #tau_rms    = np.append(tau_rms , t_rms_ch)
        #sigma_median = np.append(sigma_median, s_med_ch)
        #sigma_rms    = np.append(sigma_rms , s_rms_ch)
        # histogram(s_ch, t_ch, num_ch, percent_ch, x_lim, y_lim, 'ch', k,length_ch[j])    
        res_ch.append([t_med_ch, t_rms_ch, s_med_ch, s_rms_ch, s_hat_med, s_hat_rms,K_med, K_rms])

res_jav = []

# JAVELIN  RESULTS  
for j in range(0,len(results_jav)):
    title = 'jav' 
    print 'JAVELIN RESULTS'
    print '\nUsing ', results_jav[j]
    print 'Its length is ', length_jav[j]
    
    name_jav,  sigma_jav, tau_jav , sigma_hat_jav, K_jav = load_javelin(results_jav[j])
    jav_ind_list = err_rows_extract(name_jav, -1, ke[j])    
    prior = check_prior(results_jav[j])

    for k in range(0,ke[j]):
        print '\n Error : k=', k
        prior_err_flag = prior +', err' +  str(k+1)
        #x_lim, y_lim, z_lim = set_limits(sigma_jav, tau_jav, sigma_hat_jav, jav_ind_list[k])
        #s_jav, t_jav, s_hat_jav, num_jav, percent_jav = load_xyz(sigma_jav, tau_jav, sigma_hat_jav, jav_ind_list[k], x_lim, y_lim, z_lim)
     
        s_jav = sigma_jav[jav_ind_list[k]]
        t_jav = tau_jav[jav_ind_list[k]]
        s_hat_jav = sigma_hat_jav[jav_ind_list[k]]
        K = K_jav[jav_ind_list[k]]
      
        t_med_jav , t_rms_jav = give_stats(t_jav,tau_true)
        print '\n Javelin: log(tau_fit / tau_true)'
        print ' median: ', t_med_jav,'rms:', t_rms_jav
        
        s_med_jav , s_rms_jav = give_stats(s_jav,sigma_true)
        print '\n log(sigma_fit / sigma_true)'
        print ' median: ', s_med_jav,'rms:', s_rms_jav
        
        s_hat_med, s_hat_rms = give_stats(s_hat_jav, sigma_hat_true)
        print ' \n log(sigma_hat_fit) - log(sigma_hat_true)'
        print  'median: ' , s_hat_med, 'rms: ', s_hat_rms
        
        K_med, K_rms = give_stats(K,K_true)        
        print ' \n log(K_fit) - log(K_true)'
        print  'median: ' , K_med, 'rms: ', K_rms
        
        # histogram(s_jav, t_jav, num_jav, percent_jav, x_lim, y_lim, 'jav', k, length_jav[j], prior)
        res_jav.append([prior_err_flag, t_med_jav, t_rms_jav, s_med_jav,s_rms_jav, s_hat_med, s_hat_rms, K_med, K_rms])
        
        
# SAVE RESULTS TO TXT 
results_ch_arr = np.asarray(res_ch)        
np.savetxt('drw_chelsea_results_stats.txt',results_ch_arr,delimiter=" ", fmt="%s")        

results_jav_arr = np.asarray(res_jav)
np.savetxt('drw_javelin_results_stats.txt',results_jav_arr,delimiter=" ", fmt="%s")

print 'I have saved the results to two files  : now use  ',\
'qso_drw_Chelsea_comparison_latex.py to make the latex output'