# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 14:48:12 2014

@author: suberlak

Plotting stacked chain posterior distributions to investigate their gaussianity


Chains are named  :
 ch_DRWtest_LC547_err2.dat_chain.dat 

"""

import numpy as np 
from math import  isinf

dir_in = ['qso_drw_chains/','qso_drw_chains/no_prior/','qso_drw_medium_LC_chains/with_prior/','qso_drw_medium_LC_chains/no_prior/']

dir_out= ['qso_drw_analysis/','qso_drw_medium_analysis/']

out = dir_out[1]
dir_with_prior = dir_in[2]
dir_no_prior   = dir_in[3]

with_prior = 'no'
chain_list = 'chains_err_all.ls'

if with_prior == 'yes' : 
    filename = dir_with_prior +  chain_list
    files=np.loadtxt(filename,dtype=str)
    post = '_with_prior'
else : 
    dir_chains = dir_in[3]
    filename = dir_no_prior +  chain_list
    files=np.loadtxt(filename,dtype=str)
    post = '_no_prior'
    


#########################
#     LOADING CHAINS    #
#########################

log10_sigma = np.zeros((1,3),dtype=float)
log10_tau = np.zeros((1,3),dtype=float)


log10_sigma_err1 = np.empty(0,dtype=float)
log10_tau_err1   = np.empty(0,dtype=float)

log10_sigma_err2 = np.empty(0,dtype=float)
log10_tau_err2   = np.empty(0,dtype=float)

log10_sigma_err3 = np.empty(0,dtype=float)
log10_tau_err3   = np.empty(0,dtype=float)

#np.empty(0,dtype=float)
for j in range(len(files)):   # len(files)

    # LOAD THE CHAIN  
    fchain = dir_chains+files[j]
    flatchain= np.genfromtxt(fchain)
    drw_name = files[j][3:-14]
    error =  files[j][-15]          
    print  '\n',drw_name, 'error', error
    
    # RETRIEVE SIGMA AND TAU COLS
    sigma = np.exp(flatchain[:,0]) 
    tau = np.exp(flatchain[:,1])
    x=np.log10(sigma)
    y=np.log10(tau)
    
    # REMOVE THE INFINITE VALUES
    xinf = np.asarray(map(isinf,x),dtype=bool)
    yinf = np.asarray(map(isinf,y),dtype=bool)
    ttlinf = xinf + yinf
    # ttlwh = np.where(ttlinf == True)  list of good indices
    gi = -ttlinf  # good_indices 
    non_inf = len(np.where(gi == True)[0])
    
    print 'Out of ', len(x),' rows, we have ', non_inf, ' of those that do not',\
    ' have any infinities, and only those are used '
    
    # STACK INTO THREE SEPARATE FILES DEPENDING ON ERROR 
    if error == '1' : 
        log10_sigma_err1 = np.append(log10_sigma_err1,x[gi])
        log10_tau_err1 = np.append(log10_tau_err1,y[gi]) 
    if error == '2' :
        log10_sigma_err2= np.append(log10_sigma_err2,x[gi])
        log10_tau_err2= np.append(log10_tau_err2,y[gi]) 
    if error == '3' : 
        log10_sigma_err3 =np.append(log10_sigma_err3,x[gi])
        log10_tau_err3 = np.append(log10_tau_err3,y[gi]) 
    

print 'We have for err1 ', len(log10_sigma_err1), 'values for sigma results, and ',\
      len(log10_tau_err1), 'values for tau results '
      

print 'We have for err2 ', len(log10_sigma_err2), 'values for sigma results, and ',\
      len(log10_tau_err2), 'values for tau results '

print 'We have for err3 ', len(log10_sigma_err3), 'values for sigma results, and ',\
      len(log10_tau_err3), 'values for tau results '

###########################
# SAVING STACKED CHAINS #
###########################
if len(log10_sigma_err1) > 0: 
    arr1 = [log10_sigma_err1, log10_tau_err1]
    fname = out + 'log10_sigma_log_10_tau_err1_stacked_chains'+post+'.npy'
    np.save(fname,arr1)
    print 'Saving ', fname 

if len(log10_sigma_err2) > 0: 
    arr2 = [log10_sigma_err2, log10_tau_err2]
    fname = out + 'log10_sigma_log_10_tau_err2_stacked_chains'+post+'.npy'
    np.save(fname,arr2)
    print 'Saving ', fname 
    
if len(log10_sigma_err3) > 0: 
    arr3 = [log10_sigma_err3, log10_tau_err3]
    fname = out + 'log10_sigma_log_10_tau_err3_stacked_chains'+post+'.npy'
    np.save(fname,arr3)
    print 'Saving ', fname 

