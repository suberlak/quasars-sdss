# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 17:01:19 2014

@author: suberlak

A new program that does not require Javelin to run 

Loading chains bit .

REQUIREMENTS:

Before running you need to update the chain list file (in terminal)


NOTE:

The chain name has to have the form  

chain='ch_xxxxxxxxx.dat_chain.dat'  
because we are taking the star name from  chain[3:-14]


OUTPUT: 

A text file like javelin_chain_results_new8_sigma_tau.txt    ,  with columns

QSO_name           sigma_l,       sigma_m,       sigma_h,       tau_l,        tau_m,       tau_h            
222857.83-010734.8 0.147325389744 0.178587132035 0.227150026128 70.2886684856 116.70354396 209.732145982

NOTE : 

THE SCRIPT TAKES QUASAR NAME FROM THE INPUT FILE TO IDENTIFY A CHAIN . 
FOR SDSS QSO'S IT IS quasar_name = files[j][5:]    
BUT FOR CRTS IT IS DIFFERENT  !!!  


"""
import numpy as np 

dir_choice = ['stars_CRTS_err_rms_chains/','stars_CRTS_err_w_chains/','stars_CRTS_analysis/' ]

err = 1  # or 1 

if err ==0 : 
    err_txt = 'err_rms'
else: 
    err_txt = 'err_w'
    
dir_in = dir_choice[err]
dir_out = dir_choice[2]


'''
NOTE : must make a chain_list_  ... .ls  file before running the program!
in stars_CRTS_chains/  run :
ls ch_*.dat_chain.dat > chain_list.ls

'''
filename = dir_in + 'chain.list'
files=np.loadtxt(filename,dtype=str)

# initialise storing vecfiles_rtors

sigma_l =  np.empty(0,dtype=float)
sigma_m =  np.empty(0,dtype=float)
sigma_h =  np.empty(0,dtype=float)
tau_l =  np.empty(0,dtype=float)
tau_m =  np.empty(0,dtype=float)
tau_h =  np.empty(0,dtype=float)
files_read = np.empty(0,dtype=str)

# load multiple chains 
for j in range(len(files)):   #
    fchain = dir_in+files[j]
    flatchain= np.genfromtxt(fchain)
    flatchain_whole = np.copy(flatchain)
    ndim = flatchain.shape[1]
    hpd = np.zeros((3,ndim))
    chain_len = flatchain.shape[0]
    
    """ 
    The chain consists of two columns if we are fitting tau and sigma, which are natural 
    logs of the true values - thus at the end we need to tae 
    pct1sig are points at which we probe the chain, scaled to the length of the chain 
    for the chain of length 5000, such values will be at 
    positions 800,2500, 4200.  
     """
    pct1sig = chain_len * np.array([0.16,0.50,0.84])  
    medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
    vars = ["sigma", "tau"]
    set_verbose=True  
    
    for i in xrange(ndim):
        vsort = np.sort(flatchain[:,i])  # sorts the array along either sigma or tau dimension 
        hpd[:,i] = vsort[medlowhig] # picks out values at the positions for the 
                                    # points at 15%, 50%, and 84% of the maximum posterior distribution
                        
    
    ln_sigma_lmh = hpd[:,0] 
    ln_tau_lmh = hpd[:,1]
    
    print 'HPD of sigma', np.exp(ln_sigma_lmh)
    print 'HPD of tau', np.exp(ln_tau_lmh)

    exp_sigma = np.exp(ln_sigma_lmh)
    exp_tau = np.exp(ln_tau_lmh)
    
    sigma_l = np.append(sigma_l,exp_sigma[0])
    sigma_m = np.append(sigma_m,exp_sigma[1])
    sigma_h = np.append(sigma_h,exp_sigma[2])
    tau_l = np.append(tau_l, exp_tau[0])
    tau_m = np.append(tau_m, exp_tau[1])
    tau_h = np.append(tau_h, exp_tau[2])
    star_name = files[j][3:-14]           
    print 'Star name', star_name
    files_read=np.append(files_read,star_name)
    
## save all the information to output file
print 'We have retrieved chain results data for ', len(sigma_l) ,' CRTS stars out of ', len(files)
fout = dir_out + 'javelin_CRTS_stars_'+err_txt+'_chain_results.txt'
DAT= np.column_stack((files_read, sigma_l, sigma_m, sigma_h, tau_l, tau_m, tau_h))
print 'We saved the results to ', fout 

# sort the DAT column accoring to QSO names 
newDAT=DAT[DAT[:,0].argsort()]

np.savetxt(fout,newDAT, delimiter=" ", fmt="%s")
