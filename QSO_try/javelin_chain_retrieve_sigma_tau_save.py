# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 12:20:35 2014

@author: suberlak

This program uses javelin to run! 


"""

import numpy as np 

from javelin.zylc import get_data 
from javelin.lcmodel import Cont_Model


files=np.loadtxt('chains_list_new7.ls',dtype=str)

# initialise storing vecfiles_rtors

sigma_l =  np.empty(0,dtype=float)
sigma_m =  np.empty(0,dtype=float)
sigma_h =  np.empty(0,dtype=float)
tau_l =  np.empty(0,dtype=float)
tau_m =  np.empty(0,dtype=float)
tau_h =  np.empty(0,dtype=float)
files_read = np.empty(0,dtype=str)

# load all markov chains from the file list 

for i in range(len(files)):
    cont.load_chain(files[i])
    cont.get_hpd()
    
    sigma_lmh = cont.hpd[:,0] 
    tau_lmh = cont.hpd[:,1]
    
    
    
    print 'HPD of sigma', np.exp(sigma_lmh)
    print 'HPD of tau', np.exp(tau_lmh)
    
    exp_sigma = np.exp(sigma_lmh)
    exp_tau = np.exp(tau_lmh)
    
    sigma_l = np.append(sigma_l,exp_sigma[0])
    sigma_m = np.append(sigma_m,exp_sigma[1])
    sigma_h = np.append(sigma_h,exp_sigma[2])
    tau_l = np.append(tau_l, exp_tau[0])
    tau_m = np.append(tau_m, exp_tau[1])
    tau_h = np.append(tau_h, exp_tau[2])
    quasar_name = files[i][5:23]
    files_read=np.append(files_read,quasar_name)
#    
## save all the information to output file
#
fout = 'new7_partial_qso_sigma_tau_stack.txt'
DAT= np.column_stack((files_read, sigma_l, sigma_m, sigma_h, tau_l, tau_m, tau_h))

# sort the DAT column accoring to QSO names 
newDAT=DAT[DAT[:,0].argsort()]

np.savetxt(fout,newDAT, delimiter=" ", fmt="%s")


    
#np.savetxt(fout2, np.column_stack((mjd,mag,error,Nobs,\
#    chi2)),fmt='%11.4f')

# solved with http://stackoverflow.com/questions/16621351/how-to-use-python-numpy-savetxt-to-write-strings-and-float-number-to-an-ascii-fi  