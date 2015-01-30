# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:53:08 2015

@author: astronomy

Does the same as sf_yumi.py  EXCEPT  fitting anything - it just saves 
the mag_difference_squared,  the mjd_difference, and the mag_diff_error 
 ( made by adding the corresponding errors in quadrature )

Based on sf_yumi,  I read-in files  in the directory, and
plot the SF  : 
- for each file, I calculate the delta_mag vs delta_time  for each pair of 
  measurements
- From all those files together I bin the data into bins of with time_bin , 
  and for each bin I calculate rms from the <delta_mag>: this is SF
- I plot the SF vs delta_time 
"""

import os
import numpy as np 
import sys

# Read in the LC files 
args = sys.argv
if len(args) > 1 : 
    inDir = args[1]
    outDir = args[2]
    outfile = args[3]
    objType = args[4]
    
    if objType == 'qso' :  start = 4 ; end = -8
    if objType == 'star' :  start=4 ; end = -4
    


if len(args) == 1 :
    qso_or_star = 'qso'

    if qso_or_star == 'star':
        inDir =  './stars_CRTS_proc_err_w_good_TRY/' 
        outfile = 'SF_CRTS_stars_master.txt'
        start= 4
        end=  -8
        
    if qso_or_star == 'qso' :
        inDir = '../QSO_CRTS_proc_err_w_good/'
        outfile = 'SF_CRTS_quasars_master.txt'
        start= 4
        end=  -4
        
    if qso_or_star == 'test':
        inDir = './qso_TEST/'
        outfile = 'SF_CRTS_quasars_TEST.txt'
        start= 4
        end = -4

# regardless of choice, outDir stays the same 
    outDir = './sf_TRY/' 

if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 


inFiles = os.listdir(inDir)

# initialize big storage arrays... 
avg_err_hold   = np.empty(0,dtype=float)
avg_mag_hold   = np.empty(0,dtype=float)
ID_hold        = np.empty(0,dtype=int)
tau_hold       = np.empty(0,dtype=float)
delflx_hold    = np.empty(0,dtype=float)
delflx2_hold   = np.empty(0,dtype=float)
delflxerr_hold = np.empty(0,dtype=float)

for i in range(len(inFiles)):
    file = str(inFiles[i])
    print '\nFile ',i, 'out of',  len(inFiles)
    # load the mjd, flux, and flux_error  from the read file 
    mjd,flx4,err = np.loadtxt(inDir+'%s' % (file),usecols=(0,1,2),unpack=True)
    
    # check for any nan's...
    add = np.where(np.isnan(flx4) == True)[0]
   
    # Delete all rows with nan's - there should be none in files I provide.. 
    mjd = np.delete(mjd,add); flx4 = np.delete(flx4,add); err = np.delete(err,add)
        
    # Sort according to mjd's 
    ind = np.argsort(mjd)
    mjd = mjd[ind]; flx4 = flx4[ind]; err = err[ind]
    
     # Calculate tau, mag difference and err difference for each pair (k,j), 
     # where tau=t_j-t_k (j>k) 
    
    delflx2  = []; delflx = [];  delflxerr = []; tau = []
    for j in range(len(mjd)-1):
        for k in range(j+1):     
            tau.append(mjd[j+1]-mjd[k])  # j from 1 and k<j
            noise2 = err[k]**2+err[j+1]**2 
            delflx.append((flx4[k]-flx4[j+1]))
            delflx2.append((flx4[k]-flx4[j+1])**2.0)
            delflxerr.append((noise2)**0.5)  
            
    # Change lists to arrays..         
    tau = np.array(tau); delflx = np.array(delflx)
    delflx2 = np.array(delflx2);   delflxerr = np.array(delflxerr) 
    
    # sort according to tau
    int0 = np.argsort(tau)
    tau = tau[int0]; delflx = delflx[int0]
    delflx2 = delflx2[int0]; delflxerr = delflxerr[int0]
    
    # Make arrays that repeat the LC stats 
     
    ID_arr = np.array(len(tau)*[file[start:end]])
    avg_mag_arr = np.ones_like(tau, dtype=float) * np.mean(flx4)
    avg_err_arr = np.ones_like(tau,dtype=float) * np.mean(err)    
     # I assume that the filename structure is 
                          #  out_IDIDIDIDID.dat.txt   
    
    ###  NO FITTING HERE !!   ###  
     
    # save the results of calculations from one file ....
    # the code below does exactly the same as Yumi's WriteColumns routine,
    # but at a fraction of time...
     
    #np.savetxt(outDir+'tau%s.dat' % (file) ,np.column_stack((tau)), delimiter="\n", fmt="%s")
    #np.savetxt(outDir+'sf%s.dat' % (file) ,np.column_stack((np.sqrt(delflx2))), delimiter="\n", fmt="%s")
    
    avg_err_hold   = np.append(avg_err_hold, avg_err_arr)
    avg_mag_hold   = np.append(avg_mag_hold, avg_mag_arr)
    ID_hold        = np.append(ID_hold, ID_arr)
    tau_hold       = np.append(tau_hold, tau)
    delflx_hold    = np.append(delflx_hold, delflx)
    delflx2_hold   = np.append(delflx2_hold, np.sqrt(delflx2))
    delflxerr_hold = np.append(delflxerr_hold,delflxerr)
    
DATA = np.column_stack((delflx_hold,tau_hold, avg_mag_hold, avg_err_hold, ID_hold))    
np.savetxt(outDir +     outfile, DATA, delimiter =' ', fmt="%s")

print 'Saved ', outDir ,     outfile