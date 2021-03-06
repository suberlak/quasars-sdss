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

UPDATE : 
03/18/2015
The difference between sf_load.py   and sf_load_NEW.py , is that instead
of making one big 'master' file, I decided to switch into making more 
smaller files, that may be much easier to handle 

03/23/2015
Added saving error_ij  (err_m_i and err_m_j  added in quadrature) , removed
storing <mag> and <mag_err> per LC, since this data is already stored as one 
number per LC in SDSS-CRTS cross-matched catalog 

NOTE : 
the program can be called   sf_load_NEW.py   inDir  outDir  outRoot  objType

"""

import os
import numpy as np 
import sys

# Read in the LC files : 
# if the parameters are provided , from the user input
args = sys.argv
if len(args) > 1 : 
    inDir = args[1]
    outDir = args[2]
    outRoot = args[3]
    objType = args[4]
    
    if objType == 'qso'  : start = 4 ; end = -8
    if objType == 'star' : start = 4 ; end = -4
    

# Read in the LC files  if there is no user input 
if len(args) == 1 :
    choice = 'star'

    if choice == 'star':
        inDir =  '../stars_CRTS_proc_err_w_good/' 
        #inDir =  './stars_CRTS_processed_err_w_good/'
        outRoot = 'SF_CRTS_stars_'
        start= 4
        end=  -8
        
    if choice == 'qso' :
        #inDir = '../QSO_CRTS_proc_err_w_good/'
        inDir = '../QSO_CRTS_processed_err_w/'  
        outRoot = 'SF_CRTS_quasars_'
        start= 4
        end=  -4
        
    if choice == 'test':
        inDir = './qso_TEST/'
        outRoot = 'SF_CRTS_quasars_TEST_'
        start= 4
        end = -4

    # regardless of choice, outDir stays the same 
    outDir = './sf_TRY/' 

if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 

# How many rows can the output file have ?
file_length = 1e6  # 1 000 000 rows 

# Load LC files 
inFiles = os.listdir(inDir)

# initialize storage arrays... 
avg_err_hold   = np.empty(0,dtype=float)
avg_mag_hold   = np.empty(0,dtype=float)
ID_hold        = np.empty(0,dtype=int)
tau_hold       = np.empty(0,dtype=float)
delflx_hold    = np.empty(0,dtype=float)
delflx2_hold   = np.empty(0,dtype=float)
delflxerr_hold = np.empty(0,dtype=float)

counter = 0
file_count = 0

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

    avg_err_hold   = np.append(avg_err_hold, avg_err_arr)
    avg_mag_hold   = np.append(avg_mag_hold, avg_mag_arr)
    ID_hold        = np.append(ID_hold, ID_arr)
    tau_hold       = np.append(tau_hold, tau)
    delflx_hold    = np.append(delflx_hold, delflx)
    delflx2_hold   = np.append(delflx2_hold, np.sqrt(delflx2))
    delflxerr_hold = np.append(delflxerr_hold,delflxerr)
    
    file_count += 1 
    
    # only save if reached the required max file length... 
    if (len(tau_hold) > file_length)  :
        
        # save the SF data and LC key identifiers 
        
        ##### FILE SAVING #######
        DATA = np.column_stack((delflx_hold,tau_hold, delflxerr_hold, ID_hold))    
        outfile = outDir + outRoot + str(counter) + '.txt'
        print 'Saving...', outfile 
        np.savetxt(outfile, DATA, delimiter =' ', fmt="%s")
        
        counter += 1 
        
        # erase the storage arrays, and reset the file counter
        avg_err_hold   = np.empty(0,dtype=float)
        avg_mag_hold   = np.empty(0,dtype=float)
        ID_hold        = np.empty(0,dtype=int)
        tau_hold       = np.empty(0,dtype=float)
        delflx_hold    = np.empty(0,dtype=float)
        delflx2_hold   = np.empty(0,dtype=float)
        delflxerr_hold = np.empty(0,dtype=float)
        file_count = 0
    
    # handle the case if we had a small number of lightcurves, that we didn't
    # even reach the maximum file size limit 
    if (file_count == len(inFiles)) and (len(tau_hold) < file_length ):
        print "We didn't reach the output SF file size limit... "
         # save the SF data and LC key identifiers 
        
        ##### FILE SAVING #######
        DATA = np.column_stack((delflx_hold,tau_hold, delflxerr_hold,  ID_hold))    
        outfile = outDir + outRoot + str(counter) + '.txt'
        print 'Saving all to ...', outfile 
        np.savetxt(outfile, DATA, delimiter =' ', fmt="%s")
        