# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 18:33:40 2014

@author: suberlak

A program to create .cfg files,  as well as .csh files, each of which 
calls  javelin_drw_condor.py , in  a following fashion:

python  program.py  LC_name  dir_in  dir_chains  With_Prior?

for example : 

program =      /astro/users/suberlak/Desktop/SDSS/javelin_drw_condor.py
LC_name =      DRWtest_LC210_err1.dat
dir_in =       /astro/users/suberlak/Desktop/SDSS/qso_drw_upd/
dir_chains =   /astro/users/suberlak/Desktop/SDSS/qso_drw_chains/
With_Prior? = No   (or Yes )

Thus -  USE FULL ADDRESS ! 

This is because Condor jobs will be run from Condor directory, and 
stored elsewhere 

"""

import subprocess
import os
import stat 
import numpy as np 

root='/astro/users/suberlak/Desktop/SDSS/'
dir_in = ['QSO_CRTS_processed_err_rms/','QSO_CRTS_processed_err_w/']
dir_cond = ['QSO_CRTS_condor_err_rms/', 'QSO_CRTS_condor_err_w/']
dir_out = ['QSO_CRTS_err_rms_chains/','QSO_CRTS_err_w_chains/']

d = 1

filelist = dir_in[d] +'out.list'
names = np.loadtxt(filelist, dtype=str)
with_prior = ['Yes','No']


################################################################
# LOOP OVER NAMES,  TO MAKE A SEPARATE CONDOR RUN FOR EACH DRW #
################################################################

for i in range(2,len(names)):  #len(names)
     cshfile = root+dir_cond[d] + str(i) +'_'+ names[i] + '.csh'
     file = open(cshfile,"w")
     file.write('#!/bin/csh\n')
     file.write('source /astro/users/suberlak/.cshrc\n')
     s2 = 'echo python ' +root+'javelin_drw_condor.py'+' '+ \
         names[i] + ' ' + root+dir_in[d] +' '+ root+dir_out[d] +' '+ with_prior[d]  +'\n'   
     file.write(s2)
     s = 'python ' +root+'javelin_drw_condor.py'+' '+ \
         names[i] + ' ' + root+dir_in[d] +' '+ root+dir_out[d] +' '+ with_prior[d]     +'\n'   
     file.write(s)
     file.close()
     
     cfgfile = root+dir_cond[d] + str(i) +'_'+names[i] + '.cfg'
     file = open(cfgfile,"w")
     file.write('Notification = never\n\n')
     file.write('getenv = true\n\n')
     execs = 'Executable = ' + cshfile+'\n\n'
     file.write(execs)
     initial = 'Initialdir =' +root+dir_cond[d]+'\n\n'
     file.write(initial)
     file.write('Universe = vanilla\n\n')
     log = 'Log = '+root+dir_cond[d] + 'run_'+str(i)+'_log.txt\n\n' 
     out = 'Output = ' + root+dir_cond[d] + 'run_'+str(i)+'_out.txt\n\n' 
     err = 'Error =' +  root+dir_cond[d] + 'run_'+str(i)+'_err.txt\n\n' 
     file.write(log)
     file.write(out)
     file.write(err)
     file.write('Queue\n') 
     file.close()
     
     # source my cshrc file to enable Javelin run 
     # subprocess.check_call(['source', ' ~/.cshrc' ])
     
     # Change file permissions of the csh file 
     st = os.stat(cshfile)
     os.chmod(cshfile, st.st_mode | 0111)
     
     # Submit the fitting job to Condor 
     subprocess.check_call(['condor_submit', cfgfile])
     