# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 17:45:14 2015

@author: suberlak

A program to create .cfg files,  as well as .csh files, each of which 
calls  mcmc.py , in  a following fashion:

python  program.py bin_name  dir_in  dir_chains  

for example : 

program =      /astro/users/suberlak/Desktop/SDSS/SF_plotting/mcmc.py
LC_name =      DRWtest_LC210_err1.dat
dir_in =       /astro/users/suberlak/Desktop/SDSS/SF_plotting/QSO_bins_17-19_NEW_code/
dir_chains =   /astro/users/suberlak/Desktop/SDSS/SF_plotting/QSO_bins_17-19_chains/


Thus -  USE FULL ADDRESS ! 

This is because Condor jobs will be run from Condor directory, and 
stored elsewhere 

"""

import subprocess
import os


root='/astro/users/suberlak/Desktop/SDSS/SF_plotting/'
dir_in ='StarR_bins_18.5-19/'# 'QSO_bins_18.5-19_corr/' #'StarR_bins_18.5-19_corr/' #'StarB_bins_18.5-19_corr/'
#'QSO_bins_18.5-19/' # 'StarR_bins_18.5-19/' #'StarB_bins_18.5-19/'

# 'QSO_bins_17-19_corr/'   # 'StarR_bins_17-19_corr/' # 'StarB_bins_17-19_corr/'#
#'StarR_bins_17-19/' #'StarB_bins_17-19/' #, 'QSO_bins_17-19_NEW_code/'

# check if condor dir exists, and if not, make it 
dir_cond = dir_in[:-1]+'_condor/'
if not os.path.exists(dir_cond):
    os.makedirs(dir_cond)

# check if output chain directory exists, and if not, make it 
dir_out = dir_in[:-1]+'_chains/'
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
            

b = 0
names = os.listdir(dir_in)  

################################################################
# LOOP OVER NAMES,  TO MAKE A SEPARATE CONDOR RUN FOR EACH DRW #
################################################################

for i in range(len(names)):  #len(names)
     cshfile = root+dir_cond + names[i] + '.csh'
     file = open(cshfile,"w")
     file.write('#!/bin/bash\n')
     file.write('source /astro/users/suberlak/.bash_profile\n')
     s2 = 'echo python ' +root+'mcmc.py'+' '+ \
         names[i] + ' ' + root+dir_in +' '+ root+dir_out +' \n'   
     file.write(s2)
     s = 'python ' +root+'mcmc.py'+' '+ \
         names[i] + ' ' + root+dir_in +' '+ root+dir_out +' \n'   
     file.write(s)
     file.close()
     
     cfgfile = root+dir_cond + names[i] + '.cfg'
     file = open(cfgfile,"w")
     file.write('Notification = never\n\n')
     file.write('getenv = true\n\n')
     execs = 'Executable = ' + cshfile+'\n\n'
     file.write(execs)
     initial = 'Initialdir =' +root+dir_cond+'\n\n'
     file.write(initial)
     file.write('Universe = vanilla\n\n')
     log = 'Log = '+root+dir_cond + 'run_'+str(i)+'_log.txt\n\n' 
     out = 'Output = ' + root+dir_cond + 'run_'+str(i)+'_out.txt\n\n' 
     err = 'Error =' +  root+dir_cond + 'run_'+str(i)+'_err.txt\n\n' 
     rel = 'periodic_release = ((JobStatus == 5) && (time() - EnteredCurrentStatus) >  1800)\n\n'
     file.write(log)
     file.write(out)
     file.write(err)
     file.write(rel)
     file.write('Queue\n') 
     file.close()
     
     # source my cshrc file to enable Javelin run 
     # subprocess.check_call(['source', ' ~/.cshrc' ])
     
     # Change file permissions of the csh file 
     st = os.stat(cshfile)
     os.chmod(cshfile, st.st_mode | 0111)
     
     # Submit the mcmc job to Condor 
     subprocess.check_call(['condor_submit', cfgfile])
     