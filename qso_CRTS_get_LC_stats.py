# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 21:39:30 2014

@author: astronomy

Generate stats for  CRTS stars, to figure out why the log(tau) vs log(sigma)
plot has such a sharp cut at log(sigma ) = 0 .

Read in the list of processed stars, and read in each lightcurve. 
Calculate mjd span, mag_scatter,  <magnitude>, <error>, N_lines.. 
Save including the LC name  ! 

"""
import numpy as np
import sys 


args = sys.argv
ch = int(args[1])

#ch = 0
dir_in = ['QSO_CRTS_processed_err_w/', 'stars_CRTS_processed_err_w/']
dir_out = ['QSO_CRTS_analysis/', 'stars_CRTS_analysis/']
list_name = 'out.list'  # made automatically by   qso_crts_preprocessing.py, or stars_crts_preprocessing.py


list_file = dir_in[ch]+list_name
lc_names = np.loadtxt(list_file, dtype='str')

mjd_span = np.zeros_like(lc_names, dtype=float)
mag_rms= np.zeros_like(lc_names, dtype=float)
mag_mean=np.zeros_like(lc_names, dtype=float)
err_mean=np.zeros_like(lc_names, dtype=float)
N_lines=np.zeros_like(lc_names, dtype=float)

for i in range(len(lc_names)): #
    data = np.loadtxt(dir_in[ch]+lc_names[i], dtype='float')
    mjd = data[:,0]
    mag = data[:,1]
    err = data[:,2]
    mjd_span[i] = int(mjd.max()-mjd.min())
    mag_mean[i] = np.mean(mag)
    mag_rms[i] = np.sqrt(np.mean(np.power( mag_mean[i] - mag, 2.0)))  
    err_mean[i] = np.mean(err)
    N_lines[i] = len(mag)
    print 'Getting stats...', str((float(i) / len(lc_names))*100.0)[:5],'%'
    
    
DAT = np.column_stack((lc_names, mjd_span, mag_rms, mag_mean, err_mean, N_lines))
fout = dir_out[ch] + 'LC_stats.txt'
np.savetxt(fout, DAT, delimiter=" ", fmt="%s")

print 'Saved stats such as  lc_names,  mjd_span, mag_rms, mag_mean, err_mean, N_lines, for ',\
 len(lc_names), ' objects to ', fout