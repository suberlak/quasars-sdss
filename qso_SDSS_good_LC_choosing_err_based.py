# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 19:59:20 2014

@author: astronomy


Take CRTS lightcurves  (both stars and quasars) , and calculate the median error,
and take np.percentile(error), and only add the object name to the list of 
'good' ones,  if the median(error) < 0.1 , and 90% of points are < 0.2 in error.

"""
import numpy as np
import sys 

# note  : the one on stars has to be run on drizzle , because I don't locally 
# store all the CRTS stars!  

#args = sys.argv
#ch = int(args[1])

ch = 0
dir_in = ['QSO_SDSS_JAV/']
dir_out = ['QSO_SDSS_analysis/']
list_name = 'r_band.ls'  # made automatically by   qso_crts_preprocessing.py, or stars_crts_preprocessing.py


list_file = dir_in[ch]+list_name
lc_names = np.loadtxt(list_file, dtype='str')
good_lc = np.empty(0,dtype=str)
for name in lc_names :
    data = np.loadtxt(dir_in[ch]+name, dtype='float')
    err = data[:,2]
    if  (np.median(err) < 0.1 and np.percentile(err,90)) < 0.2 : 
        good_lc = np.append(good_lc,name)
        
print 'We accepted ', len(good_lc), 'out of ', len(lc_names)
print 'Thus the percentage of accepted LCs based on error criterion:',\
str(( len(good_lc) / float(len(lc_names)) ) * 100.0)[:5], '%'

#DAT = np.column_stack((good_lc))
np.savetxt(dir_out[ch]+'good_err_LC.txt',good_lc, delimiter=" ", fmt="%s")