# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 19:59:20 2014

@author: astronomy


Take SDSS lightcurves , ignore empty files and those less than 10 lines long, 
and calculate the median error (with np.percentile(error,50) ), 
and only add the object name to the list of  'good' ones,  
if the median(error) < 0.1 , and 90% of points are < 0.2 in error.

Thus created list is stored in QSO_SDSS_analysis  , and helps to decide 
which quasars should be plotted, when comparing the SDSS Chelsea results , 
and results of Javelin SDSS runs. 

"""
import numpy as np
import sys 
import os 


#args = sys.argv
#ch = int(args[1])

ch = 0
dir_in = ['QSO_SDSS_JAV/']
dir_out = ['QSO_SDSS_analysis/']
list_name = 'r_band.ls'  


list_file = dir_in[ch]+list_name
lc_names = np.loadtxt(list_file, dtype='str')

#
#  IGNORE EMPTY  LIGHTCURVES AND THOSE SHORTER THAN 10 LINES 
#

gi = np.ones_like(lc_names, dtype=bool) # good indices 

for i in range(len(lc_names)):
    if os.stat(dir_in[ch]+lc_names[i])[6]==0:  gi[i] = False  # bad if empty
    address=dir_in[ch]+lc_names[i]
    data=np.loadtxt(address)
    if data.size < 30.0 :  gi[i] = False  # bad if shorter than 10 lines 

num_too_short=np.sum(np.logical_not(gi))
num_long=np.sum(gi)

print 'Out of', len(lc_names), 'files, we have', num_long, 'files that were longer than 10 lines', \
'and therefore, ', num_too_short, 'almost empty files'

print 'Files that were too short :', lc_names[np.logical_not(gi)]


# substitute the name , choosing only good ones... 
lc_names = lc_names[gi]

print '\nPerforming calculations on files with more than one measurement...'

#
# CHECK WHICH FILES FULFILL THE ERROR CRITERIA 
# 

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