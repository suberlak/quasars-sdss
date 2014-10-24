# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 17:18:35 2014

@author: suberlak

For the variability analysis of the QSO's some objects have only one measurement, 
or no measurement at all, which is useless for our analysis. Thus restrict the 
further code to only those objects that have more than two measurements. 

This is particularly essential for   qso_lc_stats.py   code, which is used to pull 
all the stats I could think of for every quasar in the sample. 
So removing from the sample those that have two or less measurements seems 
sensible.  

Input : file.list  with the list of ALL quasars in  a directory  

Output : list of quasars with more than two measurements. 

Note: this code is mostly derived from  qso_var.py , taking the file - checking 
part. 
"""


import numpy as np
np.seterr(invalid='raise')
import warnings 
import os 

directory='QSO_S82/'
names_raw=np.loadtxt(directory+'file.list',dtype=str)


# check which files are empty
# and count how many nonempty files meet that criterion

cond_notempty=np.empty_like(names_raw,dtype=bool)
for i in range(len(names_raw)):
    if os.stat(directory+names_raw[i])[6]!=0: 
        cond_notempty[i] = True
    else:
        cond_notempty[i] = False 

num_empty=np.sum(np.logical_not(cond_notempty))
num_notempty=np.sum(cond_notempty)

#print 'Files that were empty:', names_raw[np.logical_not(cond_notempty)]

# pick out only thosre 
names=names_raw[cond_notempty]

# of all notempty files 
# check if there are any files with only one line of measurement



for i in range(0,len(names)): 
    address=directory+names[i]
    data=np.loadtxt(address)
    print i
    if data.size < 8.0 :
         cond_multline[i] = False
    else:
         cond_multline[i] = True

num_multline = np.sum(cond_multline)
num_singleline = np.sum(np.logical_not(cond_multline))

# print 'After  checking for single-line measurements, we have '
print 'Out of', len(names_raw), 'files, we have', num_notempty, 'not-empty files', \
'and therefore, ', num_empty, 'empty files'

print 'We have',  num_multline, ' objects with more than two measurements, and',\
 num_singleline, ' objects with less than that two measurements of QSO variability  '
 
 
good_names = names[cond_multline]
# bad_names = 
filename = directory+'out_good.list'
np.savetxt(filename,good_names,fmt="%s")
# np.savetxt()