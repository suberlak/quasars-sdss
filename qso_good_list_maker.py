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

dirs = ['QSO_S82/','stars_CRTS_try/','stars_CRTS/0/','stars_SDSS_try/','stars_SDSS/0/']
d = dirs[2]
names_raw=np.loadtxt(d+'file.list',dtype=str)


# check which files are empty
# and count how many nonempty files meet that criterion

cond_notempty=np.empty_like(names_raw,dtype=bool)
for i in range(len(names_raw)):
    if os.stat(d+names_raw[i])[6]!=0: 
        cond_notempty[i] = True
    else:
        cond_notempty[i] = False 

num_empty=np.sum(np.logical_not(cond_notempty))
num_notempty=np.sum(cond_notempty)

#print 'Files that were empty:', names_raw[np.logical_not(cond_notempty)]

'''

Pick a subset of those files that were not empty , check whether they have at least 
two lines.

That way we do not perform any cuts yet, only ensuring stability of our code that would
calculate statistics. 

Selecting of eg. only those objects that are longer than 10 observations  happens
at the javelin- run stage .

'''

names=names_raw[cond_notempty]
length=len(names)

cond_multline=np.empty_like(names,dtype=bool)

for i in range(0,len(names)): 
    address=d+names[i]
    data=np.loadtxt(address)
    print 'Processing obj', i+1,' of ',length, 'from ', d, ': ' ,names[i]
    if data.size < 8.0 :
         cond_multline[i] = False
    else:
         cond_multline[i] = True

num_multline = np.sum(cond_multline)
num_singleline = np.sum(np.logical_not(cond_multline))

# print 'After  checking for single-line measurements, we have '
print '\nOut of', len(names_raw), 'files, we have', num_notempty, 'not-empty files', \
'and therefore, ', num_empty, 'empty files'

print '\nWe have',  num_multline, ' objects with more than two measurements, and',\
 num_singleline, ' objects with less than that two photometric measurements '
 
 
good_names = names[cond_multline]
# bad_names = 
filename = d+'out_good.list'
print '\nObjects ready to go for fitting and stats retrieval are listed in ',filename
np.savetxt(filename,good_names,fmt="%s")
# np.savetxt()