# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 18:53:02 2014

@author: suberlak
"""

"""
Runs javelin on files on file.list.  Eg.

231413.86+002435.3.dat
231416.82-011411.5.dat
231425.13-010902.9.dat

Those are quasars, with eg. 
MJD          MAG  ERR
55828.27720 20.17 0.16 0
55828.28184 20.05 0.15 0

It saves the chains as 

new3_231413.86+002435.3_chain.dat

The follow-up program to read the chains is  
javelin_chain_retrieve_sigma_tau_save.py

"""
import numpy as np
from javelin.zylc import get_data 
from javelin.lcmodel import Cont_Model
import os 
import random

directory='../QSO/'
names_raw=np.loadtxt(directory+'file.list',dtype=str)

cond_notempty=np.empty_like(names_raw,dtype=bool)
for i in range(len(names_raw)):
    if os.stat(directory+names_raw[i])[6]!=0: 
        cond_notempty[i] = True
    else:
        cond_notempty[i] = False 

num_empty=np.sum(np.logical_not(cond_notempty))
num_notempty=np.sum(cond_notempty)

print 'Out of', len(names_raw), 'files, we have', num_notempty, 'not-empty files', \
'and therefore, ', num_empty, 'empty files'

# print 'Files that were empty:', names_raw[np.logical_not(cond_notempty)]

names_list=names_raw[cond_notempty]

"""
Restrict to only those lightcurves that have at least 50 observed points 
"""
cond_length = np.empty_like(names_list,dtype=bool)

for i in range(len(names_list)):  #
    test = np.loadtxt(directory+names_list[i])
    #print names_list[i], len(test)
    if len(test) > 50:
        cond_length[i] = True
    else:
        cond_length[i] = False

names_updated = names_list[cond_length]

print 'Test for length 2' 
for i in range(len(names_updated)):  #
    test = np.loadtxt(directory+names_updated[i])
    # print names_updated[i], len(test)

# So now the names_updated has only long files  

# random.shuffle(names_updated)


for i in range(len(names_updated)):  # 
    filename= directory+names_updated[i]
    print 'Working file', filename 
    data = get_data(filename,names=["Continuum"])
    cont=Cont_Model(data)
    start=len(directory)
    chain_name = 'new7_'+filename[start:start+18]+'_chain.dat'
    cont.do_mcmc(fchain=chain_name)
    
#cont.show_hist(bins=100)
