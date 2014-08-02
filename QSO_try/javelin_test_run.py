# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 18:53:02 2014

@author: suberlak
"""
import numpy as np
from javelin.zylc import get_data 
from javelin.lcmodel import Cont_Model
import os 

names_raw=np.loadtxt('file.list',dtype=str)

cond_notempty=np.empty_like(names_raw,dtype=bool)
for i in range(len(names_raw)):
    if os.stat(names_raw[i])[6]!=0: 
        cond_notempty[i] = True
    else:
        cond_notempty[i] = False 

num_empty=np.sum(np.logical_not(cond_notempty))
num_notempty=np.sum(cond_notempty)

print 'Out of', len(names_raw), 'files, we have', num_notempty, 'not-empty files', \
'and therefore, ', num_empty, 'empty files'

# print 'Files that were empty:', names_raw[np.logical_not(cond_notempty)]

names=names_raw[cond_notempty]

#filename= names[0]
filename = '231408.02-011355.4_extended.dat'
print 'Working file', filename 


data = get_data(filename,names=["Continuum"])
cont=Cont_Model(data)

cont.do_mcmc(fchain="new_chain_for_extended_curve.dat")
cont.show_hist(bins=100)