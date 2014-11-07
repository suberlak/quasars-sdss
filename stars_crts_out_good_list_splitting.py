# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 12:28:54 2014

@author: suberlak

split the long lits into several shorter lists 

"""

import numpy as np
dir_choice=['stars_CRTS/0/','QSO_SDSS_JAV/']
d=dir_choice[1]

# for stars  :  ls=np.loadtxt(d+'out_good.list',dtype=str)

# for quasars : 
band = 'u_band.ls'
ls = np.loadtxt(d+band,dtype=str)

length=1000

j=0

for i in range(len(ls)/length + 1):
    filename = d+band+'_pt'+str(i)
    print filename
    print j
    
    if j+length < len(ls): 
        part_ls = ls[j:j+length]
        print 'saving part_ls from', j ,'to', j+length
        np.savetxt(filename,part_ls,fmt="%s")
    else:
        part_ls = ls[j:]
        print 'saving part_ls from', j ,'to the end, i.e. pos',len(ls)-1
        np.savetxt(filename,part_ls,fmt="%s")
    j += length    