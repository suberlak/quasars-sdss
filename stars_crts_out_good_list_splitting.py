# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 12:28:54 2014

@author: suberlak

split the long lits into several shorter lists 

"""

import numpy as np
dir_choice=['stars_CRTS/0/','QSO_SDSS_JAV/', 'QSO_SDSS_JAV/MEAN_SUB/', 'qso_drw_upd/', 'qso_drw_long_LC/', 'qso_drw_medium_LC/']
d=dir_choice[5]

# for stars  :  ls=np.loadtxt(d+'out_good.list',dtype=str)

# for quasars : pre = 'u_band.ls'

# for drw
pre = 'drw_err3.list'
ls = np.loadtxt(d+pre,dtype=str)

length=201

j=0

for i in range(len(ls)/length + 1):
    filename = d+pre+'_pt'+str(i)
    print '\n',filename
    #print j
    
    if j+length < len(ls): 
        part_ls = ls[j:j+length]
        print 'saving part_ls from', j ,'to', j+length
        np.savetxt(filename,part_ls,fmt="%s")
    else:
        part_ls = ls[j:]
        print 'saving part_ls from', j ,'to the end, i.e. pos',len(ls)-1
        np.savetxt(filename,part_ls,fmt="%s")
    j += length    