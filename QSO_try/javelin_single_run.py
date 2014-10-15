# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 14:32:39 2014

@author: suberlak
"""

import numpy as np
from javelin.zylc import get_data 
from javelin.lcmodel import Cont_Model

names_raw=np.loadtxt('file.list',dtype=str)
chosen = names_raw[3]
lc = np.loadtxt(chosen)

if len(lc) > 50 : 
    filename= chosen
    data = get_data(filename,names=["Continuum"])
    cont=Cont_Model(data)
    chain_name = 'new4_'+filename[:18]+'_chain.dat'
    cont.do_mcmc(fchain=chain_name)
else: print 'Too short lightcurve to proceed' 


