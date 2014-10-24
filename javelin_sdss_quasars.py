# -*- coding: utf-8 -*-
"""
Created on Sat Oct 18 21:16:19 2014

@author: suberlak

Runs javelin on  SDSS Quasars, a selected band 

(filenames start as u_ , g_, etc.)

from QSO_SDSS_JAV

"""
import numpy as np
from javelin.zylc import get_data 
from javelin.lcmodel import Cont_Model

dir_choice=['QSO_try2/', 'QSO_S82/', 'QSO_SDSS_JAV/','QSO_SDSS_chains/']

dir_input=dir_choice[2]
dir_output=dir_choice[3]
names=np.loadtxt(dir_input+'u_band.ls',dtype=str)

for i in range(2):  # len(names)
    filename=dir_input+names[i]
    print '\nWorking file', filename, i+1, ' out of ', len(names) 
    data = get_data(filename,names=["Continuum"])
    cont=Cont_Model(data)
    start=len(dir_input)
    chain_name = dir_output+'ch_'+filename[start:]+'_chain.dat'
    print chain_name
    cont.do_mcmc(fchain=chain_name)
    
    
    
