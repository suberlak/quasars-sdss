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

band= 'z_band'
names=np.loadtxt(dir_input+band+'.ls',dtype=str)

    
    
"""
Restrict to only those lightcurves that have at least 50 observed points 

BEGINNING OF CODE 

"""
cond_length = np.empty_like(names,dtype=bool)
min_length = 10

for i in range(len(names)):  #
    test = np.loadtxt(dir_input+names[i])
    print names[i], len(test)
    if len(test) > min_length:
        cond_length[i] = True
    else:
        cond_length[i] = False

names_upd = names[cond_length]

print 'Test for length 2' 
for i in range(len(names_upd)):  #
    test = np.loadtxt(dir_input+names_upd[i])
    print names_upd[i], len(test)


too_short_num = len(np.where(cond_length == False)[0])

print 'We have', too_short_num,' quasars whose lightcurves were shorter than ', min_length

file_short_lc = dir_input+band+'_too_short_lc.ls'
cond_short = -cond_length
files_list = names[cond_short]
np.savetxt(file_short_lc ,files_list,  fmt="%s")

print 'Their names were saved to a file', file_short_lc
"""
END OF CODE 
"""




for i in range(len(names_upd) ):  # 
    filename=dir_input+names_upd[i]
    print '\nWorking file', filename, i+1, ' out of ', len(names_upd) 
    data = get_data(filename,names=["Continuum"])
    cont=Cont_Model(data)
    start=len(dir_input)
    chain_name = dir_output+'ch_'+filename[start:]+'_chain.dat'
    print chain_name
    cont.do_mcmc(fchain=chain_name)
    
    
    
