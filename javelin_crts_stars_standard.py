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

dir_choice=['stars_CRTS_try/','stars_CRTS/0/','stars_CRTS_chains/']

dir_input=dir_choice[1]
dir_output=dir_choice[2]

pre = 'out_good.list_pt10'

names=np.loadtxt(dir_input+pre,dtype=str)



"""
Restrict to only those lightcurves that have at least 10 observed points 

BEGINNING OF CODE 

"""
cond_length = np.empty_like(names,dtype=bool)
min_length = 10

for i in range(len(names)):  #
    test = np.loadtxt(dir_input+names[i])
    print i, names[i], len(test)
    if len(test) > min_length:
        cond_length[i] = True
    else:
        cond_length[i] = False

names_upd = names[cond_length]

#print '\nTest for length 2' 
#for i in range(len(names_upd)):  #
#    test = np.loadtxt(dir_input+names_upd[i])
#    print names_upd[i], len(test)


too_short_num = len(np.where(cond_length == False)[0])

print '\nWe have', too_short_num,' stars whose lightcurves were shorter than ', min_length
print '\nThus of all analysed ', len(names), 'stars, ', len(names_upd),' qualify for javelin fitting.' 


file_short_lc = dir_input+pre+'_stars_too_short.ls'
cond_short = -cond_length
files_list = names[cond_short]
np.savetxt(file_short_lc ,files_list,  fmt="%s")

print 'Their names were saved to a file', file_short_lc
"""
END OF CODE 
"""


for i in range(len(names_upd)):  # len(names_upd)
    filename=dir_input+names_upd[i]
    print '\nWorking file', filename, i+1, ' out of ', len(names_upd), ' from ', pre
    data = get_data(filename,names=["Continuum"])
    cont=Cont_Model(data)
    start=len(dir_input)
    chain_name = dir_output+'ch_'+filename[start:]+'_chain.dat'
    print chain_name
    cont.do_mcmc(fchain=chain_name)
    

print '\nJavelin performed', i+1, 'fittings, and it should have done ', len(names_upd)