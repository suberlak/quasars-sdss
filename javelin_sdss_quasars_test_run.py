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

dir_choice=['QSO_SDSS_JAV/','QSO_SDSS_chains/test/', 'QSO_SDSS_JAV/','QSO_SDSS_chains/']

dir_input=dir_choice[0]
dir_output=dir_choice[1]
band= 'u_band.ls_pt6'



names=np.loadtxt(dir_input+band,dtype=str)


"""
Restrict to only those lightcurves that have at least 10 observed points 

BEGINNING OF CODE 

"""
cond_length = np.empty_like(names,dtype=bool)
min_length = 10

for i in range(len(names)):  # 
    test = np.loadtxt(dir_input+names[i])
    #print names[i], len(test)
    if len(test) > min_length:
        cond_length[i] = True
    else:
        cond_length[i] = False

names_upd = names[cond_length]

#print 'Test for length 2' 
#for i in range(len(names_upd)):  #
#    test = np.loadtxt(dir_input+names_upd[i])
    #print names_upd[i], len(test)


too_short_num = len(np.where(cond_length == False)[0])



if too_short_num > 0 : 
    print '\nWe have', too_short_num,' quasars whose lightcurves were shorter than ', min_length
    
    file_short_lc = dir_input+band+'.ls_test_too_short_lc.ls'
    cond_short = -cond_length
    files_list = names[cond_short]
    np.savetxt(file_short_lc ,files_list,  fmt="%s")
    
    print 'Their names were saved to a file', file_short_lc
else:
    print '\nAll quasars in this list are longer than', min_length
    

"""
END OF CODE 
"""

"""
How I check where to resume running it, from the file where it crashed:
crashed_file = i_2104983.txt

np.where(names == 'i_2104983.txt')
names[3162]
names[3161] :  the file that was fine,  just before the one where it crashed 
np.where(names_upd == 'i_2104791.txt')

-->  3119 --> start at location  3119    

"""
#

#
for i in range(0,len(names_upd)):  # 
    filename=dir_input+names_upd[i]
    print '\nWorking file', filename, i+1, ' out of ', len(names_upd) , ' from ', band
    data = get_data(filename,names=["Continuum"])
    cont=Cont_Model(data)
    start=len(dir_input)
    chain_name = dir_output+'ch_'+filename[start:]+'_chain.dat'
    print chain_name
    cont.do_mcmc(set_prior=False, fchain=chain_name)
    
    

print '\nJavelin performed', i+1, 'fittings, and it should have done ', len(names_upd)