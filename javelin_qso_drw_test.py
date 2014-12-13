# -*- coding: utf-8 -*-
"""
Created on Sat Oct 18 21:16:19 2014

@author: suberlak

Runs javelin on  SDSS Quasars, a selected band 

(filenames start as u_ , g_, etc.)

from QSO_SDSS_JAV

"""
from time import clock
import numpy as np
from javelin.zylc import get_data 
from javelin.lcmodel import Cont_Model

dir_choice=['qso_drw_medium_LC/','qso_drw_medium_LC_chains/no_prior/',\
'qso_drw_medium_LC_chains/with_prior/',\
'qso_drw_medium_LC_chains/no_prior/','qso_drw_upd/',\
'qso_drw_long_LC/', 'qso_drw_long_LC_chains/',\
'qso_drw_chains/', 'qso_drw_chains/no_prior/', 'qso_drw_long_LC/conv_time_test/',\
'qso_drw_long_LC_chains/conv_time_test/']

dir_input=dir_choice[0]  # dir with LC's 
dir_output=dir_choice[2] # dir with chains 

pre = 'drw_err2.list_pt4'
#pre = 'drw_err.list'

names=np.loadtxt(dir_input+pre,dtype=str)

LC_length = np.empty(0,dtype=float)

for i in range(len(names)):
    filename  =dir_input+names[i]
    D = np.loadtxt(filename, dtype='str')
    # print len(D)
    LC_length = np.append(LC_length,len(D))
    
"""
No need to do any restrictions, since I know that all lc's are longer then 10 pts

"""
def bench(secs):
  print 'Time elapsed: ' + str(secs) + 's'
  
time_lasted = np.empty(0,dtype=float)
for i in range(187,len(names)): 
    start_time = clock()
    
    filename=dir_input+names[i]
    print '\nWorking file', filename, i+1, ' out of ', len(names), ' from ', pre
    data = get_data(filename,names=["Continuum"])
    cont=Cont_Model(data)
    start=len(dir_input)
    chain_name = dir_output+'ch_'+filename[start:]+'_chain.dat'
    print chain_name,  'WITH JAVELIN  PRIOR'
    cont.do_mcmc(set_prior=True,fchain=chain_name)
#    
    end_time = clock()
    delta_t = end_time-start_time
    bench(delta_t)
    time_lasted = np.append(time_lasted, delta_t)
    
print '\nJavelin performed', i+1, 'fittings, and it should have done ', len(names)
cols = np.column_stack((LC_length,time_lasted))
out=dir_output + 'lc_length_vs_fitting_time_'+pre+'.txt'
np.savetxt(out,cols,fmt='%s')
