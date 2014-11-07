# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 15:36:41 2014

@author: suberlak
"""

import numpy as np

dirs = ['QSO_SDSS_JAV/','QSO_SDSS_chains/']
bands= ['u','g','r','i','z']
ch = 4
too_short_list = dirs[0]+bands[ch]+'_band_too_short_lc.ls'
all_chains_list = dirs[1]+'chain_list_'+bands[ch]+'.ls'
 
input_names_list = dirs[0]+bands[ch]+'_band.ls'
names=np.loadtxt(input_names_list,dtype=str)

'''
Check for those before removing the spurious files with 
xargs rm < chain_list_mistake_i.ls 
executed in  QSO_SDSS_chains/
'''

chains = np.loadtxt(all_chains_list,dtype=str)
short_names = np.loadtxt(too_short_list,dtype=str)
cond = np.zeros_like(chains, dtype=bool)  # will tell which chain comes from a too short quasar 

for i in range(len(chains)):
    qso_chain = chains[i][5:-14]   # cuts the  ch_i_ ... .txt_chain.dat  part
    for j in range(len(short_names)):
       qso_short = short_names[j][2:-4]
       if qso_chain == qso_short : 
           print qso_chain, qso_short
           cond[i] = True  # for those that are too short 
           
           
number = len(np.where(cond == True)[0])

print '\nWe have done a mistake of running javelin on too short chains for', number,' objects for  ' ,bands[ch], ' band' 

num_expected = len(names) - len(short_names)
print 'We have', len(chains), 'chains , but there should be only ', num_expected
print 'Check:', len(chains),'-',num_expected,'=', len(chains)-num_expected
filename = dirs[1]+ 'chain_list_mistake_'+bands[ch]+'.ls'
print '\nNames of chains that were run by accident are stored in a file:', filename

np.savetxt(filename, chains[cond], fmt="%s")

'''
Check that xargs removed all the unnecessary chains 
'''

new_chains_list = dirs[1]+'chain_list_new_'+bands[ch]+'.ls'
new_chains =np.loadtxt(new_chains_list,dtype=str)
cond = np.zeros_like(new_chains, dtype=bool)  # will tell which chain comes from a too short quasar 

for i in range(len(new_chains)):
    qso_chain = new_chains[i][5:-14]   # cuts the  ch_i_ ... .txt_chain.dat  part
    for j in range(len(short_names)):
       qso_short = short_names[j][2:-4]
       if qso_chain == qso_short : 
           print qso_chain, qso_short
           cond[i] = True  # for those that are too short 
           
           
number = len(np.where(cond == True)[0])

print '\nNow after removing all the unnecessary chains with xargs  there are ', number,' mistakenly made chains for  ' ,bands[ch], ' band' 

num_expected = len(names) - len(short_names)
print 'We have', len(new_chains), 'chains , and  there should be only ', num_expected
print 'Check:', len(new_chains),'-',num_expected,'=', len(new_chains)-num_expected
print 'Mission accomplished'