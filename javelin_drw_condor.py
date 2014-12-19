# -*- coding: utf-8 -*-
"""
Created on Sat Oct 18 21:16:19 2014

@author: suberlak


Runs a single  lightcurve with javelin - suited for being submitted as condor job,
to take as an argument the name of the lightcurve, which would come from the wrapper.

The condor_wrapper.py can also store lightcurve length or anything else, if necessary. 

javelin_drw_condor.py  DRWtest_lc_name   dir_in    dir_out  no

--> the last argument specifies whether we want to run it with or without a prior 
--> list of arguments sys.argv includes by default as the first entry a name of the 
    python script  run 
"""
from time import clock
from javelin.zylc import get_data 
from javelin.lcmodel import Cont_Model
import sys 

###############
#    INPUT    #
###############

args = sys.argv
lc_name = args[1]
dir_input=args[2]  # dir with LC's 
dir_output=args[3] # dir with chains 
prior_on=args[4]

print 'Fitting with Javelin LC', lc_name , ' from ', dir_input, 'saving the chain to',\
 dir_output, 'javelin prior = ', prior_on

################
# TIME KEEPING #
################

def bench(secs):
  print 'Time elapsed: ' + str(secs) + 's'
  
#################
# JAVELIN  RUN  #
#################
start_time = clock()

filename = dir_input + lc_name    
print '\nWorking file', filename
data = get_data(filename,names=["Continuum"])
cont=Cont_Model(data)
start=len(dir_input)
chain_name = dir_output+'ch_'+filename[len(dir_input):-4]+'_chain.dat'  # -4 cuts off the .txt
print 'Chain name will be ', chain_name

if prior_on == ('Yes' or 'yes') : 
    print  'FITTING WITH JAVELIN  PRIOR'
    cont.do_mcmc(set_prior=True,fchain=chain_name)
if prior_on == ('No' or 'no') :
    print 'NO JAVELIN  PRIOR'
    cont.do_mcmc(set_prior=False,fchain=chain_name)
    
end_time = clock()
delta_t = end_time-start_time
bench(delta_t)

