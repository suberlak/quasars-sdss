# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 17:01:19 2014

@author: suberlak

Loading chains bit 

"""
import numpy as np 

files=np.loadtxt('chains_list_new7.ls',dtype=str)

# initialise storing vecfiles_rtors

sigma_l =  np.empty(0,dtype=float)
sigma_m =  np.empty(0,dtype=float)
sigma_h =  np.empty(0,dtype=float)
tau_l =  np.empty(0,dtype=float)
tau_m =  np.empty(0,dtype=float)
tau_h =  np.empty(0,dtype=float)
files_read = np.empty(0,dtype=str)

# load a single chain 

fchain = files[0]

flatchain= np.genfromtxt(fchain)
flatchain_whole = np.copy(flatchain)
ndim = flatchain.shape[1]
hpd = np.zeros((3,ndim))
chain_len = flatchain.shape[0]

""" 
The chain consists of two columns if we are fitting tau and sigma, which are natural 
logs of the true values - thus at the end we need to tae 
pct1sig are points at which we probe the chain, scaled to the length of the chain 
for the chain of length 5000, such values will be at 
positions 800,2500, 4200.  
 """
pct1sig = chain_len * np.array([0.16,0.50,0.84])  
medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
vars = ["sigma", "tau"]
set_verbose=True  

for i in xrange(ndim):
    vsort = np.sort(flatchain[:,i])  # sorts the array along either sigma or tau dimension 
    hpd[:,i] = vsort[medlowhig] # picks out values at the positions for the 
                                # points at 15%, 50%, and 84% of the maximum posterior distribution
    if set_verbose :
        print("HPD of %s"%vars[i])
        if i < 2 :
            # tau and sigma are stored as natural logs - other variables may not 
            print("low: %8.3f med %8.3f hig %8.3f"%tuple(np.exp(hpd[:,i])))
        else :
            print("low: %8.3f med %8.3f hig %8.3f"%tuple(hpd[:,i]))
                    

sigma_lmh = hpd[:,0] 
tau_lmh = hpd[:,1]

print 'HPD of sigma', np.exp(sigma_lmh)
print 'HPD of tau', np.exp(tau_lmh)