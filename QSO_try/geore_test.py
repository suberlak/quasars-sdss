# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 15:26:22 2014

@author: suberlak
"""


# read in the files

import numpy as np
from george import kernels
import george

#fname = '231408.02-011355.4.dat'   #Quasar lightcurve : fails at 363-th minor
#fname = 'DRW_test_new_1_short.txt' #DRW shortened mock lightcurve : works fine
fname = 'con.dat'  # The javelin example continuum lightcurve 

data = np.loadtxt(fname)

t=mjd = data[:,0]
y=mags = data[:,1]
errs = data[:,2]

# Kernel setup 

kernel =( 0.06**2.0 )* kernels.ExpKernel(700.0) + 0.06 
#print type(kernel)

# Optimization

gp = george.GP(kernel, mean = np.mean(y))
gp.compute(t, yerr=errs)


print(gp.lnlikelihood(y))
print(gp.grad_lnlikelihood(y))

# Now it works all the way until here , when I added the line that 
# K[362,362] = K[361,361]

# And it breaks right here  
parameters,results = gp.optimize(t,y)

x = np.linspace(min(t), max(t))
mu, cov = gp.predict(y, x)
std = np.sqrt(np.diag(cov))

plt.plot(x,mu)