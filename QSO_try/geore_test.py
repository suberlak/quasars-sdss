# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 15:26:22 2014

@author: suberlak
"""


# read in the files

import numpy as np
from george import kernels
import george

fname = '231408.02-011355.4.dat'
data = np.loadtxt(fname)

t=mjd = data[:,0]
y=mags = data[:,1]
errs = data[:,2]

# Kernel setup 

kernel = 17.0**2 * kernels.ExpSquaredKernel(57.0**2)
print type(kernel)

# Optimization

gp = george.GP(kernel, mean = np.mean(y))
gp.kernel(gp._x[:, None], gp._x[None, :])
#gp = george.HODLRGP(kernel,mean=np.mean(y))
gp.compute(t,yerr=errs)

print(gp.lnlikelihood(y))
#print(gp.grad_lnlikelihood(y))

parameters,results = gp.optimize(t,y)