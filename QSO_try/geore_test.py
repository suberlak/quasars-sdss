# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 15:26:22 2014

@author: suberlak
"""


# read in the files
import matplotlib.pyplot as plt
import numpy as np
from george import kernels
import george
from scipy.linalg import cho_solve 

#fname = '231408.02-011355.4.dat'   #Quasar lightcurve : fails at 363-th minor
#fname = '231413.86+002435.3.dat' # QSO fails at 94-th minor
#fname = '231425.13-010902.9.dat' # QSO fails at 49-th minor
#fname = '231428.54-011331.9.dat'   # strange result
#fname = '231416.82-011411.5.dat'  # strange result
#fname = '231435.45+010138.9.dat'   # strange result
#fname = '231436.60+003405.9.dat' # QSO fails at 72-th minor
#fname = '231444.74-011554.5.dat' # QSO fails at 62-th minor
#fname = '231445.49-010735.6.dat' # QSO fails at 203-th minor
#fname = '231446.49-000717.5.dat' # strange result
#fname='231453.63+004151.5.dat' # strange result
fname = 'DRW_test_new_1_short.txt' #DRW shortened mock lightcurve : works fine
#fname = 'con.dat'  # The javelin example continuum lightcurve 

data = np.loadtxt(fname)

t=mjd = data[:,0]
y=mags = data[:,1]
errs = data[:,2]

# Kernel setup 

kernel =( 0.2 )* kernels.ExpKernel(188.0) 
# This means that  we have two hyperparameters 
# ConstantKernel(0.2, ndim=1) * ExpKernel(188.0, ndim=1)

# Optimization

gp = george.GP(kernel, mean = np.mean(y))
gp.compute(t, yerr=errs)


print(gp.lnlikelihood(y))
print(gp.grad_lnlikelihood(y))

# Now it works all the way until here , when I added the line that 
# K[362,362] = K[361,361]

# And it breaks right here  
parameters,results = gp.optimize(t,y, yerr=errs)

print('Sigma is',parameters[0], ' and Tau is ', parameters[1] )

x = np.linspace(min(t), max(t))
mu, cov = gp.predict(y, x)
std = np.sqrt(np.diag(cov))

plt.plot(x,mu)

# Additional parts about lnlikelihood : 

# residuals
r = gp._check_dimensions(y)[gp.inds] - gp.mean(gp._x)

# compute lnlike

lnlike = gp._const - 0.5*np.dot(r.T, cho_solve(gp._factor, r))