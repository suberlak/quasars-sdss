# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 16:06:07 2014

@author: suberlak
"""

# Another program to figure out what scipy is doing with my matrix, given that 
# printing anything does not seem to work at all...

from scipy.linalg import cholesky, cho_solve, decomp_cholesky
#import decomp_cholesky


import numpy as np
from george import kernels
import george

fname = '231408.02-011355.4.dat'
data = np.loadtxt(fname)

t=mjd = data[:,0]
y=mags = data[:,1]
errs = data[:,2]

# Kernel setup 

kernel =( 0.06**2.0 )* kernels.ExpKernel(7.0) + 0.06 
gp = george.GP(kernel, mean = np.mean(y))



# First pick each x-coordinate , and assign the position number of  the row
# this is done by  parse_samples()  method  
# We can call the GP class since above we have made an instance of it

x,ins = gp.parse_samples(t,sort=False)


yerr = errs 

# another method for diagonal errors is to use  where  TINY= 1.25e-12
# yerr = float(TINY) * np.ones(len(x))  



# Calculate covariance matrix : call the kernel   class which  has an 
# instance kernel above 

K = kernel(x[:,None],x[None,:])

# taking care of the diagonal entries  ... 

K[np.diag_indices_from(K)] += yerr ** 2

# calculating the Cholesky decomposition of the Covariance matrix


#   Upper-  triangular Cholesky factor of K  
# I call here the method directly from scipy.linalg , to see what are the used parameters

# And I pulled them from George's  basic.py 

# So we factor the matrix... 
upper = decomp_cholesky.cholesky(K, lower=False, overwrite_a=False, check_finite=True)
factor=upper

#And we calculate the log-determinant
scale=0.5*np.log(2*np.pi)
const = -(np.sum(np.log(np.diag(factor))) + scale*len(x))
_x,inds  =gp.parse_samples(x)
r = gp._check_dimensions(y)[ins] - gp.mean(x)


# Calculating the ln-likelihood function as a function of the kernel parameters
#mean=np.mean(y)
#r = gp._check_dimensions(y)[ins] # - mean(x)  - I don't get what is calculated here 
# what is mean(x) ? Mean is a number, but there is this function mean which does 
# some really strange stuff... 

# ll= const - 0.5*np.dot(r.T, cho_solve(factor, r) ) 



# Now another problem lies in  the optimization algorithm - it also uses the 
# gp.compute()  routine, which has this somehow erroneous calculation of K ... 
