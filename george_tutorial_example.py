# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 15:39:50 2014

@author: suberlak
"""

#  Based on the tutorial from 
#  dan.iel.fm/george/current/user/hyper/

import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt


data = sm.datasets.get_rdataset("co2").data
t = np.array(data.time)
y = np.array(data.co2)

plt.plot(t,y)

plt.show()

from george import kernels 
reload(kernels)
k1 = 66.0**2 * kernels.ExpSquaredKernel(67.0**2)
k2 = 2.4**2 * kernels.ExpSquaredKernel(90**2) * kernels.ExpSine2Kernel(2.0 / 1.3**2, 1.0)
k3 = 0.66**2 * kernels.RationalQuadraticKernel(0.78, 1.2**2)
k4 = 0.18**2 * kernels.ExpSquaredKernel(1.6**2) + kernels.WhiteKernel(0.19)
kernel = k1 + k2 + k3 + k4

import george
gp = george.GP(kernel, mean=np.mean(y))
gp.compute(t)
print(gp.lnlikelihood(y))
print(gp.grad_lnlikelihood(y))

p, results = gp.optimize(t, y)

x = np.linspace(max(t), 2025, 2000)
mu, cov = gp.predict(y, x)
std = np.sqrt(np.diag(cov))

