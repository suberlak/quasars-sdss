# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 13:50:01 2014

@author: suberlak

Quick LC length plot

"""

import matplotlib.pyplot as plt

import numpy as np
from pylab import *
from scipy.optimize import curve_fit

dir_choice=['qso_drw_upd/','qso_drw_long_LC/', 'qso_drw_long_LC_chains/','qso_drw_chains/', 'qso_drw_chains/no_prior/', 'qso_drw_long_LC/conv_time_test/','qso_drw_long_LC_chains/conv_time_test/']

dir_input=dir_choice[5]
dir_output=dir_choice[6]

#pre = 'drw_err2.list_pt1'
pre = 'drw_err.list'

names=np.loadtxt(dir_input+pre,dtype=str)

LC_length = np.empty(0,dtype=float)

for i in range(len(names)):
    filename  =dir_input+names[i]
    D = np.loadtxt(filename, dtype='str')
    print len(D)
    LC_length = np.append(LC_length,len(D))
    
n = LC_length
t = [845,1069,1433,1979,2764,3847,5316,7320,9744, 12863,16694,21471,27308,34002,41887] 

# I think my time counting was flawed in that the start point was always referring 
# to the same point: the length of the input dir, which overwrote the start=clock()
# in the javelin loop.  But the overwriting happened before the cont.do_mcmc, so 
# there is just an offset of 31 secs in each measurement.  

dirr = 'qso_drw_long_LC/conv_time_test/'
offset = len(dirr)

t=  np.array(t)
y=(t-714.0-offset) /  60.0   # convert to minutes

x=n[:len(y)]

#def func(x,a,b):
#    return a*np.exp(b*x)
#
## popt, pcov = curve_fit(func, x, y)
#xx = np.linspace(0,1200,1000)
#yy = func(xx,2.0,0.0085)

plt.clf()
fig = plt.figure(figsize=[8,6])
figure = fig.add_subplot(111)
figure.plot(x,np.log(y),'ko')
figure.set_ylabel('ln (Fitting time [min])',fontsize=15)
figure.set_xlabel('Lightcurve length [number of points]',fontsize=15)
plt.title('DRW lightcurve JAVELIN fitting time',fontsize=15)
plt.xlim(250,max(x)+20)
y1, y2=figure.get_ylim()
x1, x2=figure.get_xlim()
ax2=figure.twinx()
ax2.set_ylim(np.exp(y1),np.exp(y2))
#ax2.plt(range(np.exp(y2)),np.log(np.ones(100)))
ax2.set_yscale("log")
ax2.set_ylabel('Fitting time [min]',fontsize=15)

fout = dir_output + 'LC_length_vs_time_test.png'
plt.savefig(fout)


time_day = 20.0 * 24.0 * 60.0    #  20 times 1440 mins 
how_many = time_day / y

plt.clf()
fig = plt.figure(figsize=[8,6])
figure= fig.add_subplot(111)
figure.plot(x,np.log(how_many),'*')
plt.xlim(min(x)-10,max(x)+10)
plt.title('Feasibility of using JAVELIN for long lightcurves',fontsize=15)
plt.xlabel('Lightcurve length [number of points]',fontsize=15)
plt.ylabel('ln (Number of possible fits per day)',fontsize=15)
y1, y2=figure.get_ylim()
x1, x2=figure.get_xlim()
ax2=figure.twinx()
ax2.set_ylim(np.exp(y1) / 1000.0,np.exp(y2)/1000.0)
ax2.set_yscale("log")
ax2.set_ylabel('Thousands of possible fits per day',fontsize=15)


fout = dir_output + 'LC_length_vs_number_of_fits_per_day.png'
plt.savefig(fout)

data = np.column_stack((x,y))
np.savetxt(dir_output+'LC_length_vs_fitting_time_minutes.txt',data,delimiter=" ", fmt="%s")