# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 19:50:39 2014

@author: astronomy

Plot the distributions  of  chi2  (not divided by N-1 because on some days there
may have been only one measurement )   , as well as mean and rms. 

This takes all days from all sources, and plots it   

"""

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec
# make appropriate arrays 

directory = 'QSO_CRTS_processed_N2/' #  QSO_try/'
names=np.loadtxt(directory+'out.list',dtype=str)

#mjd =  np.empty(0,dtype=float)
mag = np.empty(0,dtype=float)
err_rms =  np.empty(0,dtype=float)
err_weights=  np.empty(0,dtype=float)
N =  np.empty(0,dtype=float)
chisq =  np.empty(0,dtype=float)

for obj in names : 
    address=directory+obj 
    data=np.loadtxt(address) 
    #mjd = np.append(mjd,data[:,0])
    mag = np.append(mag,data[:,1])
    err_weights = np.append(err_weights,data[:,2])
    err_rms = np.append(err_rms,data[:,3])
    N = np.append(N,data[:,4]) 
    chisq = np.append(chisq,data[:,5])




nbins =100
plt.clf()
fig1 = plt.figure()
    
# Define the canvas to work on   and the  grid  
fig1 = plt.figure(figsize=[10,8])
gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)    
        
# First histogram  : Chelsea results 
H, xedges,yedges = np.histogram2d(err_weights,err_rms,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)

ax1 = fig1.add_subplot(gs[:,:90])   
pcObject1 = ax1.pcolormesh(xedges, yedges, Hmasked)

plt.ylabel('Error from daily rms',fontsize=15)
plt.xlabel('Error from weighted CRTS errors',fontsize=15)
plt.xlim((0,0.6))
plt.ylim((0,3))
plt.title('CRTS quasars : calculated  errors')
# Add the colorbar  
axC = fig1.add_subplot(gs[:,95:])
cbar = fig1.colorbar(pcObject1,ax=ax1, cax=axC, orientation='vertical')
cbar.ax.set_ylabel('Number of days with given error value ')
plt.savefig('qso_CRTS_errors_testing.png')