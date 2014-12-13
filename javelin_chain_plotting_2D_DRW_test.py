# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 17:01:19 2014

@author: suberlak

Plotting the Markov Chains from Javelin

INPUT : 

1) list of chains for a given band 
2) results of comparison of Chelsea's results to Javelin results 

OUTPUT : 

1) Plots of ln(tau)  vs ln(sigma) for each chain, together with overplotted 
crosshairs marking the position of 50% tau and sigma for javelin  (blue) and Chelsea (red).

NOTES:

The 50% hpd for Javelin was used to compare the results with Chelsea, so if it 
differs significantly from the plotted maximum of the chain distribution,
then 50% of hpd does not represent well the Javelin fit, and one needs to use 
a different value of sigma and tau to compare with Chelsea's fits.


"""
import numpy as np 
import matplotlib.pyplot as plt 
from math import  isinf

dir_choice = ['qso_drw_chains/','qso_drw_chains/w_prior_figs/']

dir_in = dir_choice[0]
dir_out = dir_choice[1]


'''
NOTE : must make a chain_list_  ... .ls  file before running the program!
in QSO_SDSS_chains/  run :
ls ch_u_*.txt_chain.dat > chain_list_u.ls

'''
filename = dir_in + 'chains_err_all.ls'
files=np.loadtxt(filename,dtype=str)


# load multiple chains 

for j in range(1):   #len(files)
    fchain = dir_in+files[j]
    flatchain= np.genfromtxt(fchain)
    fig1 = plt.figure()
    
    # plot log10  of the distr...  The distr is stored as  ln.... 
    # " a two-column txt file with the first column log(sigma) and the second one log(tau), both natural logs"
    # need to convert to log10  ...
    sigma = np.exp(flatchain[:,0]) 
    tau = np.exp(flatchain[:,1])
    x=np.log10(sigma)
    y=np.log10(tau)
    
    xinf = np.asarray(map(isinf,x),dtype=bool)
    yinf = np.asarray(map(isinf,y),dtype=bool)
    ttlinf = xinf + yinf
    # ttlwh = np.where(ttlinf == True)  list of good indices
    gi = -ttlinf  # good_indices 
    non_inf = len(np.where(gi == True)[0])
#    print 'Out of ', len(x),' rows, we have ', non_inf, ' of those that do not',\
#    ' have any infinities, and only those are used for plotting '
    
    plt.plot(x[gi],y[gi],'.r')
    nbins =100
    
    H, xedges,yedges = np.histogram2d(x[gi],y[gi],bins=nbins)
    
    # finding the maximum of the 2D distribution    
    a,b = np.where(H == H.max())    
    x_max = xedges[a[0]]
    y_max = yedges[b[0]]
    
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    fig2 = plt.figure()
    plt.pcolormesh(xedges, yedges, Hmasked)
    
    title = 'Chain plot for '+files[j][3:-14]
    plt.title(title)
    plt.axhline(y_max,color='r',lw=2)
    plt.axvline(x_max,color='r',lw=2)
    plt.xlabel(r'$\log_{10}{\,\sigma}$',fontsize=16)
    plt.ylabel(r'$\log_{10}{\,\tau}$',fontsize=16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    fname3=dir_out+'chain_2D_'+files[j][3:-14]+'.png'
    plt.savefig(fname3)
    print 'Saving ', fname3

print 'We plotted ', j, 'chains, out of ', len(files), ' to sample the chain distribution'

