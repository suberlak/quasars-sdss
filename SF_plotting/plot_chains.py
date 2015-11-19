# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 14:04:43 2015

@author: suberlak
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['ytick.labelsize'] = 25
rcParams['xtick.labelsize'] = 25
rcParams['axes.labelsize'] = 35
rcParams['axes.linewidth'] = 3
rcParams['font.size'] = 25
rcParams.update({'figure.autolayout': False})

# a file where a mean tau per bin was saved 
taufile = np.genfromtxt('Stars_with_corr_mean_tau_all_bins.txt')
tau1 = taufile[:,0]
tau2 = taufile[:,1]


#pre = ['QSO_bins_17-19','StarB_bins_17-19','StarR_bins_17-19' ]






    
#fig = plt.figure(figsize=[12,6])
#ax = fig.add_subplot(111)



def print_panel(pre, ax,  corr, tau1, col ):
    ''' Function to print a panel of SF given
    an array of file prefixes to use...
    '''
    
    for i in range(len(pre)) : 
    
        data = np.genfromtxt(pre[i]+corr+'_chains_res.txt')
        bin_N = data[:,0]
        #mu = data[:,1]
        sig = data[:,2]
        
        # grab numbers of bins that we have
        # to use them as indices need to convert 
        # to integers and subtract 1 (because bins 
        # start from  1 )
        bin_int = [int(N)-1 for N in bin_N]
        
        ax.scatter(np.log10(tau1[bin_int]), sig, color=col[i], alpha=0.5)

col = [ 'blue', 'red', 'black']
pre_obj = ['StarB', 'StarR', 'QSO']
pre_mag = ['17-19', '18.5-19']
pre_corr = ['','_corr']




# y limits for ALL PANELS 
y_top  = 0.45
y_bott = 0

# x limits for ALL PANELS 
x_left = 0.5
x_right = 3.7

lh_w   = 1.0  # horizontal line thickness 
lh_st  = '--' # horizontal line style 
lh_al  = 0.5  # horizontal line alpha parameter 

ax_cnt = 0

for corr in pre_corr:
    
    fig,ax = plt.subplots(2,1, figsize=(12,9), sharex=True)
    fig.subplots_adjust(hspace=0)
    
    plt.gcf().subplots_adjust(bottom=0.15)
    axs = ax.ravel()

    print 'corr =', corr 
    for mag in pre_mag:
        pre = [obj + '_bins_'+mag for obj in pre_obj]        
        print pre
        print ax_cnt 
        ax = axs[ax_cnt]

        print_panel(pre,ax,corr, tau1, col )
        
        ax.set_ylabel(r'$SF $')
        ax.set_ylim(bottom=y_bott, top=y_top)
        ax.set_xlim(left=x_left, right=x_right)
        ax.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)    
        ax.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
        ax.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
        ax.set_yticks([0,0.1,0.2,0.3,0.4])
        ax.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
        text = r'$ \mathrm{mag:} \,'+ mag+'$'
        ax.text(x=0.65, y=0.3,s = text )
        ax_cnt += 1 
    ax_cnt = 0

    ax = axs[1]
    ax.set_xlabel(r'$log_{10} (\Delta _{t})$ [days]')
    #plt.tight_layout()
    name = 'SF_mag_17-19_18.5_19'+corr+'_mcmc.png'
    print 'saving the print as', name
    plt.savefig(name)