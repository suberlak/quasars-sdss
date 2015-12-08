# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:44:24 2015

@author: suberlak
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os 
import re 
from sf_bin_statistics_exact_one_plot_L import read_file
from sf_bin_statistics_exact_one_plot_L import p_sigma_mu

pre_obj = ['StarB', 'StarR', 'QSO']
pre_mag = ['17-19', '18.5-19']
pre_corr = ['','_corr']

files_sig = {}

# Compare across all objects and 
# all mag ranges and 
# corrections

for obj in pre_obj : 
    for mag in pre_mag:
        for corr in pre_corr :
            pre = obj+'_bins_'+mag+corr
            fname = pre +'_N_max_sig_full_sig.txt'
            
            
            # check if the file exists
            if os.path.exists(fname) : 
                num_lines = sum(1 for line in open(fname))
            
                if num_lines == 200 : 
                    
                    #print 'Reading in the results ...\n'
                    d = np.genfromtxt(fname)
                    N = d[:,0]
                    sig_max = d[:,1]
                    sig_full = d[:,2]
                    
                    a = np.where(sig_max - sig_full != 0)
                    if np.size(a) > 0:
                        print '\nFor ',  fname 
                        print 'There is some discrepancy...'
                        ind = a[0]
                        print 'Affected bins: ', N[ind]
                        for i in range(len(ind)):
                            print sig_max[ind], sig_full[ind]
                        files_sig[fname] = [obj,mag,corr,N[ind], sig_max[ind], sig_full[ind]]
            # if not, go ahead with the comparison        
            else:
                print 'This file does not exist... ', fname
                
# Plot log(L) for those particular bins...
                
def chi2plotMarginalAstro( lnL, name, sig_full, sig_max,sig_mcmc, mu_mcmc,
                           mu_s=100, sig_s=40, 
                          sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2]):
 
 
     # Calculate sigma, mu grids 
     sigma = np.linspace(sig_lim[0], sig_lim[1], sig_s)
     mu = np.linspace(mu_lim[0], mu_lim[1], mu_s)
      
     # initialize the figure 
     fig = plt.figure(figsize=(8, 8))
     fig.subplots_adjust(left=0.08, bottom=0.15, right=0.95, top=0.90, wspace=0.29, hspace=0.46)
    
     lnL[lnL < -10] = -10  # truncate for clean plotting
     
     ## lnL image
     ax = fig.add_axes([0.35, 0.35, 0.45, 0.6], xticks=[], yticks=[])
     ax.set_title('ln(L) image', fontsize=14)
     # pretty color map
     plt.imshow(lnL, origin='lower', cmap=plt.cm.RdYlGn, aspect='auto')
     # extent : left, right , bottom, top 
     # colorbar
     cax = plt.axes([0.82, 0.35, 0.02, 0.6])
     cb = plt.colorbar(cax=cax)
     cb.set_label(r'$ln(L)$', fontsize=14)
     plt.clim(np.min(lnL), np.max(lnL))
     
     # mark max of logL 
     ind = np.where(lnL == np.max(lnL))
     # mu : ind[1][0]  ,  sig: ind[0][0]
     # mu : 0-100   ,  sig : 0:40
     ax.plot(ind[1][0],ind[0][0], 'x', color='white', alpha=0.99, lw=45)
     
       # mark mcmc results ...
     ind_mm = np.interp(mu_mcmc, mu, np.arange(len(mu)))
     ind_sm = np.interp(sig_mcmc, sigma, np.arange(len(sigma)))
     ax.plot(ind_mm,ind_sm, '+', color='white', alpha=0.99, lw=45)

     # mark sig max : solid line , ought be same as the cross
     # just to confirm...
     ind_s =  np.interp(sig_max, sigma, np.arange(len(sigma)))   #np.where(abs(sigma - sig_max) < 0.0001)
     ax.axhline(y=ind_s, color='white', alpha=0.99)
     
    
     # mark sig full : dot-dash line 
     ind_f = np.interp(sig_full, sigma, np.arange(len(sigma)))  # np.where(abs(sigma - sig_full) < 0.0001)
     ax.axhline(y=ind_f, color='white', ls = '-.',alpha=0.99)
     
     # compute marginal projections
     p_mu = lnL.sum(0)  # mu
     p_sigma = lnL.sum(1)  # sigma 
     
     ax1 = fig.add_axes([0.35, 0.1, 0.45, 0.23], yticks=[])
     ax1.plot(mu, p_mu, '-k')
     ax1.set_xlabel(r'$\mu$', fontsize=12)
     ax1.set_ylabel(r'$p(\mu)$', fontsize=12)
     ax1.set_xlim(np.min(mu), np.max(mu))
     
     ax2 = fig.add_axes([0.15, 0.35, 0.18, 0.6], xticks=[])
     ax2.plot(p_sigma, sigma, '-k')
     ax2.set_xlabel(r'$p(\sigma)$', fontsize=12)
     ax2.set_ylabel(r'$\sigma$', fontsize=12)
     ax2.set_xlim(ax2.get_xlim()[::-1])  # reverse x axis
     ax2.set_ylim(np.min(sigma), np.max(sigma))
     
     
     if (name is None):
         plt.show()
     else:
         print 'saving plot to:', name
         plt.savefig(name, bbox_inches='tight')
 
     return p_mu, p_sigma



for key in files_sig.keys():
    
    obj = files_sig[key][0]
    mag = files_sig[key][1]
    corr = files_sig[key][2]
    pre = obj+'_bins_'+mag+corr
    
    for i in range(len(files_sig[key][3])):
        N = int(files_sig[key][3][i])
        bin_N = str(N).zfill(3)
        File= pre+'/'+obj+'_bin_'+bin_N+'_xi_ei'+mag+'.txt'
        print '\n'+File
        xi, ei= read_file(File)
        mu_bin, sigma_bin, logL_bin, mu, sigma = p_sigma_mu(xi,ei)
        name = obj+mag+corr+'_bin_'+bin_N+'.png'
        sig_max, sig_full= files_sig[key][4][i], files_sig[key][5][i]
        print('Sig max = %f , Sig full = %f ' % (sig_max, sig_full))
        
        # read in the mcmc result for that bin...
        mcmc_file = obj+'_bins_'+mag+corr+'_chains_res.txt'
        print 'Using mcmc result from ', mcmc_file
        mcmc_res = np.genfromtxt(mcmc_file)
        sig_mcmc = mcmc_res[:,2][mcmc_res[:,0] == N]
        mu_mcmc= mcmc_res[:,1][mcmc_res[:,0] == N]
        
        p_mu, p_sigma = chi2plotMarginalAstro(logL_bin, name, sig_full, sig_max, 
                                              sig_mcmc, mu_mcmc)