# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:51:08 2015

@author: suberlak

A test code to check capability of using monte-carlo...

Based on sf_plot_NEW_panels_only.py

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os 
import re 

rcParams['ytick.labelsize'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['axes.labelsize'] = 30
rcParams['axes.linewidth'] = 2
rcParams['font.size'] = 20
# Read in xi : delflux and ei : associated error
# points for bin 1 for CRTS QSO 
def read_file(File):
    
    datatable = np.genfromtxt(File)
    xi = datatable[:,0]  # data points
    ei = datatable[:,1]  # error points 
    return xi, ei
    
def split_bin(array, length):
    ''' Split an input array into chunks of desired length
    The chunks will all have the same length, the part that
    is left (non-divisible by length) will form the last chunk'''
    
    length = max(1, length)
    return [array[i:i + length] for i in range(0, len(array), length)]
    

    
def p_sigma_mu(xi, ei, mu_s=100, sig_s=40,  sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2]):
    '''
    Make sigma, mu space and calculate the log-likelihood. p_sigma is 
    marginalization of L along one axis, and p_mu is marginalization 
    along another axis . 
    
    Parameters:
    -----------
    xi : measurement points 
    ei : error points
    mu_s :  number of grid points for mu space
    sig_s : number of grid points for sigma space
    sig_lim : limits of sigma space
    mu_lim : limits of mu space 

    Returns:
    ---------
    mu_bin : a value for mu calculated for that bin
    sig_bin : a value for sigma calculated for that bin
    '''
    # Calculate sigma, mu grids 
    sigma = np.linspace(sig_lim[0], sig_lim[1], sig_s)
    mu = np.linspace(mu_lim[0], mu_lim[1], mu_s)
      
    
    # Calculate log-likelihood  
    ndim = len(np.broadcast(sigma[:, np.newaxis], mu).shape)

    xi = xi.reshape(xi.shape + tuple(ndim * [1]))
    ei = ei.reshape(ei.shape + tuple(ndim * [1]))

    s2_e2 = sigma[:, np.newaxis] ** 2 + ei ** 2
    logL =  -0.5 * np.sum(np.log(s2_e2) + (xi - mu) ** 2 / s2_e2,  -1 - ndim)    
    
    logL -= logL.max()
    L = np.exp(logL)
    
    p_sigma = L.sum(1)
    p_sigma /= (sigma[1] - sigma[0]) * p_sigma.sum()
    
    p_mu = L.sum(0)
    p_mu /= (mu[1] - mu[0]) * p_mu.sum()
    
    # Change into numpy arrays 
    p_mu = np.array(p_mu)
    p_sigma = np.array(p_sigma)
    
    # find sigma and mu using mode 
    sig_max = sigma[p_sigma == max(p_sigma)][0]
    mu_max = mu[p_mu== max(p_mu)][0]
    sig_bin, mu_bin = sig_max, mu_max


    return mu_bin, sig_bin, logL, mu, sigma 

# Split into chunks, calculate mu, sigma per chunk and average 


def compare_sigma(pre):
    ''' A function to calculate sigma in two ways. 
    One way is full mode calculation, 
    by finding a max of 1D p(sigma) distribution, marginalizing 
    log-likelihood : sig_full
    Another way is a maximum of a 2D log-likelihood : sig_max
    
    Commented out : splitting a bin in chunks , deemed 
    unnecessary and introducing another level of uncertainty 
    
    Parameters:
    -----------
    
    Returns:
    --------
    logL_bin : the log-likelihood value (2D), 40x100, using sigma, mu space
               as defined in p_sigma_mu()
               
    N : the bin number for which a calculation is made 
    max_sig : sigma from 1D max of p(sigma)
    full_sig: sigma from a 2D max of full log-likelihood space
    '''
    
    #ch_sig = []
    full_sig = []
    max_sig = [] 
    N = []
    
    #pre = 'StarB_bins_17-19' # 'StarB_bins_17-19'
    bin_dir = pre+'/'
    bins = os.listdir(bin_dir)
    start = '_bin_'
    end = '_xi_'
    
    for i in range(len(bins)): #how_many
        
        File =  bin_dir+bins[i]
        result = re.search('%s(.*)%s' % (start,end),bins[i]).group(1)
        N_bin = float(result)
        print '\nFor bin N=', N_bin
        xi, ei= read_file(File)
        N.append(N_bin)

        
        try : 
            # Full calculation for a single bin 
            mu_bin, sigma_bin, logL_bin, mu, sigma = p_sigma_mu(xi,ei)
            full_sig.append(sigma_bin)
            print('Full  : mu =%f and sigma=%f'%(mu_bin, sigma_bin))
        
        except MemoryError:
            # Chunks calculation : split bin into chunks
            print 'Oops! Memory Error.. Splitting into chunks'
            chunk_size = 20000 
            print 'Using chunk size = ', chunk_size
            split_xi = split_bin(xi,chunk_size)
            split_ei = split_bin(ei, chunk_size)
            mu_chunks, sigma_chunks = [],[]
            
            for i in range(len(split_xi)):
                mu_ch, sigma_ch, logL, mu, sigma = p_sigma_mu(split_xi[i],split_ei[i])
                # print('Indiv : mu_chunk =%f and sigma_chunk=%f'%(mu, sigma))
                mu_chunks.append(mu_ch) , sigma_chunks.append(sigma_ch)
            mu_bin, sigma_bin = np.average(mu_chunks), np.average(sigma_chunks)    
            max_sig.append(sigma_bin)    
            print('Chunks : mu =%f and sigma=%f'%(mu_bin, sigma_bin))
            
        # Find max of logL 
        ind = np.where(logL_bin == np.max(logL_bin))
        mu_max, sigma_max = mu[ind[1]], sigma[ind[0]]
        max_sig.append(sigma_max)
        print('Max  : mu =%f and sigma=%f'%(mu_max, sigma_max))
     
    # save the results ...
    d = np.column_stack((N,max_sig, full_sig))
    fname = pre+'_N_max_sig_full_sig.txt'
    print 'Saved results of two ways for finding sigma to ...', fname
    np.savetxt(fname, d, delimiter = ' ',fmt='%s' )
    
    return logL_bin, N, max_sig, full_sig

pre_obj = ['StarB', 'StarR', 'QSO']
pre_mag = ['17-19', '18.5-19']
pre_corr = ['','_corr']

for obj in pre_obj : 
    for mag in pre_mag:
        for corr in pre_corr :
            pre = obj+'_bins_'+mag+corr
            fname = pre +'_N_max_sig_full_sig.txt'
            
            
            # check if the file was not already saved...
            if os.path.exists(fname) : 
                num_lines = sum(1 for line in open(fname))
                if num_lines == 200 : 
                    print fname
                    print 'This one was already compared, moving on...\n'
            # if not, go ahead with the comparison        
            else:
                print 'Running comparison for', fname
                logL_bin, N, max_sig, full_sig = compare_sigma(pre)
                
                
def chi2plotMarginalAstro( lnL,  mu_s=100, sig_s=40, sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2]):
 
 
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
     cb.set_label(r'$lnL(X_{centroid}, C_{mod})-lnLmax$', fontsize=14)
     plt.clim(np.min(lnL), np.max(lnL))
     
     # mark max of logL 
     ind = np.where(lnL == np.max(lnL))
     
     #ax.plot(sigtrue, 1000.0, 'o', color='red', alpha=0.75)
     # mark ML solution: (sigmaML, CmodML)
     ax.plot(sigma[ind[0]],mu[ind[1]], 'x', color='white', alpha=0.99, lw=35)
     
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
     ax2.set_ylabel(r'$\sigma)$', fontsize=12)
     ax2.set_xlim(ax2.get_xlim()[::-1])  # reverse x axis
     ax2.set_ylim(np.min(sigma), np.max(sigma))
     
     name = None
     if (name is None):
         plt.show()
     else:
         print 'saving plot to:', name
         plt.savefig(name, bbox_inches='tight')
 
     return p_mu, p_sigma
 
#p_mu, p_sigma = chi2plotMarginalAstro(logL_bin)
    
#fig = plt.figure(figsize=(12,4))
#ax = fig.add_subplot(111)
#fig.subplots_adjust(hspace=0)
#ax = axs.ravel()

#x = np.arange(1,bin_max+1)
#ax.set_ylabel(r'SF ($=\sigma$)')
#ax.set_xlabel('Bin number')
#ax.plot(x,ch_sig, label='chunks', lw=2)
#ax.plot(x,full_sig, label='full', lw=2)
#ax.legend(framealpha=0.7, loc = 'upper left')
#title = 'TEST_QSO_mags_17-19_err_0.3_no_corr_bins_1-'+str(bin_max)+'.png'
#plt.savefig(title)

