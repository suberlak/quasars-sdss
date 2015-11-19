# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:51:08 2015

@author: suberlak

A test code to check capability of using monte-carlo...

Based on sf_plot_NEW_panels_only.py

"""
import numpy as np

# Read in xi : delflux and ei : associated error
# points for bin 1 for CRTS QSO 
File = 'QSO_bin_1_xi_ei.txt'
datatable = np.genfromtxt(File)

xi = datatable[0]  # data points
ei = datatable[1]  # error points 

def gaussgauss_logL(xi, ei, mu, sigma):
    ''' Verbatim from AstroML  
    Equation 5.63: gaussian likelihood with gaussian errors
    
    Parameters:
    -----------
    xi : measurement points 
    ei : error points
    mu : a grid of mu values where p_mu is evaluated
    sigma : a grid of sigma values where p_sigma is evaluated
    
    Returns:
    ---------
    A Gaussian log-likelihood with gaussian errors on the grid of mu, sigma
    '''
    ndim = len(np.broadcast(sigma, mu).shape)

    xi = xi.reshape(xi.shape + tuple(ndim * [1]))
    ei = ei.reshape(ei.shape + tuple(ndim * [1]))

    s2_e2 = sigma ** 2 + ei ** 2
    return -0.5 * np.sum(np.log(s2_e2) + (xi - mu) ** 2 / s2_e2,  -1 - ndim)
              
def p_sigma_mu(xi, ei, mu_s=100, sig_s=40, sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2]):
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
    # I assume sigma and mu range as I think they are for my distribution 
    sigma = np.linspace(sig_lim[0], sig_lim[1], sig_s)
    mu = np.linspace(mu_lim[0], mu_lim[1], mu_s)

    logL = gaussgauss_logL(xi, ei, mu, sigma[:, np.newaxis])
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


    return mu_bin, sig_bin
        
mu_bin, sigma_bin = p_sigma_mu(xi,ei)
print('mu =%f and sigma=%f'%(mu_bin, sigma_bin))