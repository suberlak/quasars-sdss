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

xi = datatable[:,0]  # data points
ei = datatable[:,1]  # error points 

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


    return mu_bin, sig_bin

# Split into chunks, calculate mu, sigma per chunk and average 
chunk_size = 10000 
split_xi = split_bin(xi,chunk_size)
split_ei = split_bin(ei, chunk_size)
mu_chunks, sigma_chunks = [],[]

for i in range(len(split_xi)):
    mu, sigma = p_sigma_mu(split_xi[i],split_ei[i])
    print('Indiv : mu_chunk =%f and sigma_chunk=%f'%(mu, sigma))
    mu_chunks.append(mu) , sigma_chunks.append(sigma)
mu_bin, sigma_bin = np.average(mu_chunks), np.average(sigma_chunks)    
    
print('Chunks : mu =%f and sigma=%f'%(mu_bin, sigma_bin))

# Do the calculation on the entire bin in one go 
mu_bin, sigma_bin = p_sigma_mu(xi,ei)
print('Full  : mu =%f and sigma=%f'%(mu_bin, sigma_bin))