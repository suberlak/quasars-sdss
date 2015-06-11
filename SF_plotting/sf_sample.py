# -*- coding: utf-8 -*-
"""
Created on Thu May 28 10:55:03 2015

@author: suberlak

A program to read-in a sample of QSO, Stars red, Stars blue, 
selected from my master files according to  SDSS_r < 20 ,  
CRTS_err < 0.3  ,  and log10(tau) < 1.7 . 

Note that SDSS_r is the SDSS_r magnitude, and 

Further, I cut QSO     -1.5 < delflux < 1.5 , and 
              stars    -1.0 < delflux < 1.0 

Here I am reading in those files,  and plotting  them in various ways. 

"""

import numpy as np 
import matplotlib.pyplot as plt 
from astroML.plotting import scatter_contour 
from astroML.stats import median_sigmaG
from scipy.stats import binned_statistic

# pip install astroML --user

# READ IN THE FILES  

def read_file(obj):
    '''
    By default, want to use SDSS_r filter data. However, accidentally 
    I saved SDSS_g filter data for Stars for Sample_4  . Sample_4 for QSO has 
    SDSS_r mags.   So  I made a Sample_6 , with identical criteria (hence 
    the same number of lines), but storing SDSS_r mag instead.  Thus, if need
    be,  for stars I can plot delflx / delflxerr vs SDSS_r or SDSS_g 
    
    '''


    if obj == 'StarB_g' : 
        File = 'Sample_4_Stars_blue_314248_lines.txt'
    if obj == 'StarB':
        File = 'Sample_6_Stars_blue_314248_lines.txt'
        
    if obj == 'StarR_g' :
        File = 'Sample_4_Stars_red_373306_lines.txt'
    if obj == 'StarR':
        File = 'Sample_6_Stars_red_373306_lines.txt'
        
    if obj == 'QSO' :
        File = 'Sample_6_QSO_721283_lines.txt'
        
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    data = {}
    print 'Zipping ', obj, ' data from ', File, ' ...'
    for label, column in zip(colnames, datatable.T):
        data[label] = column

    return data
    
starB_gMag = read_file('StarB_g')
starB_rMag = read_file('StarB')
starR_gMag = read_file('StarR_g') 
starR_rMag = read_file('StarR')
QSO = read_file('QSO') 


# DEFINE GAUSSIAN TO BE CONVOLVED 

def gaussian(x,mu,sigma):
    exponent = -(x-mu)**2.0 / (2.0 * (sigma ** 2.0))
    f = (1.0 / (np.sqrt(2.0*np.pi)*sigma)) * np.exp(exponent)    
    return f

def gaussgauss_logL(xi, ei, mu, sigma):
    ''' Verbatim from AstroML  '''
    """Equation 5.63: gaussian likelihood with gaussian errors"""
    ndim = len(np.broadcast(sigma, mu).shape)

    xi = xi.reshape(xi.shape + tuple(ndim * [1]))
    ei = ei.reshape(ei.shape + tuple(ndim * [1]))

    s2_e2 = sigma ** 2 + ei ** 2
    return -0.5 * np.sum(np.log(s2_e2) + (xi - mu) ** 2 / s2_e2,  -1 - ndim)

def p_sigma_mu(xi, ei, mu_s, sig_s,  get_sigma = False, sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2]):
    '''
    Instead of calculating approximate mu and sigma with approximate_mu_sigma(),
    I calculate their distribution p_sigma, p_mu, and choose the best value...
    The exact method, more time-consuming than the approximate_mu_sigma() but 
    better 
    '''
    # I assume sigma and mu range as I think they are for my distribution 
    sigma = np.linspace(sig_lim[0], sig_lim[1], sig_s)
    mu = np.linspace(mu_lim[0], mu_lim[1], mu_s)
    
    if get_sigma == False :
        logL = gaussgauss_logL(xi, ei, mu, sigma[:, np.newaxis])
        logL -= logL.max()
        L = np.exp(logL)
        
        p_sigma = L.sum(1)
        p_sigma /= (sigma[1] - sigma[0]) * p_sigma.sum()
        
        p_mu = L.sum(0)
        p_mu /= (mu[1] - mu[0]) * p_mu.sum()
        
        return p_mu, p_sigma
        
    if get_sigma == True :
        return mu , sigma
     
def approximate_mu_sigma(xi, ei, axis=None):
    """Estimates of mu0 and sigma0 via equations 5.67 - 5.68"""
    if axis is not None:
        xi = np.rollaxis(xi, axis)
        ei = np.rollaxis(ei, axis)
        axis = 0

    mu_approx, sigmaG = median_sigmaG(xi, axis=axis)
    e50 = np.median(ei, axis=axis)
    var_twiddle = (sigmaG ** 2 + ei ** 2 - e50 ** 2)
    sigma_twiddle = np.sqrt(np.maximum(0, var_twiddle))

    med = np.median(sigma_twiddle, axis=axis)
    mu = np.mean(sigma_twiddle, axis=axis)

    zeta = np.ones_like(mu)
    zeta[mu != 0] = med[mu != 0] / mu[mu != 0]

    var_approx = zeta ** 2 * sigmaG ** 2 - e50 ** 2
    sigma_approx = np.sqrt(np.maximum(0, var_approx))

    return mu_approx, sigma_approx     
     
def get_sigma_mu(xi,ei, approx=True, y_34 = 'mode', return_p = False, return_sigma=False, mu_s=100, 
             sig_s=40, sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2]):
    '''
    y_34 defines how sigma and mu are calculated based on p_sigma and p_mu 
    distributions : available   mode, exp, median (note : need to increase
    sigma and mu calculation grid to get median right !)
    
    return sigma : if set to True, it returns mu, sigma  linspaces. 
    return_p : if set to True, then it returns mu, sigma, as well as  p_mu, p_sigma 
               distributions
    mu_s : number specifying the number of grid points in mu space 
    sig_s : number specifying the number of grid points in sigma space 
    check_y34 : boolean variable (True or False), specifying whether 
    sig_lim : a list of two values: lower and upper limits on sigma linspace 
    mu_lim : a list of two values: lower and upper limits on mu linspace 
    
    '''
    #
    # RETURN ONLY SIGMA, MU LINSPACE ....
    #
    
    if return_sigma == True : 
        mu, sigma = p_sigma_mu(xi,ei, mu_s, sig_s, True, sig_lim, mu_lim)
        return mu, sigma
        
    # 
    # APPROXIMATE WAY 
    #
    if approx == True : 
        mu_i, sigma_i = approximate_mu_sigma(xi, ei, mu_s=mu_s, sig_s = sig_s)
        return mu_i, sigma_i
    
    #
    # EXACT WAY  
    #
    if approx == False : 
       
        p_mu,p_sigma = p_sigma_mu(xi,ei, mu_s, sig_s, False, sig_lim, mu_lim)
        # Change into numpy arrays 
        p_mu = np.array(p_mu)
        p_sigma = np.array(p_sigma)
        #  Retrieve sigma, mu linspace 
        mu, sigma = p_sigma_mu(xi,ei, mu_s, sig_s,  True, sig_lim, mu_lim )    
        
        # find sigma 1, 2 or 3 for this bin.... 
        
        if y_34 == 'mode' : 
            sig_max = sigma[p_sigma == max(p_sigma)][0]
            mu_max = mu[p_mu== max(p_mu)][0]
            sig_bin, mu_bin = sig_max, mu_max
            print  ' Here I am , and my sig_lim are ', sig_lim
        if y_34 == 'exp' :
            exp_val_sig = (sigma[1]-sigma[0])*np.sum(p_sigma*sigma)
            exp_val_mu = (mu[1]-mu[0])*np.sum(p_mu*mu)
            sig_bin, mu_bin = exp_val_sig, exp_val_mu
            
        if y_34 == 'median' : 
            
            # use only if enough gridding !  
            if sig_s > 100 : 
                delta_sig = sigma[1]-sigma[0]
                suma = 0
                i=0
                # integrating until the sum is 1/2 
                while suma < 0.5 : 
                    suma += delta_sig * p_sigma[i]
                    i += 1 
                sig3 = sigma[i-1]
                sig_bin = sig3
                
            if mu_s > 40 :
                
                delta_mu = mu[1]-mu[0]
                suma = 0
                i=0
                # integrating until the sum is 1/2 
                while suma < 0.5 : 
                    suma += delta_mu * p_mu[i]
                    i += 1 
                mu3 = mu[i-1]
                mu_bin = mu3
            else :
                print 'Too bad ! '
                print 'Your sigma gridding is too small to accurately calculate median....'
                print 'Calculating it with np.median instead of the proper definition '
                print ' (of integrating delta_sigma* p(sigma) until = 0.5 '
                
                sig_3 = np.median(p_sigma)
                mu_3 = np.median(p_mu)
                sig_bin, mu_bin = sig_3, mu_3
      
        # Two exit ways... 
        
        if return_p == True :
            return mu_bin, sig_bin, p_mu, p_sigma
            
        if return_p == False : 
            return mu_bin, sig_bin



# Calculate statistics for the sample ... 

def get_stats(y, y_err, mu_sig_sample):
 
    y_34    = mu_sig_sample['y_34']
    approx  = mu_sig_sample['approx'] 
    mu_s    = mu_sig_sample['mu_s']
    sig_s   = mu_sig_sample['sig_s']
    sig_lim = mu_sig_sample['sig_lim']
    mu_lim  = mu_sig_sample['mu_lim']
    
        
    # Calculate mu, sigma with method from Fig. 5.7 AstroML
    # it calculates it either approximate or full method, 
    # depending on choice of approx. 
    
    # Set higher gridding for sigma and mu,  and   
    # decrease the range of sigma and mu, just like in p_distributions_sample()
    
    mu_bin, sigma_bin = get_sigma_mu(y,y_err, approx, y_34, sig_s=sig_s, mu_s=mu_s,
                                 sig_lim = sig_lim, mu_lim = mu_lim)
    #mu_app, sigma_app = approximate_mu_sigma(xi, ei)
    SF = sigma_bin
    #mu_mean = np.mean(y)    
        
    return SF, mu_bin 

       
def get_histogram(xdata, nbins):
    # Calculate unnormalized histogram, divided by the N_sample at the end 
    N_hist_ttl = float(len(xdata))
    hist, bin_edges = np.histogram(xdata, bins=nbins, density=False) 
    # density = False instead of density = True 
    # ensures I do the normalisation my way : then I need to calculate the area 
    # under the curve, and scale up the gaussian so that they are normalised with 
    # respect to each other , but not to unit area  
    
    bin_cen = (bin_edges[:-1] + bin_edges[1:])/2 
    bin_width = (bin_edges[-1] - bin_edges[0]) / float(nbins)
    bin_width = bin_edges[1]-bin_edges[0]
    hist  = hist / (N_hist_ttl * bin_width)
    
    area = np.sum(bin_width * hist)
    
    print 'N_hist_ttl=', N_hist_ttl
    print 'Delta_bin = ', bin_width
    
    
    # Calculate normalised histogram, so that the INTEGRAL under the curve=1
    # This allows me to overplot the Gaussian too! 
    hist_n, edges = np.histogram(xdata, bins=nbins, density=True)
    # Just in case, also calculate the histogram normalized to 1  
    
    return hist, hist_n, bin_cen, area        
    
    
def combine_gaussians(xi=starB_rMag['delflx'], ei=starB_rMag['delflxerr'], hist_xlim, SF, mu ):
    xgrid = np.linspace(hist_xlim[0], hist_xlim[1], 10000)
    N = float(len(xi))
           
    
    model2 = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt(sigma ** 2.0 + ei[i]**2.0)
        Gi =  gaussian(xgrid,mu, sig_com)
        model2 += Gi  
    model2= (1.0 / N) * model2   
    
    return xgrid, model2
    
# PLOT THINGS...
# Cols are # delflx  delflxerr  obj_ID  SDSS_r_mMed


fig1 = plt.figure(figsize=(12,12))
ax1 = fig1.add_subplot(221)

# Use SDSS g magnitudes...
#ax1.scatter( starB_gMag['SDSS_g_mMed'], starB_gMag['delflx'] / starB_gMag['delflxerr']  )
scatter_contour(starB_gMag['SDSS_g_mMed'], starB_gMag['delflx'] / starB_gMag['delflxerr'],
                threshold=400, log_counts=True, ax=ax1,histogram2d_args=dict(bins=40),
                plot_args=dict(marker=',', linestyle='none', color='black'))
ax1.set_xlabel('SDSS g magnitudes')
ax1.set_ylabel(r'$\chi= \Delta (CRTS \, mag) \, / \,(CRTS \, err)$')
ax1.set_title('Blue Stars')

# Use SDSS r magnitudes... 
fig2 = plt.figure(figsize=(12,12))
ax1 = fig2.add_subplot(221)
scatter_contour(starB_rMag['SDSS_r_mMed'], starB_rMag['delflx'] / starB_rMag['delflxerr'],
                threshold=400, log_counts=True, ax=ax1,histogram2d_args=dict(bins=40),
                plot_args=dict(marker=',', linestyle='none', color='black'))


# PLOT  HISTOGRAMS and GAUSSIANS for r magnitudes.... 

# SETTINGS THAT CONTROL THE APPEARANCE OF HISTOGRAM, AND SIGMA, MU LINSPACES 
mu_sig_sample={}
colnames = ['bins_hist', 'hist_xlim', 'mu_s', 'sig_s','mu_lim', 'sig_lim','y_34', 'approx']
datatable = [200, [-1.5,1.5], 2*40, 2*100, [-0.007,0.007 ], [0.17, 0.19],'mode', True]

for label, column in zip(colnames, datatable):
    mu_sig_sample[label] = column  
    
    
# BLUE Stars   mag: 17-18 
mask = (starB_rMag['SDSS_r_mMed'] < 18) * (starB_rMag['SDSS_r_mMed'] >17)
delflx = starB_rMag['delflx'][mask]
delflxerr = starB_rMag['delflxerr'][mask]
hist, hist_n, bins, area = get_histogram(xdata=delflx, nbins=200)
SF, mu = get_stats(y=delflx, y_err=delflxerr, mu_sig_sample, mu_sig_generic) 
xgrid, model = combine_gaussians(xi=delflx, ei=delflxerr, 
                                 hist_xlim=[-1.0, 1.0], SF=SF, mu=mu )


ax2 = fig2.add_subplot(222)
ax2.plot(bins, hist,ls='steps', lw=2, label='data')
ax2.plot(xgrid, model, lw=2, label=r'$\sigma=0, \, \mu=%d, \, f_{c}=1.0$' % mu )

# BLUE Stars   mag: 18-18.5 

# BLUE Stars   mag: 18.5-19


