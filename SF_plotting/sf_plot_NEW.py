# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 13:01:59 2015

@author: suberlak

A complete rewrite of sf_plotting.py  routine. This program 

a)  Reads in the CRTS-SDSS cross-matched catalog
b)  Performs cuts based on those catalogs, saving the id_cut array
c)  Load the 'master' files (outpuf of sf_load_NEW.py), from sf_TRY dir,  and :
    - choose which rows have tau and del_mag  for object ids that are in id_cut
    - add those points to tau_hold, del_mag hold arrays
    - recalculate bins, and bin statistics: rms, std, mean, etc. 
    - replot tau vs del_mag (three lines plot)
    - replot tau vs rms_del_mag (Structure Function plot)
    
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib.mlab as mlab
from matplotlib.ticker import FuncFormatter
matplotlib.rcParams['font.size'] = 17
from scipy.stats import binned_statistic
from scipy.signal import argrelextrema
import sys
from astroML.plotting import scatter_contour
from astroML.stats import median_sigmaG
import time 
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

##############################
# READING IN CATALOG DATA    #
##############################


# Load catalogs : output of sf_CRTS_SDSS_matching_NEW.py 
# This part takes only a second to run 

def get_qso_catalog():
    File = 'CRTS_SDSS_cross_matched_qso_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    qso_catalog = {}
    print 'Zipping CRTS-SDSS quasars catalog...'
    for label, column in zip(colnames, datatable.T):
        qso_catalog[label] = column
    
    qso_names = np.genfromtxt('CRTS_SDSS_cross_matched_qso_names.txt', dtype=str)    
    for i in range(len(qso_names)):
        qso_names[i] = qso_names[i][4:-4]
    print 'Read in ', len(qso_catalog['redshift']), ', quasars from CRTS'
    return  colnames, qso_catalog, qso_names
    
def get_stars_catalog():
    File = 'CRTS_SDSS_cross_matched_stars_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    stars_catalog = {}
    print 'zipping CRTS-SDSS stars catalog...'
    for label, column in zip(colnames, datatable.T):
        stars_catalog[label] = column
        
    return  colnames, stars_catalog

cols1, qso_cat, qso_names = get_qso_catalog() 
cols2 , star_cat= get_stars_catalog()


# Perform cuts 
def cut_qso(qso_cat=qso_cat, qso_names=qso_names, mMin=-9, mMax=19, 
            mErrMin = -9, mErrMax = 0.2, zMin = 0, zMax=2 ):

    mask_mag = (qso_cat['CRTS_avg_m'] > mMin) * (qso_cat['CRTS_avg_m'] < mMax) 
    mask_err = (qso_cat['CRTS_avg_e'] > mErrMin) * (qso_cat['CRTS_avg_e'] < mErrMax)
    mask_z = (qso_cat['redshift'] > zMin ) * (qso_cat['redshift']<zMax)
    mask = mask_mag * mask_err * mask_z
    qso_id = qso_names[mask]
    print '\n These cuts reduced the number of qso  in the sample from', \
          len(qso_cat['redshift']), ' to ', len(qso_id)
    return  qso_id, mask 

def cut_stars(star_cat=star_cat, mMin=-9, mMax=19, mErrMin = -9, 
              mErrMax = 0.2, gi_Min = -1, gi_Max=1  ):

    mask_mag = (star_cat['CRTS_M'] > mMin) * (star_cat['CRTS_M'] < mMax) 
    mask_err = (star_cat['CRTS_Merr'] > mErrMin) * (star_cat['CRTS_Merr'] < mErrMax)
    SDSS_gi = star_cat['g_mMed'] - star_cat['i_mMed']
    mask_color = (SDSS_gi > gi_Min ) * (SDSS_gi < gi_Max)
    mask = mask_mag * mask_err * mask_color
    star_id_f = star_cat['crts_id'][mask]
    # convert floats to strings without comma and zeros
    star_id = np.array(["{:.0f}".format(name) for name in star_id_f])
    print '\n These cuts reduced the number of stars  in the sample from', \
          len(star_cat['CRTS_M']), ' to ', len(star_id)
    return  star_id


########################
# STATISTICS FUNCTIONS #
########################

# FROM  http://www.astroml.org/book_figures/chapter5/fig_posterior_gaussgauss.html
        #book-fig-chapter5-fig-posterior-gaussgauss


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
        
def get_sigma_mu(xi,ei, approx=True, y_34 = 'mode', return_p = False, return_sigma=False, mu_s=100, 
                 sig_s=40, check_y34 = True,sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2]):
    '''
    y_34 defines how sigma and mu are calculated based on p_sigma and p_mu 
    distributions : available   mode, exp, median (note : need to increase
    sigma and mu calculation grid to get median right !)
    
    return sigma : if set to True, it returns mu, sigma  linspaces. 
    return_p : if set to True, then it returns mu, sigma, p_mu, p_sigma 
    
    
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
        
        
#######################
# PLOTTING FUNCTIONS  #
#######################


plt.clf()    
def get_histogram(xdata, nbins,N_hist_ttl):
    # Calculate unnormalized histogram, divided by the N_sample at the end 

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
    

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'
        
def gaussian(x,mu,sigma):
    exponent = -(x-mu)**2.0 / (2.0 * (sigma ** 2.0))
    f = (1.0 / (np.sqrt(2.0*np.pi)*sigma)) * np.exp(exponent)    
    return f
  
  
def model_sf(t, sf_inf=0.25, tau = 1.0):
    br = 1.0-np.exp(-t/tau)
    sf = sf_inf * np.power(br,0.5)
    return sf


def get_plotted_quantities(data, nbins, mu_sig_sample, mu_sig_generic, err_factor=1.0 ): #approx=True, y_34 = 'mode', hist_xlim=[-2,2]
    '''
    Create input for sf_plot_panels()
    data : the input data for a given bin, consisting of all del_mag, tau, err per bin
    nbins : number of bins for the entire panel 
    '''    
    delflx = data[0]
    tau = data[1]
    delflxerr = data[2]*err_factor #  for error experiment 
    
    # UNPACK PARAMETERS 
    hist_xlim = mu_sig_sample['hist_xlim']
    bins_hist = mu_sig_sample['bins_hist']     
    y_34      = mu_sig_generic['y_34']
    approx    = mu_sig_generic['approx']
    
    # Save PARAMETERS of this function, eg. for reliable testing ... 
    pars = {}    
    colnames = ['nbins', 'mu_sig_sample', 'mu_sig_generic']
    datatable = [nbins,  mu_sig_sample, mu_sig_generic]
    for label, column in zip(colnames, datatable):
        pars[label] = column  
    
    # Select SUBSAMPLE , and treat all those points as one bin
    # Also, make sure that only points plotted on final histogram bounds are actually 
    # taken into calculation 
    
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx
    
    tau_sm = tau[mask]
    N_hist_ttl = float(len(tau_sm))
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    
    
    # For delta_Mag / err : need to limit the range fo histogram, etc. 
    #flx_err = delflx_sm / delflxerr_sm
    #flx_err = flx_err[(flx_err<5.0) * (flx_err>-5.0)]
    #N_flx_err = float(len(flx_err))
    
    # Calculate histograms... (values of counts are normalized : 
    #                          n_bin / N_in_sample)
    # I calculate histograms. Some may be obsolete : but 
    # I calculate histograms for delflx, delflx_err, tau, delflx/err ...
    print 'Calculating histograms content...'
    hist_tau , hist_tau_n, bin_tau_cen , area_tau = get_histogram(xdata=tau_sm, nbins=bins_hist,
                                           N_hist_ttl=N_hist_ttl)
    hist_delflx , hist_delflx_n, bin_delflx_cen, area_delflx= get_histogram(xdata=delflx_sm, nbins=bins_hist,
                                                N_hist_ttl=N_hist_ttl)
    hist_err, hist_err_n, bin_err_cen, area_err = get_histogram(xdata=delflxerr_sm, 
                                                      nbins=bins_hist,
                                                      N_hist_ttl=N_hist_ttl)
 #   hist_flx_err, hist_flx_n, bin_flx_err, area_flx_err = get_histogram(xdata= flx_err,
  #                                               nbins=bins_hist, N_hist_ttl = N_flx_err)                  
    
    # Calculate all the STATISTICS for the SUBSAMPLE... 
    # sigma_std,   sigma_G,  mean,  SF , etc.. 
    
    st = gauss_stats(tau_sm, delflx_sm, delflxerr_sm, mu_sig_sample, mu_sig_generic)
    # Calculate the convolved gaussian for the SUBSAMPLE : 
    # sum of gaussians over all points in my sample  ...
    
    sigma = st['SF']
    mu = st['mu_bin']
    
    sigma1= st['SF1']
    mu1 =st['mu_bin1']
    
    print 'For my histograms I am using sigma and mu calculated with ', y_34, \
    ' setting. My err0  sigma=', sigma, ' mu=', mu, ' and the 1.3*err0 ', \
    ' sigma1 = ', sigma1, ' mu1=', mu1
     
    xi = delflx_sm
    ei = delflxerr_sm
    xgrid = np.linspace(hist_xlim[0],hist_xlim[1], 1000)  # step of 0.01
    
    # Calculate the standard error model... 
    # I sum Gaussians evaluated for each point in the SUBSAMPLE, 
    # And Sigma I use for each Gaussian is a sum of sigma evaluated 
    # for the entire SUBSAMPLE using y_34 setting (mode, exp_value, median) from
    # p(sigma), p(mu)  via  gauss_stats(),  and the intrinsic point error 
    
    model = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt(sigma ** 2.0 + ei[i]**2.0)
        Gi =  gaussian(xgrid,mu, sig_com)
        model += Gi  
      
    # Normalize the sum 
    N = float(len(delflx_sm))
    model = (1.0/N) * model
    
    # Calculate augmented error model...
    err_model=1.3
    model1 = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt(sigma1 ** 2.0 + (err_model*ei[i])**2.0)
        Gi =  gaussian(xgrid,mu1, sig_com)
        model1 += Gi  
              
    # Normalize the sum 
    model1 = (1.0/N) * model1   
    
    
    #master_acc_list = data[3]
    
    # Define functions for bin statistics 
    rms_robust = lambda x : 0.7414 *(np.percentile(x,75) - np.percentile(x,25))
    rms_std = lambda x : np.std(x)
        
    nbins = nbins # ensure uniform sampling in all statistics (same bins...)
    
    # Calculate bin statistics 
    print 'Calculating bin statistics'
    
    # Pull out some tau to plot means : common to all panels 
    binned_tau = binned_statistic(tau, tau, statistic='mean', bins=nbins)
    mean_tau = binned_tau[0]
    # Take N from each bin... 'count' function works like a regular histogram
    binned_count = binned_statistic(tau, tau, statistic='count', bins=nbins)
    bin_count = binned_count[0]
    #bin_names = np.arange(1,len(binned_count[2]))
    
     # Calculate median preprocessed photometric error per bin 
    binned_err_median = binned_statistic(tau, delflxerr, statistic='median', bins=nbins) 
    err_median = binned_err_median[0]
    
    # checking for empty bins : either mean or some custom function, but not
    # count! If statistic='count', then check for 0's , and not for nan's/ 
    non_empty_bins = np.bitwise_not(np.isnan(mean_tau))
    
    # reassign number of points in a bin and  tau position 
    
    bin_count = bin_count[non_empty_bins]
    mean_tau = mean_tau[non_empty_bins]
    err_median = err_median[non_empty_bins]
      
    
    ####
    ####  Panel 1 : Standard Deviation 
    ####
    
    stdev_binned = binned_statistic(tau, delflx, statistic = rms_std, 
                                              bins=nbins)
    
    
    bin_stdev = stdev_binned[0][non_empty_bins]  
    bin_number = stdev_binned[2]  
     # since each point belongs to some bin : len(bin_number) =len(delflx)

   
    # error on standard deviation in the bin     
    err_stdev = bin_stdev / np.sqrt(2.0*(bin_count - 1.0))
    
    #####
    ##### Panel 2  : Gaussian rms  
    #####
    bin_sigma_G = binned_statistic(tau, delflx, statistic = rms_robust, 
                                      bins=nbins)[0][non_empty_bins]

    # error on Gaussian estimate of rms in the bin 
    err_sigma_G = bin_sigma_G* 1.06 / np.sqrt(bin_count)
    
    
    ##### Panel 3 : SF  , Panel 4 : mu_approx   
    
    # Perform the calculation analoguous to AstroML fig. 5.8: 
    #
    # APPROXIMATE WAY 
    #
    # I loop over bins: each mu_i, sigma_i, is an approximate value for that 
    # calculated for delflx and delflxerr in a given bin  
    
    if approx == True : 
        
        mu_approx = []
        sigma_approx=[]
        for N in np.unique(bin_number):
            xi = delflx[bin_number == N]
            ei = delflxerr[bin_number == N]
            mu_i, sigma_i = approximate_mu_sigma(xi, ei)
            mu_approx.append(mu_i)
            sigma_approx.append(sigma_i)
    
        SF = np.array(sigma_approx)
        mu_plot = np.array(mu_approx)
        
        
        err_SF = SF * 1.06 / np.sqrt(bin_count)
        err_mu_plot = bin_stdev / np.sqrt(bin_count)
    
    #
    # EXACT WAY  
    #
    # Calculating log-Likelihood for a range of values and 
    # Choosing the best one 
    # this code comes verbatim from fig. 5.8 on AstroML, but I make 
    # arrays of sigma and mu values more relevant to our problem 

    if approx == False : 
        SF = []
        mu_plot = []
        
        # Calculate sigma, mu for each bin 
        for N in np.unique(bin_number):
            xi = delflx[bin_number == N]
            ei = delflxerr[bin_number == N]
            mu_bin, sigma_bin = get_sigma_mu(xi,ei, approx, y_34)
            SF.append(sigma_bin)
            mu_plot.append(mu_bin)
        
        SF = np.array(SF)
        mu_plot = np.array(mu_plot)
            
        SF = SF.reshape(nbins,)
        mu_plot = mu_plot.reshape(nbins,)
        
        err_SF = SF * 1.06 / np.sqrt(bin_count)
        err_mu_plot = bin_stdev / np.sqrt(bin_count)

    # lines to add if checking for p_mu, p_sigma distr : 
        # 'p_mu', 'p_sigma'
        # p_mu, p_sigma 
    # histogram data not used anymore... 
    # labels 'hist_tau', 'hist_tau_n', 'bin_tau_cen',   'hist_delflx_n',
    # 'hist_err' , 'hist_err_n',  'bin_err_cen', 'area_err' 'area_tau',
    #  'hist_flx_err', 'hist_flx_n', 'bin_flx_err', 'area_flx_err', 

    plot_data = {}
    print ' passing on the  plot data...'
    colnames = ['mean_tau', 'bin_stdev', 'err_stdev', 'bin_sigma_G', 'err_sigma_G',
                'mu_plot', 'SF', 'err_SF', 'err_mu_plot', 'err_median', 
                'hist_delflx' ,'bin_delflx_cen', 'st',  'area_delflx',
                'model', 'model1', 'xgrid', 'pars']
    datatable = [mean_tau, bin_stdev, err_stdev, bin_sigma_G, err_sigma_G, 
                 mu_plot, SF, err_SF, err_mu_plot, err_median, 
                 hist_delflx , bin_delflx_cen, st,  area_delflx,
                 model, model1, xgrid, pars]
                 
    for label, column in zip(colnames, datatable):
        plot_data[label] = column    
    
    return plot_data

def gauss_stats(tau,y, y_err, mu_sig_sample, mu_sig_generic):
    '''
    Used to calculate all necessary statistics for the subsample 
    log(tau) < 1.7 , and  his_lim[0] < delflx < hist_lim[1], 
    tau, y, y_err passed   on from get_plotted_quantities 
    '''
    # UNPACK parameters...
    
    y_34    = mu_sig_generic['y_34']
    approx  = mu_sig_generic['approx'] 
    mu_s    = mu_sig_sample['mu_s']
    sig_s   = mu_sig_sample['sig_s']
    sig_lim = mu_sig_sample['sig_lim']
    mu_lim  = mu_sig_sample['mu_lim']
    
    # Define functions for bin statistics 
    rms_robust = lambda x : 0.7414 *(np.percentile(x,75) - np.percentile(x,25))
    rms_std = lambda x : np.std(x)
        
    # Calculate statistics on histograms... 
    sigma_stdev = binned_statistic(tau,y, statistic=rms_std, bins=1)[0][0]
    sigma_G = binned_statistic(tau,y, statistic=rms_robust, bins=1)[0][0]
    
    # Calculate mu, sigma with method from Fig. 5.7 AstroML
    # it calculates it either approximate or full method, 
    # depending on choice of approx. 
    
    # Set higher gridding for sigma and mu,  and   
    # decrease the range of sigma and mu, just like in p_distributions_sample()
    
    #mu_s = 2*40 ;   sig_s = 2*100
    #sig_lim = [0.17, 0.19] ;  mu_lim = [-0.007,0.007 ]
    mu_bin, sigma_bin = get_sigma_mu(y,y_err, approx, y_34, sig_s=sig_s, mu_s=mu_s,
                                 sig_lim = sig_lim, mu_lim = mu_lim)
    #mu_app, sigma_app = approximate_mu_sigma(xi, ei)
    SF = sigma_bin
    mu_mean = np.mean(y)
    
    # increase the error... 
    mu_bin1, sigma_bin1 = get_sigma_mu(y,y_err*1.3, approx, y_34,sig_s=sig_s, mu_s=mu_s,
                                 sig_lim = sig_lim, mu_lim = mu_lim)
    #mu_app, sigma_app = approximate_mu_sigma(xi, ei)
    SF1 = sigma_bin1
  
    
    # fit a gaussian
    # http://stackoverflow.com/questions/7805552/fitting-a-histogram-with-python
    # best fit of data
    #(mu_gauss, sigma_gauss) = norm.fit(y)
    #(mu_gauss_2, sigma_gauss_2) = norm.fit(flx_err)
    print 'Calculated sigma_stdev, sigma_G, mu_bin, SF_bin  for ', len(tau), 'points'
    print ' Returning statistics on the chosen subsample (log(tau)<1.7) of the  plot data...'
    stat_data = {}    
    colnames = ['sigma_stdev', 'sigma_G', 'SF', 'mu_bin','SF1', 'mu_bin1',
                 'mu_mean', 'y_34']
    datatable = [sigma_stdev, sigma_G, SF, mu_bin, SF1, mu_bin1,
                 mu_mean, y_34 ]
    for label, column in zip(colnames, datatable):
        stat_data[label] = column    
    
    return stat_data  
           
def sf_plot_panels(qso_data,star_data_blue, star_data_red, sample, choice, 
                   nbins, bins_hist, err_factor, approx, y_34):  
    '''
    NEW : instead of sf_plotting, this routine is more versatile, as it 
    calls external function  get_plotted_quantities()  to actually 
    calculate the things to be plotted : means of delmag  per bin, etc. 
    
    It plots the four panels, first getting quantities to plot for stars, 
    and then for quasars, and then plotting them altogether . It also 
    creates a separate figure that plots mu_approx on a linear scale 
    '''               
    # Pass on the parameters used in calculating mu, sigma, for all points 
    # across the entire tau range 
    qso_plot=[]   
    'Returning ...' 
    
    return qso_plot
    
    mu_sig_generic={}   
    colnames = ['y_34', 'approx']
    datatable = [y_34, approx]
    for label, column in zip(colnames, datatable):
        mu_sig_generic[label] = column 
    
    # Define parameters used in calculating mu, sigma, and histograms 
    # for the subsample. Those parameters may be different for stars or quasars
    # so they are DEFINED HERE ONCE  . 
    # This  is passed on to  get_plotted_quantities(),gauss_stats(), get_sigma_mu(), etc
    
    # DEFINE HOW MU AND SIGMA IS CALCULATED FOR QUASARS, AND DEFINE HISTOGRAM
    # BOUNDARIES 
    
    mu_sig_sample={}
    colnames = ['bins_hist', 'hist_xlim', 'mu_s', 'sig_s','mu_lim', 'sig_lim']
    datatable = [bins_hist, [-1.5,1.5], 2*40, 2*100, [-0.007,0.007 ], [0.17, 0.19]]
    for label, column in zip(colnames, datatable):
        mu_sig_sample[label] = column  
    ##############################################
    # CALCULATE PLOTTED QUANTITIES  : QUASARS 
    ############################################## 
        
    qso_plot  = get_plotted_quantities(qso_data, nbins, mu_sig_sample, mu_sig_generic, err_factor )   
    #return qso_plot
    
    # DEFINE HOW MU AND SIGMA IS CALCULATED FOR RED STARS, AND DEFINE HISTOGRAM
    # BOUNDARIES
    
    mu_sig_sample={}
    colnames = ['bins_hist', 'hist_xlim', 'mu_s', 'sig_s','mu_lim', 'sig_lim']
    datatable = [bins_hist, [-1.0,1.0], 2*40, 2*100, [-0.01,0.01 ], [0.08, 0.10]]
    for label, column in zip(colnames, datatable):
        mu_sig_sample[label] = column 
    ##############################################
    # CALCULATE PLOTTED QUANTITIES  : RED STARS  
    ##############################################
              
    star_plot = get_plotted_quantities(star_data_blue, nbins, mu_sig_sample, mu_sig_generic, err_factor)
    #return star_plot
    
    # DEFINE HOW MU AND SIGMA IS CALCULATED FOR RED STARS, AND DEFINE HISTOGRAM
    # BOUNDARIES
    
    mu_sig_sample={}
    colnames = ['bins_hist', 'hist_xlim', 'mu_s', 'sig_s','mu_lim', 'sig_lim']
    datatable = [bins_hist, [-1.0,1.0], 2*40, 2*100, [-0.005,0.005 ], [0.02, 0.05]]
    for label, column in zip(colnames, datatable):
        mu_sig_sample[label] = column 
        
    ##############################################
     # CALCULATE PLOTTED QUANTITIES  : BLUE STARS 
    ##############################################
    
    star_plot1 = get_plotted_quantities(star_data_red, nbins, mu_sig_sample, mu_sig_generic, err_factor)
    #return star_plot1
    
    # ####################################
    # Plot  the standard deviation  vs tau 
    # ####################################
    
    # nbins above, check also how many stars are in the sample :
    N_qso = len(np.unique(qso_data[3]))
    N_star_b = len(np.unique(star_data_blue[3]))
    N_star_r = len(np.unique(star_data_red[3]))
    
    # set all plot parameters
    lh_w   = 1.0  # horizontal line thickness 
    lh_st  = '--' # horizontal line style 
    lh_al  = 0.5  # horizontal line alpha parameter 
    # dot size 
    p_size = 10
    p_al   = 0.5 
    # y limits for sigma, sigma_G, SF panels 
    y_top  = 0.45
    y_bott = -0.05
    # y limits for mu approx 
    y_mu_top = 0.1
    y_mu_bott = -0.1
    # x limits for ALL PANELS 
    x_left = 0.5
    x_right = 3.7
    # colors for quasars, blue and red stars 
    col1 = 'black'
    col2 = 'blue'
    col3   = 'red'
    
    # histograms:
    lw_h = 2
    plt.clf()
    fig1 = plt.figure(figsize=(12,16))

    print 'Plotting Standard Deviation vs Delta t ... '    
    ax1 = fig1.add_subplot(411)
    
    ax1.set_title('SF '+str(N_qso)+' qso, '+str(N_star_b)+' blue, and '+ 
      str(N_star_r)+' red  stars, '+ str(nbins)+  ' bin means', fontsize=20)    
    
    # Fiducial DRW     
    # Fitting to QSOs , and y-data from panel 3 (SF)
    # Used for Panel 1 and Panel 2 
    
    xdata = qso_plot['mean_tau']
    sf = qso_plot['SF']    
    popt, pcov = curve_fit(model_sf, xdata, sf)
    y = model_sf(xdata, sf_inf=popt[0], tau = popt[1]) # tau 1 year in days 
    
    err_sig = qso_plot['err_stdev']
    sf_folded = np.sqrt((y ** 2.0)+ (err_sig ** 2.0) )
    ax1.plot(np.log10(xdata), sf_folded , lw=3, c = 'orange', ls='--')
    
    # quasars ...
    ax1.scatter(np.log10(qso_plot['mean_tau']), qso_plot['bin_stdev'], s=p_size, 
                alpha=p_al, c = col1)
    ax1.errorbar(np.log10(qso_plot['mean_tau']), qso_plot['bin_stdev'],
                 qso_plot['err_stdev'], linestyle='None', c = col1  )
                 
    # blue stars ..         
    ax1.scatter(np.log10(star_plot['mean_tau']), star_plot['bin_stdev'], s=p_size, 
                alpha=p_al, c = col2)
    ax1.errorbar(np.log10(star_plot['mean_tau']), star_plot['bin_stdev'],
                 star_plot['err_stdev'], linestyle='None', c = col2  )  

    # red stars ..         
    ax1.scatter(np.log10(star_plot1['mean_tau']), star_plot1['bin_stdev'], s=p_size, 
                alpha=p_al, c = col3)
    ax1.errorbar(np.log10(star_plot1['mean_tau']), star_plot1['bin_stdev'],
                 star_plot1['err_stdev'], linestyle='None', c = col3  )  
    
         
    ax1.set_ylabel(r'$\sigma_{stdev}$',fontsize=20)  
    ax1.tick_params( axis='x', which='both',  bottom='off', 
                    top='off', labelbottom='off') 
    ax1.set_ylim(bottom=y_bott, top=y_top)
    ax1.set_xlim(left=x_left, right=x_right)
    ax1.set_yticks([0,0.1,0.2,0.3,0.4])
    ax1.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
    ax1.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.grid(axis='x')
 
    
    # ####################################
    # Plot the Gaussian sigma  vs tau  
    # ####################################
    
    print 'Plotting sigma_Gaussian  vs Delta t ... '
    ax2 = fig1.add_subplot(412)
    
    # Fiducial DRW 
    # seems relevant... http://stackoverflow.com/questions/26058792/correct-fitting-with-scipy-curve-fit-including-errors-in-x
    # 
    
    err=qso_plot['err_median']
    sf_folded  = np.sqrt((y**2.0) + (err ** 2.0) )
    ax2.plot(np.log10(xdata), sf_folded , lw=3, c = 'orange', ls='--')
    
    # quasars 
    ax2.scatter(np.log10(qso_plot['mean_tau']), qso_plot['bin_sigma_G'],
                s=p_size, alpha=p_al,c = col1)
    ax2.errorbar(np.log10(qso_plot['mean_tau']), qso_plot['bin_sigma_G'],
                 qso_plot['err_sigma_G'], linestyle='None',c = col1)
                 
    # blue stars 
    ax2.scatter(np.log10(star_plot['mean_tau']), star_plot['bin_sigma_G'],
                s=p_size, alpha=p_al,c = col2)
    ax2.errorbar(np.log10(star_plot['mean_tau']), star_plot['bin_sigma_G'],
                 star_plot['err_sigma_G'], linestyle='None',c = col2)
                 
    # red stars 
    ax2.scatter(np.log10(star_plot1['mean_tau']), star_plot1['bin_sigma_G'],
                s=p_size, alpha=p_al,c = col3)
    ax2.errorbar(np.log10(star_plot1['mean_tau']), star_plot1['bin_sigma_G'],
                 star_plot1['err_sigma_G'], linestyle='None',c = col3)
                 
                 
    ax2.set_ylim(bottom=y_bott, top=y_top)
    ax2.set_xlim(left=x_left, right=x_right)
    ax2.set_ylabel(r'$\sigma_{G} = 0.7413 \cdot (q_{75} - q_{25}) $',fontsize=20)
    ax2.tick_params( axis='x', which='both',  bottom='off', 
                    top='off', labelbottom='off') 
    ax2.set_yticks([0,0.1,0.2,0.3,0.4])
    ax2.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
    ax2.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax2.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
    ax2.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
    ax2.grid(axis='x')
    #ax3.set_title('SF '+choice+', '+str(N_objects)+' objects, '+str(nbins)+' bin means')
 
    
    
    # ##############################################
    # Plot the SF from approximate_mu_sigma  vs tau  
    # ##############################################
    
    print ' Plotting SF vs Delta t... ' 
    ax3 = fig1.add_subplot(413)
        
    # Plot fiducial DRW, parameters calculated before panel 1
    ax3.plot(np.log10(xdata), y , lw=3, c = 'orange', ls='--')
    text = r'$ \mathrm{Model:}\ \tau=%.3f, \ SF_{\infty}=%.3f \mathrm{days}$'%(popt[0],popt[1])
    ax3.text(x=1.0, y=0.3,s = text )
    # quasars
    #return qso_plot
    ax3.scatter(np.log10(qso_plot['mean_tau']), qso_plot['SF'], s=p_size, 
                alpha=p_al,c = col1)
    ax3.errorbar(np.log10(qso_plot['mean_tau']), qso_plot['SF'], 
                 qso_plot['err_SF'], linestyle='None',c = col1)
                 
    # blue stars 
    #return star_plot
    ax3.scatter(np.log10(star_plot['mean_tau']), star_plot['SF'], s=p_size, 
                alpha=p_al,c = col2)
    ax3.errorbar(np.log10(star_plot['mean_tau']), star_plot['SF'], 
                 star_plot['err_SF'], linestyle='None',c = col2)
                 
    # red stars
    ax3.scatter(np.log10(star_plot1['mean_tau']), star_plot1['SF'], s=p_size, 
                alpha=p_al,c = col3)
    ax3.errorbar(np.log10(star_plot1['mean_tau']), star_plot1['SF'], 
                 star_plot1['err_SF'], linestyle='None',c = col3)
                 
    ax3.set_ylim(bottom=y_bott, top=y_top)
    ax3.set_xlim(left=x_left, right=x_right)
    ax3.set_ylabel(r'$SF $',fontsize=20)
    ax3.tick_params( axis='x', which='both',  bottom='off', 
                    top='off', labelbottom='off')
    ax3.grid(axis='x')
    ax3.set_yticks([0,0.1,0.2,0.3,0.4])
    ax3.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
    ax3.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)    
    ax3.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax3.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al) 
        
    # ###########################################
    # Plot mu approximate ... 
    # ###########################################
    
    ax4 = fig1.add_subplot(414)
    
    # quasars 
    ax4.scatter(np.log10(qso_plot['mean_tau']), qso_plot['mu_plot'],
                s=p_size, alpha=p_al,c = col1 )
    ax4.errorbar(np.log10(qso_plot['mean_tau']), qso_plot['mu_plot'], 
                 qso_plot['err_mu_plot'], linestyle='None',c = col1)
                
    # blue stars 
    ax4.scatter(np.log10(star_plot['mean_tau']), star_plot['mu_plot'],
                s=p_size, alpha=p_al,c = col2 )
    ax4.errorbar(np.log10(star_plot['mean_tau']), star_plot['mu_plot'], 
                 star_plot['err_mu_plot'], linestyle='None',c = col2)
                 
    # red stars 
    ax4.scatter(np.log10(star_plot1['mean_tau']), star_plot1['mu_plot'],
                s=p_size, alpha=p_al,c = col3 )
    ax4.errorbar(np.log10(star_plot1['mean_tau']), star_plot1['mu_plot'], 
                 star_plot1['err_mu_plot'], linestyle='None',c = col3)
                 
                 
    ax4.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax4.set_ylim(top=y_mu_top, bottom=y_mu_bott)
    ax4.set_xlim(left=x_left, right=x_right)
    ax4.set_yticks([-0.05,0,0.05])
    ax4.set_yticklabels(['-0.05','0.0', '0.05'])  
    ax4.set_ylabel(r'$\mu_{plot}$', fontsize=20)
    ax4.grid(axis='x')
    ax4.set_xlabel(r'$log_{10} (\Delta _{t})$ [days]',fontsize=20)
    
    title2 = 'SF_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'
    fig1.subplots_adjust(hspace=0)
    plt.savefig(title2)
    plt.show()
    
    # Make a separate  figure for the sigma_G vs delta(t) on a linear scale
    # plotting EXACTLY THE SAME  as for ax4 above, but on a linear x-scale
    # thus I am NOT TAKING np.log10  of x-quantities 
    
#    plt.clf()
#    fig2 = plt.figure(figsize=(12,4))
#    ax1 = fig2.add_subplot(111)
#    ax1.scatter(qso_plot['mean_tau'], qso_plot['bin_sigma_G'],
#                s=p_size, alpha=p_al,c = col1 )
#    ax1.errorbar(qso_plot['mean_tau'], qso_plot['bin_sigma_G'], 
#                 qso_plot['err_sigma_G'], linestyle='None',c = col1)
#                 
#    ax1.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
#    ax1.set_ylim(top=y_mu_top, bottom=y_mu_bott)
#    #ax1.set_xlim(left=x_left, right=x_right)
#    ax1.set_yticks([-0.05,0,0.05])
#    ax1.set_yticklabels(['-0.05','0.0','0.05'])  
#    ax1.set_ylabel(r'$\sigma_{G}$', fontsize=20)
#    ax1.grid(axis='x')
#    ax1.set_xlabel(r'$\Delta _{t}$ [days]',fontsize=20)
#    
#    title3 = 'Sigma_del_t_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'
#    plt.savefig(title3)
#    plt.show()
    
    #
    # ----------------------------------- 
    #
    # Make HISTOGRAMS   : plot on one subplot histograms of delflx for stars, qso
    # on another subplot, make hist of delflxerr for stars, qso 
    # ----------------------------------- 
    
    plt.clf()
    fig3 = plt.figure(figsize=(12,12))
    #fig3.suptitle('Histograms, model = '+r'$\frac{1}{N} \sum_{i=1}^{N}Gauss(\mu, \Sigma=\sqrt{\sigma ^{2} + e_{i}^{2}})$')
    
    # Subplot1 : Delta_Mag  (delflx) : quasars and model(err1)
    ax1  = fig3.add_subplot(221)  # rows, columns, subplot_nnn
    # 
    ax1.plot(qso_plot['bin_delflx_cen'],qso_plot['hist_delflx'], ls='steps', 
                      label='QSO', color='black', lw=lw_h)
    ax1.plot(qso_plot['xgrid'], qso_plot['model'], lw=lw_h, color='green',
                     label='err0')    
    ax1.set_xlabel(r'$\Delta m$')   
    ax1.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$') 
    ax1.set_xlim(left=-1.5, right=1.5)
    ax1.legend(framealpha=0.7)
    
    # Subplot2 : Delta_Mag  (delflx) : quasars and model(err3)
    
    ax2  = fig3.add_subplot(222)  # rows, columns, subplot_nnn
    # 
    ax2.plot(qso_plot['bin_delflx_cen'],qso_plot['hist_delflx'], ls='steps', 
                      label='QSO', color='black', lw=lw_h)
    ax2.plot(qso_plot['xgrid'], qso_plot['model1'], lw=lw_h, color='orange',
                     label='1.3*err0')    
    ax2.set_xlabel(r'$\Delta m$')   
    ax2.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$') 
    ax2.set_xlim(left=-1.5, right=1.5)
    ax2.legend(framealpha=0.7) 
    
    # Subplot3: Blue stars, with two error models 
    ax3  = fig3.add_subplot(223)            
    ax3.plot(star_plot['bin_delflx_cen'],star_plot['hist_delflx'], ls='steps',
             label='Blue *', color='blue',lw=lw_h)
    ax3.plot(star_plot['xgrid'], star_plot['model'], lw=lw_h, color='green',
             label='err0')
    ax3.plot(star_plot['xgrid'], star_plot['model1'], lw=lw_h, color='orange',
             label='1.3*err0')
    ax3.set_xlabel(r'$\Delta m$')   
    ax3.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$') 
    ax3.set_xlim(left=-1, right=1)
    ax3.legend(framealpha=0.7)
    
    # Subplot4: Red stars, with two error models 
    
    ax4 = fig3.add_subplot(224)
    ax4.plot(star_plot1['bin_delflx_cen'], star_plot1['hist_delflx'], ls='steps',
             label='Red *', color='red',lw=lw_h)
    ax4.plot(star_plot1['xgrid'], star_plot1['model'], lw=lw_h, color='green',
             label='err0')
    ax4.plot(star_plot1['xgrid'], star_plot1['model1'], lw=lw_h, color='orange',
             label='1.3*err0')   
             
    ax4.set_xlabel(r'$\Delta m$')   
    ax4.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$') 
    ax4.set_xlim(left=-1, right=1)
    ax4.legend(framealpha=0.7)
    
    
    
    fig3.subplots_adjust(wspace=0.3)    
    
    title = 'Hist_delflx_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'     
    plt.savefig(title)
    plt.show()
    
    
    # Old things I was plotting before :
    
    # Subplot 3 : Plot Quasar histogram 
    #ax3.plot(qso_plot['bin_delflx_cen'],qso_plot['hist_delflx'], ls='steps', 
    #                  label='Quasars', color='black', lw=lw_h)
     
    # Subplot 3 : overplot a Gaussian 
    # Gaussian 
    #bins = qso_plot['bin_delflx_cen']
    #st = qso_plot['st']
    #mu = st['mu_mean']
    #sigma = st['sigma_stdev']
    #area = qso_plot['area_delflx']
    #y = area*gaussian(bins,mu, sigma)
    #ax3.set_title('Normalized distribution'+r'$  \sigma=%.3f$' %(sigma))     
    #ax3.plot(bins, y, 'b--', lw=lw_h, label ='Gaussian')    
    
    # Subplot 4 : plot the CDF of Quasars and Gaussian CDF 
    # http://stackoverflow.com/questions/9378420/how-to-plot-cdf-in-matplotlib-in-python
#    c_gauss = np.cumsum(y)
#    c_qso = np.cumsum(qso_plot['hist_delflx'])
#    ax4.plot(bins, c_gauss, 'b--', lw=lw_h, label='Gauss')
#    ax4.plot(bins, c_qso, 'black', lw=lw_h, label='QSO')
#    ax4.set_ylabel('CDF')
#    ax4.set_xlabel(r'$\Delta m$') 
#    ax4.legend(loc=4,framealpha=0.7)
#    ax4.set_xlim(left=-2, right=2)
    
    
    return qso_plot
    


# inside the main loop : get tau, delflx from a master file, either qso or star
def add_tau_delflx(masterFiles, inDir, good_ids, i, data):
    # read in storage arrays
    delflx = data[0]  
    tau = data[1]
    err = data[2]
    master_acc_list = data[3]   
    
    # read in the i-th master file 
    master =  np.genfromtxt(inDir+masterFiles[i], dtype=str)
    master_names = master[:,3]
    unique_names = np.unique(master_names)
    
    # choose good rows 
    mask_unique = np.in1d(unique_names,good_ids)
    unique_acc = unique_names[mask_unique]
    master_mask = np.in1d(master_names, unique_acc)
    
    # accepted stars / quasars from the master files:
    master_acc = master_names[master_mask]
    print '\n We accepted', len(master_acc), ' out of ', len(master_names),\
    ' rows of master file', i
    
    # read in tau,  del_mag,  del_mag_err for quasars on the list 
    delflx = np.append(delflx, master[:,0][master_mask].astype(float))
    tau = np.append(tau, master[:,1][master_mask].astype(float))
    err = np.append(err, master[:,2][master_mask].astype(float))
    master_acc_list  = np.append(master_acc_list, master_acc)
    
    return delflx, tau, err, master_acc_list
    
def plot_both_SF(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                 good_ids_QSO, choice, nbins=200, bins_hist=50, err_factor=1.0,
                 approx=True, y_34 = 'mode'):
    inDir_S       = inDirStars
    good_ids_S_blue    = good_ids_S_blue
    good_ids_S_red    = good_ids_S_red
    inDir_Q       = inDirQSO
    good_ids_Q    = good_ids_QSO
    
    
    # Read the Stellar Master file names 
    masterFiles_S = os.listdir(inDir_S)
    delflx_S      = np.empty(0,dtype=float)
    tau_S         = np.empty(0,dtype=float)
    err_S         = np.empty(0,dtype=float)
    master_acc_list_S = np.empty(0, dtype=str)

    
    # Read the QSO Master file names 
    masterFiles_Q = os.listdir(inDir_Q)
    delflx_Q      = np.empty(0,dtype=float)
    tau_Q         = np.empty(0,dtype=float)
    err_Q         = np.empty(0,dtype=float)
    master_acc_list_Q = np.empty(0, dtype=str)
    
    # Initialize the data structures to which more and more delta_t and delta_mag
    # are addded from each consecutive master file 
    
    star_data_blue = [delflx_S, tau_S, err_S, master_acc_list_S]
    star_data_red  = [delflx_S, tau_S, err_S, master_acc_list_S]
    qso_data = [delflx_Q, tau_Q, err_Q, master_acc_list_Q]    
    
    for i in range(len(masterFiles_Q)): # len(masterFiles_Q)
        #delflx_Q, tau_Q, err_Q, master_acc_list_Q
        qso_data = add_tau_delflx(masterFiles_Q,inDir_Q, good_ids_Q, i, 
                                  qso_data)
        #delflx_S, tau_S, err_S, master_acc_list_S
        star_data_blue = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_blue, i, 
                                   star_data_blue)
        
        star_data_red = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_red, i, 
                                   star_data_red)                            
                                   
        out = sf_plot_panels(qso_data, star_data_blue, star_data_red,  i, 
                             choice, nbins, bins_hist, err_factor, approx, y_34)
    
    return out, qso_data, star_data_blue, star_data_red
    
inDirStars   = 'sf_TRY/sf_stars/'
inDirQSO = 'sf_TRY/sf_qso/'



# Experiment on how different CRTS errors affect SF of both quasars and stars 
# A step forward from plotting just quasars (as below)

# Require  Merr < 0.2 mag for quasars and stars ( big error) : same as was done 
# before 

good_ids_S_blue  = cut_stars(mMax=19, mErrMax = 0.2, gi_Min = -1, gi_Max=1)
good_ids_S_red = cut_stars(mMax=19, mErrMax = 0.2, gi_Min = 1, gi_Max=3)
good_ids_QSO, mask_qso = cut_qso(mErrMax = 0.2 , mMax = 19)

out, qso, star_b, star_r = plot_both_SF(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO, choice='1.0Eboth0.2', nbins=200, bins_hist=200,
                  err_factor=1.0, approx=True, y_34 = 'mode')
 
          
#
#
#
#
#
# ----------------------------------------------------------------
# ----------------------------------------------------------------
# ----------------------------------------------------------------            
# ----------------------------------------------------------------       
# ----------------------------------------------------------------
# CHECK CHECK CHECK  #
# ----------------------------------------------------------------          
# ----------------------------------------------------------------
# ----------------------------------------------------------------
# ----------------------------------------------------------------
# ----------------------------------------------------------------
# ----------------------------------------------------------------
# ---------------------------------------------------------------- 
#           

def get_magnitudes(good_ids_QSO=good_ids_QSO, good_ids_S_red=good_ids_S_red,
                   good_ids_S_blue = good_ids_S_blue, qso_cat = qso_cat, 
                   star_cat=star_cat, qso=qso, star_b = star_b, star_r = star_r,
                   qso_names=qso_names,pre='Sample_1_'):
    
    # get MAGNITUDES for QSO  
    delflx = qso[0]
    tau = qso[1]
    name = qso[3]
    
    hist_xlim=[-1.5,1.5]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx

    names_pre = np.unique(name[mask])    
    print '\nFor QSO, we have ', len(names_pre), 'names in the Sample'
    print 'While in the overall sample, there are', len(good_ids_QSO), 'names'
    
    mask_good = np.in1d(qso_names, good_ids_QSO) * np.in1d(qso_names, names_pre)
    acc_qso_names = qso_names[mask_good]
    print 'The overlap of those two sets are ', len(acc_qso_names), 'names'

    acc_mags_SDSS = qso_cat['M_i'][mask_good]
    acc_mags_CRTS = qso_cat['CRTS_avg_m'][mask_good]
    
    DATA = np.column_stack((acc_qso_names, acc_mags_SDSS, acc_mags_CRTS))
    outfile = pre+'QSO_'+str(len(acc_mags_SDSS))+'_objects'   
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    print 'Saving names, and mags to  to ...', outfile 
    
    print 'Mean SDSS magnitudes are ', np.mean(acc_mags_SDSS)
    print 'Median SDSS magnitudes are ', np.median(acc_mags_SDSS)
    print 'Mean CRTS magnitudes are ', np.mean(acc_mags_CRTS)
    print 'Median CRTS magnitudes are ', np.median(acc_mags_CRTS)
    
    star_names = star_cat['crts_id']
    
    # get MAGNITUDES for STARS BLUE 
    
     
    delflx = star_b[0]
    tau = star_b[1]
    name = star_b[3].astype(float)
    
    hist_xlim=[-1,1]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx

    names_pre = np.unique(name[mask]) 
    good_ids = good_ids_S_blue.astype(float)
    mask_good = np.in1d(star_names, good_ids)* np.in1d(star_names, names_pre)
    acc_star_names = star_names[mask_good]
    
    acc_mags_SDSS = star_cat['g_mMed'][mask_good]
    acc_mags_CRTS = star_cat['CRTS_M'][mask_good]
    print '\nFor Blue Stars, we have ', len(names_pre), 'names in the Sample'
    print 'While in the overall sample, there are', len(good_ids_S_blue), 'names'
    print 'The overlap of those two sets are ', len(acc_star_names), 'names'

   
    DATA = np.column_stack((acc_star_names, acc_mags_SDSS, acc_mags_CRTS))
    outfile = pre+'Stars_blue_'+str(len(acc_mags_SDSS))+'_objects'   
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    print 'Saving names, and mags  to ...', outfile 
    
    print 'Mean SDSS magnitudes are ', np.mean(acc_mags_SDSS)
    print 'Median SDSS magnitudes are ', np.median(acc_mags_SDSS)
    print 'Mean CRTS magnitudes are ', np.mean(acc_mags_CRTS)
    print 'Median CRTS magnitudes are ', np.median(acc_mags_CRTS)

    # get MAGNITUDES for STARS RED 
    
    delflx = star_r[0]
    tau = star_r[1]
    name = star_r[3].astype(float)
    
    hist_xlim=[-1,1]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx

    names_pre = np.unique(name[mask]) 
    good_ids = good_ids_S_red.astype(float)
    mask_good = np.in1d(star_names, good_ids)* np.in1d(star_names, names_pre)
    acc_star_names = star_names[mask_good]
    
    acc_mags_SDSS = star_cat['g_mMed'][mask_good]
    acc_mags_CRTS = star_cat['CRTS_M'][mask_good]
    
    print '\nFor Red Stars, we have ', len(names_pre), 'names in the Sample'
    print 'While in the overall sample, there are', len(good_ids_S_red), 'names'
    print 'The overlap of those two sets are ', len(acc_star_names), 'names'


    DATA = np.column_stack((acc_star_names, acc_mags_SDSS, acc_mags_CRTS))
    outfile = pre+'Stars_red_'+str(len(acc_mags_SDSS))+'_objects'   
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    print 'Saving names, and mags  to ...', outfile 
    
    print 'Mean SDSS magnitudes are ', np.mean(acc_mags_SDSS)
    print 'Median SDSS magnitudes are ', np.median(acc_mags_SDSS)
    print 'Mean CRTS magnitudes are ', np.mean(acc_mags_CRTS)
    print 'Median CRTS magnitudes are ', np.median(acc_mags_CRTS)
    
       
    
get_magnitudes()    
    
def pickle_sample(qso=qso, star_b = star_b, star_r = star_r, pre='Sample_1_'):
    '''
    Save the data in the sample log(tau) < 1.7 , and delflx within 
    histogram limits (for QSO, -1.5 : 1.5,   for stars -1 : 1), 
    save to a text file, and to an npz file (which may be easier to 
    read in later )
    '''    

    # QSO DATA SAVE 
    delflx = qso[0]
    tau = qso[1]
    delflxerr = qso[2]
    name = qso[3]
    
    DATA = np.column_stack((delflx,tau,delflxerr, name))
    outfile = 'QSO_'+str(len(delflx))+'_lines'   
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    #np.savez(outfile+'.npz', delflx=delflx,tau=tau,delflxerr=delflxerr, name=name)
    print 'Saving all to ...', outfile 
    
    # DEFINE THE SAMPLE (FOR QSO... )
    hist_xlim=[-1.5,1.5]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    xi = delflx_sm
    ei = delflxerr_sm
    ni = name[mask]
    
    # SAVE TO FILE .... 
    DATA = np.column_stack((xi,ei,ni))    
    outfile = pre+'QSO_'+str(len(xi))+'_lines'
    print 'Saving the sample to ...', outfile 
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    #np.savez(outfile+'.npz', xi=xi, ei=ei, ni=ni)
    
    # BLUE STARS DATA SAVE 
    
    delflx = star_b[0]
    tau = star_b[1]
    delflxerr = star_b[2]
    name = star_b[3]
    
    DATA = np.column_stack((delflx,tau,delflxerr, name))
    outfile = 'Stars_blue_'+str(len(delflx))+'_lines'  
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    print 'Saving all to ...', outfile 
     
    # DEFINE THE SAMPLE (FOR QSO... )
    hist_xlim=[-1,1]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    xi = delflx_sm
    ei = delflxerr_sm
    ni = name[mask]
    
    # SAVE TO FILE .... 
    DATA = np.column_stack((xi,ei,ni))    
    outfile = pre+'Stars_blue_'+str(len(xi))+'_lines'
    
    print 'Saving the sample to ...', outfile 
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    #np.savez(outfile+'.npz', xi=xi, ei=ei, ni=ni)
    
    # RED STARS DATA SAVE 
    
    delflx = star_r[0]
    tau = star_r[1]
    delflxerr = star_r[2]
    name = star_r[3]
    
    DATA = np.column_stack((delflx,tau,delflxerr, name))
    outfile = 'Stars_red_'+str(len(delflx))+'_lines'   
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    print 'Saving all to ...', outfile
    # DEFINE THE SAMPLE (FOR QSO... )
    hist_xlim=[-1,1]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    xi = delflx_sm
    ei = delflxerr_sm
    ni = name[mask]
    
    # SAVE TO FILE .... 
    DATA = np.column_stack((xi,ei,ni))    
    outfile = pre+'Stars_red_'+str(len(xi))+'_lines'
    print 'Saving the sample to ...', outfile 
     
    np.savetxt(outfile+'.txt', DATA, delimiter =' ', fmt="%s")
    #np.savez(outfile+'.npz', xi=xi, ei=ei, ni=ni)
    
    return  xi, ei
    
ah = pickle_sample()   
    
def p_distr_quick(qso=qso):
    
    delflx = qso[0]
    tau = qso[1]
    delflxerr = qso[2]
    name = qso[3]
    
    
    # DEFINE THE SAMPLE (FOR QSO... )
    hist_xlim=[-1.5,1.5]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    xi = delflx_sm
    ei = delflxerr_sm
    ni = name[mask]
    
    # SAVE TO FILE .... 
    DATA = np.column_stack((xi,ei,ni))    
    outfile = 'Sample_QSO_'+str(len(xi))+'_lines.txt'
    print 'Saving all to ...', outfile 
    np.savetxt(outfile, DATA, delimiter =' ', fmt="%s")
    
    
    return xi, ei
    # EXACT WAY  : increase sampling....
    
    mu_s = 40
    sig_s = 100
    
    # Calculate in three ways...  mode     
    sig_lim = [0.17, 0.22] ;  mu_lim = [-0.007,0.007 ]
    sig_ex= [] ;  mu_ex=[]
    #return xi, ei, p_mu, p_sigma
    mu_exact, sigma_exact, p_mu, p_sigma = get_sigma_mu(xi,ei, approx=False, 
                                 y_34='mode', return_p=True, sig_s=sig_s, mu_s=mu_s,
                                 sig_lim = sig_lim, mu_lim = mu_lim)
                                 
    sig_ex.append(sigma_exact); mu_ex.append(mu_exact)                                    
    print 'For log(tau) < 1.7, mode : exact sigma= ', sigma_exact  , ' mu=', mu_exact
    
    
    mu, sigma =  get_sigma_mu(xi,ei, approx=False,  return_sigma=True, 
                              sig_s=sig_s, mu_s=mu_s,
                              sig_lim=sig_lim, mu_lim=mu_lim)
 
 
 
    return xi, ei, mu, p_mu, sigma, p_sigma

#ble = p_distr_quick()    
    
def p_distributions_sample(out=out, qso=qso):
    
    delflx = qso[0]
    tau = qso[1]
    delflxerr = qso[2]
    
    # DEFINE THE SAMPLE (FOR QSO... 
    hist_xlim = out['pars']['mu_sig_sample']['hist_xlim']
    #hist_xlim=[-1.5,1.5]
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx
    #tau_sm = tau[mask]
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    xi = delflx_sm
    ei = delflxerr_sm
    
    # APPROXIMATE WAY  
    
    mu_approx, sigma_approx = approximate_mu_sigma(xi, ei) 
    print 'For log(tau) < 1.7, approx sigma= ', sigma_approx  , ' mu=', mu_approx
    
    # EXACT WAY  : increase sampling.... If I go to 300 in sig_s  I get memory error...
    mu_s = out['pars']['mu_sig_sample']['mu_s']
    sig_s = out['pars']['mu_sig_sample']['sig_s']
    #mu_s = 2*40
    #sig_s = 2*100
    
    # Calculate in three ways...  mode     
    
    sig_lim = out['pars']['mu_sig_sample']['sig_lim']
    mu_lim  =out['pars']['mu_sig_sample']['mu_lim']
    
    #sig_lim = [0.17, 0.19] ;  mu_lim = [-0.007,0.007 ]
    sig_ex= [] ;  mu_ex=[]
    
    mu_exact, sigma_exact, p_mu, p_sigma = get_sigma_mu(xi,ei, approx=False, 
                                 y_34='mode', return_p=True, sig_s=sig_s, mu_s=mu_s,
                                 sig_lim = sig_lim, mu_lim = mu_lim)
                                 
                                # get_sigma_mu(y,y_err, approx, y_34, sig_s=sig_s, mu_s=mu_s,
                                # sig_lim = sig_lim, mu_lim = mu_lim)
                                 
                                 
    sig_ex.append(sigma_exact); mu_ex.append(mu_exact)                                    
    print 'For log(tau) < 1.7, mode : exact sigma= ', sigma_exact  , ' mu=', mu_exact
    
    # Calculate in three ways...  expectation   
    
    mu_exact, sigma_exact, p_mu, p_sigma = get_sigma_mu(xi,ei, approx=False, 
                                 y_34='exp', return_p=True, sig_s=sig_s, mu_s=mu_s,
                                 sig_lim = sig_lim, mu_lim = mu_lim )
                                                
    sig_ex.append(sigma_exact); mu_ex.append(mu_exact)                                      
    print 'For log(tau) < 1.7, exp: exact sigma= ', sigma_exact  , ' mu=', mu_exact

    # Calculate in three ways...  median    
    
    mu_exact, sigma_exact, p_mu, p_sigma = get_sigma_mu(xi,ei, approx=False, 
                                 y_34='median', return_p=True, sig_s=sig_s, mu_s=mu_s,
                                 sig_lim=sig_lim, mu_lim=mu_lim)
                                                
    sig_ex.append(sigma_exact); mu_ex.append(mu_exact)                                   
    print 'For log(tau) < 1.7, median : exact sigma= ', sigma_exact  , ' mu=', mu_exact

                                   
    mu, sigma =  get_sigma_mu(xi,ei, approx=False,  return_sigma=True, 
                              sig_s=sig_s, mu_s=mu_s,
                              sig_lim=sig_lim, mu_lim=mu_lim)
 
    
    
    fig  = plt.figure(figsize=(6,12))
    ax1 = fig.add_subplot(211)
    ax1.plot(sigma, p_sigma, '-o')
    #return sigma ,p_sigma, mu, p_mu, sig_ex, mu_ex
    #for sig in sig_ex: ax1.plot(sig, p_sigma[sigma==sig], 'o')
    ax1.set_title('$Blue \, Stars,$ '+r'$\log_{10}(\tau) < 1.7 $')
    ax1.set_xlabel(r'$\sigma$')
    ax1.set_ylabel(r'$p(\sigma)$')
    ax1.text(sig_ex[0]+0.005, 300, '$\sigma_{1}=%.5f$' % sig_ex[0])
    ax1.text(sig_ex[0]+0.005, 250, '$\sigma_{2}=%.5f$' % sig_ex[1])
    ax1.text(sig_ex[0]+0.005, 200, '$\sigma_{3}=%.5f$' % sig_ex[2])
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #ax1.set_xlim(xmin=0.15 , xmax=0.25)
    
    
    ax2 = fig.add_subplot(212)
    ax2.plot(mu, p_mu, '-o')
    #for m in mu_ex: ax2.plot(m, p_mu[mu==m], 'o')
    #ax2.set_xlim(xmin=-0.05 , xmax=0.05)
    ax2.set_xlabel(r'$\mu$')
    ax2.set_ylabel(r'$p(\mu)$')
    # http://stackoverflow.com/questions/25983218/scientific-notation-colorbar-in-matplotlib
    
    # https://mkaz.com/2012/10/10/python-string-format/
    ax2.text(1.5e-3, 300, '$\mu_{1}=%.2e$' % mu_ex[0])
    ax2.text(1.5e-3, 250, '$\mu_{2}=%.2e$' % mu_ex[1])
    ax2.text(1.5e-3, 200, '$\mu_{3}=%.2e$' % mu_ex[2])
    # http://stackoverflow.com/questions/11577665/change-x-axes-scale-in-matplotlib
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax2.xaxis.major.formatter._useMathText = True
    plt.savefig('p_mu_sigma_sample_tau_1.7.png')
    plt.show()
    
    return sigma, p_sigma, mu, p_mu, sig_ex, mu_ex
    

#aaa = p_distributions_sample(out=out, qso=qso)


def combined_gaussian(out=out, qso=qso, err_factor=1.3, hist_xlim=[-1.5,1.5]):
    delflx = qso[0]
    tau = qso[1]
    delflxerr = qso[2]
    
    mask_tau =np.log10(tau)<1.7
    mask_delflx = (delflx<hist_xlim[1]) * (delflx>hist_xlim[0])
    mask = mask_tau * mask_delflx
    
    #tau_sm = tau[mask]
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    # Grab sigma and mu as calculated with the AstroML Fig. 5.7 code ... 
    sigma = out['st']['SF']
    mu = out['st']['mu_bin']
    
    sigma1= out['st']['SF1']
    mu1 =out['st']['mu_bin1']
    
    # figure out how was sigma and mu calculated from p_mu and p_sigma... 
    calc_mode = out['st']['y_34']
    
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
     
    xi = delflx_sm
    ei = delflxerr_sm
    N = float(len(delflx_sm))
    xgrid = np.linspace(hist_xlim[0], hist_xlim[1], 10000)  # step of 0.01
    
    
    # Plot individual Gaussians, to show how they are added 
    # For illustration of my model... 
    
    sum_illustrate = False
    if sum_illustrate == True :     
    
        fig = plt.figure(figsize=(6,6))
        ax1 = fig.add_subplot(111)
        title = 'Gaussians,'+\
        r'$G_{i}(\mu, \Sigma_{i}=\sqrt{\sigma^{2}+err_{i}^{2}})$'+' \n'+\
        'Combined: '+r'$\frac{1}{N}\sum_{i=1}^{N}G_{i}(xgrid)$'
        
        ax1.set_title(title, y=1.08)
        ax1.set_xlabel(r'$x_{grid} (= \Delta m)$')
        ax1.set_ylabel(r'$P(x_{i})$')
        model = np.zeros_like(xgrid)
        for i in range(len(xi)):
            sig_com = np.sqrt(sigma ** 2.0 + ei[i]**2.0)
            Gi =  gaussian(xgrid,mu, sig_com)
            model += Gi  
            if i%50 == 0 : # plot every n-th out of 36000...
                ax1.plot(xgrid, Gi)
                print 'sigma_combined=', sig_com
                
        # Normalize the sum 
        N = float(len(delflx_sm))
        model = (1.0/N) * model
        ax1.plot(xgrid, model, ls='--', lw=4, c='orange')
        
        plt.savefig('Gauss_model_on_qso.png')
        plt.show()    
        
        
        
    # Calculate with sigma=0, corr_f = 1.0
    model1 = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt(ei[i]**2.0)
        Gi =  gaussian(xgrid,mu, sig_com)
        model1 += Gi  
    model1= (1.0 / N) * model1       
    print 'For model 1, mu=', mu, 'sigma for last added gaussian=e_i', sig_com 
    
    # Calculate the standard model : sigma, mu from p(sigma), p(mu), corr_fact = 1.0
    model2 = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt(sigma ** 2.0 + ei[i]**2.0)
        Gi =  gaussian(xgrid,mu, sig_com)
        model2 += Gi  
    model2= (1.0 / N) * model2   
    print 'For model 2, mu=', mu, 'combined sigma for last added gaussian=', sig_com 
    print 'sigma from p(sigma) = ', sigma
    
     # Calculate with the sigma=0, mu=0,  and corr_factor=1.3  
    model3a = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt( (err_factor*ei[i])**2.0)
        Gi =  gaussian(xgrid,0, sig_com)
        model3a += Gi  
    model3a= (1.0 / N) * model3a
    print 'For model 3a, mu=0, combined sigma for last added G() = 1.3*ei = ', sig_com  
    
    # Calculate with the sigma=0, and corr_factor=1.3  
    model3b = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt(sigma1 ** 2.0 + (err_factor*ei[i])**2.0)
        Gi =  gaussian(xgrid,mu1, sig_com)
        model3b += Gi  
    model3b= (1.0 / N) * model3b
    
    print 'For model 3b, mu=', mu1, 'combined sigma for last added G()= ', sig_com           
    print 'sigma from p(sigma) for 1.3*err sample =', sigma1
    
    
    # Calculate area under curves...
    calc_area  = False 
    if calc_area == True :
        area = (xgrid[1]-xgrid[0])*np.sum(model)
        print 'For model, area=', area
        
        area =  (xgrid[1]-xgrid[0])*np.sum(model1)
        print 'For model1, area=', area
        
        area =  (xgrid[1]-xgrid[0])*np.sum(model0)
        print 'For model0, area=', area
        
        area =  (xgrid[1]-xgrid[0])*np.sum(model01)
        print 'For model01, area=', area
        
    # combine both models  with quasar histogram...
    hist , hist_n, bins, area = get_histogram(delflx_sm,100,N)
    
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 ,box.width*0.8, box.height])
    fig.suptitle('$QSO,$ '+r'$\log_{10}(\tau) < 1.7, $'+'\n'+str(N)+' pts, '+calc_mode)
    ax.set_xlabel(r'$\Delta m$')    
    ax.set_ylabel(r'$n(bin) / (N \cdot \Delta (bin))$')
    ax.plot(bins, hist,ls='steps', lw=2, label='data')
    ax.plot(xgrid, model1, lw=2, label=r'$\sigma=0, \, \mu=%d, \, f_{c}=1.0$' % mu )
    ax.plot(xgrid, model2, lw=2, label=r'$\sigma=%.3f, \, \mu=%.4f, f_{c}=1.0$' %(sigma,mu))
    ax.plot(xgrid, model3a, lw=2, label=r'$\sigma=0, \, \mu=0, \, f_{c}=1.3$' )
    ax.plot(xgrid, model3b, lw=2, label=r'$\sigma=%.3f, \, \mu=%.4f, f_{c}=1.3$' %(sigma1,mu1) )
    ax.set_xlim(hist_xlim[0],hist_xlim[1])
    lgd =  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))    
    plt.savefig('Gaussian_model_and_QSO_hist.png',bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()
    
    return model1, model2, model3a, N, xgrid, area, bins, hist
    
#combined = combined_gaussian()   

def median_finer_grid(out=out, sigN=100, muN=40, nbins = 2):
    '''
    I grab the quasar data, and calculate p_mu and p_sigma with much 
    finer gridding,  to see how long it takes for just one bin (and how well
    it works at all)
    Time performance: 
    http://www.pythoncentral.io/time-a-python-function/
    
    '''
    delflx = qso[0]
    tau = qso[1]
    delflxerr = qso[2]
    
    # bin the delta_mag just like in get_plotted_quantities()
    rms_std = lambda x : np.std(x)
    stdev_binned = binned_statistic(tau, delflx, statistic = rms_std, 
                                              bins=200)                                   
    bin_number = stdev_binned[2] 
    
    # bin tau, to be able to say which tau given bin number corresponds to 
    binned_tau = binned_statistic(tau, tau, statistic='mean', bins=200)
    mean_tau = binned_tau[0]
    
    p_mu1 = []
    p_sigma1 =[]
    
    a = time.time()
    for N in range(nbins): #np.unique(bin_number):
        xi = delflx[bin_number == N]
        ei = delflxerr[bin_number == N]
        p_mu_N, p_sigma_N = p_sigma_mu(xi,ei, sig_linspace=sigN, mu_linspace=muN)
        p_mu1.append(p_mu_N)
        p_sigma1.append(p_sigma_N)
        
        # arrays of p_mu and p_sigma : a distribution of p_mu and p_sigma per bin,
        # 100 values for sigma, 40 values for mu      
        p_mu = np.array(p_mu1)
        p_sigma = np.array(p_sigma1)
    
    b = time.time() - a
    print 'It took ', b, 'secs with ', N+1, 'bins, and sigN=', sigN, 'muN=', muN
    perbin = b / float(N+1)
    print 'That is ', perbin, 'secs per bin'
    # retrieve sigma and mu arrays along which p_sigma and p_mu are evaluated...
    sigma, mu = p_sigma_mu(xi,ei, sig_linspace=sigN, mu_linspace=muN, get_sigma=True)
        
    # find the median : sigma3
    
    delta_sig = sigma[1]-sigma[0]
    sig3 = []
    
    
    for N in range(nbins):
        
        suma = 0
        i=0
        while suma < 0.5 : 
            suma += delta_sig * p_sigma[N,i]
            print suma, sigma[i]
            i += 1 
        sig3.append(sigma[i-1])
    
    # Plot medians found for some bins ..
    colormap = plt.cm.spectral
    #plt.gca().set_color_cycle([colormap(j) for j in np.linspace(0, 0.9, nbins)])
    
    fig1 = plt.figure(figsize=(6,6))
    ax1 = fig1.add_subplot(111)
    ax1.set_title('Median  '+r'$\sigma_{3}$')
    ax1.set_xlabel(r'$\sigma$')
    ax1.set_ylabel(r'$p(\sigma)$')
    for N in range(nbins): 
        ax1.plot(sigma, p_sigma[N,:], color=colormap(N+5))
        ax1.plot(sig3[N],p_sigma[N,sigma==sig3[N]], 'o',color=colormap(N))

    plt.savefig('p_sigma_median_finer_grid.png')
    
    return p_mu, p_sigma, mu, sigma, sig3
 
def wrapper(func,*args, **kwargs):
    # from http://www.pythoncentral.io/time-a-python-function/
    def wrapped():
        return func(*args, **kwargs)
    return wrapped
    
#wrapped = wrapper(median_finer_grid,out=out, sigN=100, muN=40, nbins=2 )
#p_mu, p_sigma, mu, sigma,sig3 = median_finer_grid(out=out, sigN=200, muN=80, nbins=20)

  
def do_astro_ml_check(qso=qso):
    '''
    Takes the data passed on to get_plotted_quantities(), and calculates the 
    same binning.  
    The aim is to check what goes into any bins that behave abnormally, eg. 
    have weird p(mu) distribution, etc. 
    
    '''
    delflx = qso[0]
    tau = qso[1]
    delflxerr = qso[2]
    
    rms_std = lambda x : np.std(x)
     
    stdev_binned = binned_statistic(tau, delflx, statistic = rms_std, 
                                              bins=200)
                                              
    bin_number = stdev_binned[2] 
    
    binned_tau = binned_statistic(tau, tau, statistic='mean', bins=200)
    mean_tau = binned_tau[0]
    
    p_mu1 = []
    p_sigma1 =[]
    
    
    for N in np.unique(bin_number):
        xi = delflx[bin_number == N]
        ei = delflxerr[bin_number == N]
        p_mu_N, p_sigma_N = p_sigma_mu(xi,ei)
        p_mu1.append(p_mu_N)
        p_sigma1.append(p_sigma_N)
        
        # arrays of p_mu and p_sigma : a distribution of p_mu and p_sigma per bin,
        # 100 values for sigma, 40 values for mu      
        p_mu = np.array(p_mu1)
        p_sigma = np.array(p_sigma1)
        
        # retrieve sigma and mu arrays along which p_sigma and p_mu are evaluated...
        sigma, mu = p_sigma_mu(xi,ei, get_sigma=True)
        
    # find the mode of each bin distribution 
    sig_1 = []
    mu_1  = []
    for N in range(0,200):
       
        sig_max = sigma[p_sigma[N,:] == max(p_sigma[N,:])][0]
        sig_1.append(sig_max)
        
        mu_max = mu[p_mu[N,:] == max(p_mu[N,:])][0]
        mu_1.append(mu_max)   
    
    
    # check where values are abnormal : here I'm fishing out  those curves 
    # which have the biggest value at biggest mu: most likely continuing to 
    # rise... 
    
    bad_mu = np.where(mu_1 == max(mu))[0]
    fig = plt.figure(figsize=(12,6))
    fig.suptitle(r'$p_{\mu} > max(\mu)$')
    ax1 = fig.add_subplot(121)
    ax1.set_xlabel(r'$\Delta(m)$')
    ax1.set_ylabel('Counts')
    ax2 = fig.add_subplot(122)
    ax2.set_xlabel(r'$err(\Delta(m))$')
    for N in bad_mu : 
        xi = delflx[bin_number == N]
        ei = delflxerr[bin_number == N]
        hist, bin_edges = np.histogram(xi)
        
        print 'For N=', N, 'log10(tau) is ', np.log10(mean_tau[N])
        
        ax1.plot(bin_edges[:-1], hist, ls='steps', label='N='+str(N), lw=2 )
        
        hist, bin_edges = np.histogram(ei)
        ax2.plot(bin_edges[:-1], hist, ls='steps', label='N='+str(N), lw=2 )
    ax1.set_xlim(-1.5,1.5)  
    ax1.legend(framealpha=0.7)
    ax2.legend(framealpha=0.7)
    
    
    plt.savefig('bad_mu_hist_qso.png')    
    plt.show()
    return p_mu, p_sigma, mu, sigma, mu_1, sig_1

#p_mu, p_sigma, mu, sigma, mu_1, sig_1 = do_astro_ml_check()



def do_mu_sigma_check(out=out):
    '''
    Here I assume that p_mu and p_sigma distributions as calculated in 
    get_plotted_quantities() are passed as the output, and that grids are 
    40 points for mu and 200 for sigma.  Given that , I plot here the 
    distributions to visualise. I also reproduce fig. 5.8 using parts of the 
    same code that I use to make p(mu) and p(sigma), to show that it all 
    behaves as it should. Thus I can also test here sigma finding code. 
    '''
    mu = np.linspace(-0.2, 0.2, 40)
    p_mu = out['p_mu']
    
    sigma = np.linspace(0,0.5, 100)
    p_sigma = out['p_sigma']
    
    fig = plt.figure(figsize=(6,12))
    ax1 = fig.add_subplot(211)
    for N in range(0,200): ax1.plot(mu,p_mu[N,:])
    ax1.set_ylabel(r'$p_{\mu}$')
    ax1.set_xlabel(r'$\mu$')
    
    ax2 = fig.add_subplot(212)    
    for N in range(0,200):ax2.plot(sigma,p_sigma[N,:])
    ax2.set_ylabel(r'$p_{\sigma}$')
    ax2.set_xlabel(r'$\sigma$')
    
    plt.savefig('p_mu_p_sigma_plot_stars.png')
    plt.show()
    
    reproduce_fig_5_8 = False 
    if reproduce_fig_5_8 == True :
        
        # Code verbatim from http://www.astroml.org/book_figures/chapter5/fig_posterior_gaussgauss.html#book-fig-chapter5-fig-posterior-gaussgauss
        #--------------------------------------------------
        # Generate data
        np.random.seed(5)
        mu_true = 1.
        sigma_true = 1.
        N = 10
        ei = 3 * np.random.random(N)
        xi = np.random.normal(mu_true, np.sqrt(sigma_true ** 2 + ei ** 2))
        
        sigma = np.linspace(0.001, 5, 70)
        mu = np.linspace(-3, 5, 70)
        
        logL = gaussgauss_logL(xi, ei, mu, sigma[:, np.newaxis])
        logL -= logL.max()
        L = np.exp(logL)
        
        p_sigma = L.sum(1)
        p_sigma /= (sigma[1] - sigma[0]) * p_sigma.sum()
        
        p_mu = L.sum(0)
        p_mu /= (mu[1] - mu[0]) * p_mu.sum()  
        fig = plt.figure(figsize=(6,12))
        ax1 = fig.add_subplot(211)
        # plot the marginalized distribution
        ax1.plot(mu, p_mu, '-k', label='marginalized')
        ax1.set_ylabel(r'$p_{\mu}$')
        ax1.set_xlabel(r'$\mu$')
        
        ax2 = fig.add_subplot(212)    
        ax2.plot(sigma, p_sigma, '-k', label='marginalized')
        ax2.set_ylabel(r'$p_{\sigma}$')
        ax2.set_xlabel(r'$\sigma$')
        plt.savefig('p_mu_p_sigma_plot_astroML.png')
        plt.show()
        
    # find  sigma_1
    
    sig_1_arg = []
    #sig_1_max = []
    for N in range(0,200):
        max_ind = argrelextrema(p_sigma[N,:],np.greater)[0] 
        sig_1_arg.append(sigma[max_ind])
        
        #sig_max = sigma[p_sigma[N,:] == max(p_sigma[N,:])]
        #sig_1_max.append(sig_max)
        
    # find sigma_2     
    sig_2 = []
    for N in range(0,200):
        exp_val = (sigma[1]-sigma[0])*np.sum(p_sigma[N,:]*sigma)
        sig_2.append(exp_val)
        
    
#do_mu_sigma_check()
  
def do_err_check(good_ids_Q =  good_ids_QSO, inDir_Q =inDirQSO ):
    
    # plot the histogram of the entry catalog filtered by good ids     
    
    
    # plot the histogram of the filtered out master files... 
    
    masterFiles_Q = os.listdir(inDir_Q)
    
    delflx_Q      = np.empty(0,dtype=float)
    tau_Q         = np.empty(0,dtype=float)
    err_Q         = np.empty(0,dtype=float)
    master_acc_list_Q = np.empty(0, dtype=str)
    qso_data =  [delflx_Q, tau_Q, err_Q, master_acc_list_Q]
    
    i = 0
    
    qso_data = add_tau_delflx(masterFiles_Q,inDir_Q, good_ids_Q, i, 
                                  qso_data)
    return qso_data
    
#qso_input = do_err_check()


def do_check(out,qso):
    delflx= qso[0]
    tau=qso[1]
    delflxerr=qso[1]
    
    mask = np.log10(tau) < 1.7
    
    tau = tau[mask]
    delflx = delflx[mask]
    delflxerr = delflxerr[mask]
   
    
    st = gauss_stats(tau, delflx, delflxerr)

    fig = plt.figure(figsize=(6,6))    
    ax1 = fig.add_subplot(111)
 
    hist = out['hist_delflx']
    bins = out['bin_delflx_cen']
    area = out['area_delflx']
    
    ax1.plot(bins, hist, ls='steps', label='Quasars', color='black', lw=2)
    y = area * gaussian(bins, st['mu_gauss'], st['sigma_gauss'])
    ax1.plot(bins, y, 'r--', lw=2)
    
    ax1.grid(True)
    ax1.set_xlabel(r'$\Delta M$')
    ax1.set_ylabel(r'$n_{bin} / N_{sample}$')
    ax1.set_title(r'$  \sigma=%.3f$' %( st['sigma_gauss']))
    ax1.grid(True)
    
    # The CDF 
    fig2 = plt.figure(figsize=(6,6))    
    ax2 = fig2.add_subplot(111)    
    
    #dx = (max(bins)-min(bins)) / float(len(bins))
    cy = np.cumsum(y)
    cy_hist = np.cumsum(hist)
    ax2.plot(bins, cy, 'b', label='gauss', lw=2)
    ax2.plot(bins, cy_hist, 'r', label='qso', lw=2)
    ax2.legend()
    # another check : from 
    #    http://stackoverflow.com/questions/23251759/how-to-determine-what-is-the-probability-distribution-function-from-a-numpy-arra
    
    fig3 = plt.figure(figsize=(6,6))    
    ax2 = fig3.add_subplot(111)    
    
    n, bins, patches = ax2.hist(delflx, 100, normed=1)
    mu = np.mean(delflx)
    sigma = np.std(delflx)
    ax2.plot(bins, mlab.normpdf(bins, mu, sigma))
    ax2.set_xlabel(r'$\Delta M$')
    plt.show()
    
    plt.clf()
    
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(111)
    x = out['mean_tau']
    sf = out['SF']
    plt.scatter(np.log10(x), sf)         
    plt.plot(np.log10(x), model_sf(x, sf_inf=0.25, tau=365))
    plt.show()
    return st, hist, bins


#st, hist, bins = do_check(out,qso)                  
                  
# Require  Merr < 0.1 mag for quasars and stars (medium error)
#
#good_ids_S_blue  = cut_stars(mMax=20, mErrMax = 0.1, gi_Min = -1, gi_Max=1)
#good_ids_S_red = cut_stars(mMax=20, mErrMax = 0.1, gi_Min = 1, gi_Max=3)
#good_ids_QSO = cut_qso(mErrMax = 0.1 , mMax = 20)
#
#out = plot_both_SF(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
#                  good_ids_QSO, choice='both0.1')
                  
# Require  Merr < 0.05 mag for quasars and stars (small error)

#good_ids_S_blue  = cut_stars(mMax=20, mErrMax = 0.05, gi_Min = -1, gi_Max=1)
#good_ids_S_red = cut_stars(mMax=20, mErrMax = 0.05, gi_Min = 1, gi_Max=3)
#good_ids_QSO = cut_qso(mErrMax = 0.05 , mMax = 20) 
#
#out = plot_both_SF(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
#                  good_ids_QSO, choice='both0.05')

# Plot statistics on quasars : simple one-function routine 
#plot_qso_stats(qso_cat)

# Experiment on how CRTS error affects the SF, but here it's only for quasars  

# Plot smallest error Quasars  : Merr < 0.05 mag 
#good_ids = cut_qso(mErrMax = 0.05 , mMax = 20)  
#out_qso = plot_single_SF(inDir='sf_TRY/sf_qso/',choice='qso0.05', good_ids = good_ids )


# Plot small error Quasars  : Merr < 0.1 mag 
#good_ids = cut_qso(mErrMax = 0.1 , mMax = 20)  
#out_qso = plot_single_SF(inDir='sf_TRY/sf_qso/',choice='qso0.1', good_ids = good_ids )

# Plot big error Quasars : Merr < 0.2 Mag 
#good_ids = cut_qso(mErrMax = 0.2, mMax = 20 )  
#out_qso = plot_single_SF(inDir='sf_TRY/sf_qso/',choice='qso0.2', good_ids = good_ids)

#plot_SF(inDir='sf_TRY/sf_stars/',good_ids = cut_stars(), choice='star' )



# Unused  functions..... 

def plotting(tau,delflx) :   
    '''
    Simple plot routine just to show tau vs delflx 
    '''
    plt.scatter(tau, delflx,s=4)
    plt.xlabel('Time difference [days]')
    plt.ylabel(r'$\Delta$ m')
    
    
    
def plot_three_lines(delflx, mean_tau, bin_means, tau, bin_rms_std, 
                     bin_rms_robust, choice, sample, nbins):
    '''
    Quickly plot the three lines : and Zeljko meant to plot those, and +/-
    '''
                         
    print 'Plotting three lines...'
    plt.clf()
    fig1 = plt.figure(figsize=(12,6))
    ax1 = fig1.add_subplot(111)
    scatter_contour(tau, delflx, threshold=len(tau)/4500, log_counts=True, 
                    ax=ax1, histogram2d_args=dict(bins=40),
                plot_args=dict(marker=',', linestyle='none', color='black'),
                contour_args=dict(cmap=plt.cm.bone))   
    #plotting(tau,delflx)
    ax1.plot(mean_tau, bin_means, color='Gold', label='Mean', lw = 2)
    ax1.plot(mean_tau, bin_means + bin_rms_std, color='r',lw = 2, 
             label='Mean+/-RMS_std')
    ax1.plot(mean_tau, bin_means - bin_rms_std, color='r',lw = 2)
    ax1.plot(mean_tau, bin_means - bin_rms_robust, color='Magenta',lw = 2 ,
             label='Mean+/-RMS_robust')
    ax1.plot(mean_tau, bin_means + bin_rms_robust, color='Magenta',lw = 2)
    ax1.set_title(r'Flux difference vs tau, nbins ='+str(nbins)+
                 ', npoints='+ str(len(tau))) 
    ax1.set_xlabel('Time difference [days]')
    ax1.set_ylabel(r'$\Delta$ m')
    ax1.set_ylim(top=2.0, bottom=-2.0)
    ax1.legend()
    title1 = 'SF_'+choice+'tau-vs-del_mag_'+sample+'.png'
    plt.savefig(title1)
    plt.show()


def sf_plotting(tau,delflx, delflxerr, master_names , sample ,choice, 
                nbins=100):
    
    '''
    
    OLD OLD OLD OLD OLD OLD OLD OLD     
    
    Note  : it was used to just plot quasars, now DEPRECEATED and superseded 
    with  sf_plot_panels , which plots quasars and stars 
    
    A more complex plotting routine, plotting the binned time difference vs 
    mean del_mag in the bin, as well as rms  on del_mag within the bin, 
    calculated in two ways:   
    rms_std     is np.std(del_mag_in_the_bin)   
    rms_robust  is np.percentile(mag_bin,75) - np.precentile(mag_bin,25)
    
    '''
    
    # Define functions for bin statistics 
    rms_robust = lambda x : 0.7414 *(np.percentile(x,75) - np.percentile(x,25))
    rms_std = lambda x : np.std(x)
    nbins = nbins # ensure uniform sampling in all statistics (same bins...)
    
    # Calculate bin statistics 
    print 'Calculating bin statistics'
    
    # Pull out some tau to plot means : common to all panels 
    binned_tau = binned_statistic(tau, tau, statistic='mean', bins=nbins)
    mean_tau = binned_tau[0]
    # Take N from each bin... 'count' function works like a regular histogram
    binned_count = binned_statistic(tau, tau, statistic='count', bins=nbins)
    bin_count = binned_count[0]
    #bin_names = np.arange(1,len(binned_count[2]))
    
    # checking for empty bins : either mean or some custom function, but not
    # count! If statistic='count', then check for 0's , and not for nan's/ 
    non_empty_bins = np.bitwise_not(np.isnan(mean_tau))
    
    # reassign number of points in a bin and  tau position 
    
    bin_count = bin_count[non_empty_bins]
    mean_tau = mean_tau[non_empty_bins]
    
    ####
    ####  Panel 1 : Standard Deviation 
    ####
    
    stdev_binned = binned_statistic(tau, delflx, statistic = rms_std, 
                                              bins=nbins)
    
    bin_stdev = stdev_binned[0][non_empty_bins]  
    bin_number = stdev_binned[2]  
     # since each point belongs to some bin : len(bin_number) =len(delflx)

   
    # error on standard deviation in the bin     
    err_stdev = bin_stdev / np.sqrt(2.0*(bin_count - 1.0))
    
    #####
    ##### Panel 2  : Gaussian rms  
    #####
    bin_sigma_G = binned_statistic(tau, delflx, statistic = rms_robust, 
                                      bins=nbins)[0][non_empty_bins]

    # error on Gaussian estimate of rms in the bin 
    err_sigma_G = bin_sigma_G* 1.06 / np.sqrt(bin_count)
    
    
    ##### Panel 3 : SF  , Panel 4 : mu_approx   
    
    # Perform the calculation analoguous to AstroML fig. 5.8: 
    
    mu_approx = []
    sigma_approx=[]
    for N in np.unique(bin_number):
        xi = delflx[bin_number == N]
        ei = delflxerr[bin_number == N]
        mu_i, sigma_i = approximate_mu_sigma(xi, ei)
        mu_approx.append(mu_i)
        sigma_approx.append(sigma_i)

    SF = np.array(sigma_approx)
    mu_approx = np.array(mu_approx)
    
    
    err_SF = SF * 1.06 / np.sqrt(bin_count)
    err_mu_approx = bin_stdev / np.sqrt(bin_count)
    
    # ########### 
    # Plot three lines ... 
    # ###########
    #bin_means = binned_statistic(tau, delflx, statistic = 'mean', bins=nbins)[0]
    #plot_three_lines(tau=tau,delflx=delflx, ax1=ax1, mean_tau=mean_tau, 
    #                     bin_means=bin_means, bin_rms_std=bin_stdev, 
    #                     bin_rms_robust=bin_sigma_G, choice=choice,
    #                     sample=sample, nbins=nbins)
    
    
    # ####################################
    # Plot  the standard deviation  vs tau 
    # ####################################
    
    # nbins above, check also how many stars are in the sample :
    N_objects = len(np.unique(master_names))
    lh_w   = 1.0  # horizontal line thickness 
    lh_st  = '--' # horizontal line style 
    lh_al  = 0.5  # horizontal line alpha parameter 
    p_size = 10
    p_al   = 0.5 
    y_top  = 0.6
    y_bott = -0.05
    x_left = 0.5
    x_right = 3.7
    
    plt.clf()
    fig1 = plt.figure(figsize=(12,16))

    print 'Plotting Standard Deviation vs Delta t ... '    
    ax1 = fig1.add_subplot(411)
    ax1.scatter(np.log10(mean_tau), bin_stdev, s=p_size, alpha=p_al)
    ax1.errorbar(np.log10(mean_tau), bin_stdev,err_stdev, linestyle='None' )
    ax1.set_ylabel(r'$\sigma_{stdev}$',fontsize=20)
    ax1.set_title('SF '+choice+', '+str(N_objects)+' objects, '+str(nbins)+
    ' bin means', fontsize=20)
    ax1.tick_params( axis='x', which='both',  bottom='off', 
                    top='off', labelbottom='off') 
    ax1.set_ylim(bottom=y_bott, top=y_top)
    ax1.set_xlim(left=x_left, right=x_right)
    ax1.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax1.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4', '0.5'])
    ax1.axhline(y=0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.grid(axis='x')
 
    
    # ####################################
    # Plot the Gaussian sigma  vs tau  
    # ####################################
    
    print 'Plotting sigma_Gaussian  vs Delta t ... '
    ax2 = fig1.add_subplot(412)
    ax2.scatter(np.log10(mean_tau), bin_sigma_G, s=p_size, alpha=p_al)
    ax2.errorbar(np.log10(mean_tau), bin_sigma_G, err_sigma_G, linestyle='None')
    ax2.set_ylim(bottom=y_bott, top=y_top)
    ax2.set_xlim(left=x_left, right=x_right)
    ax2.set_ylabel(r'$\sigma_{G} = 0.7413 \cdot (q_{75} - q_{25}) $',fontsize=20)
    ax2.tick_params( axis='x', which='both',  bottom='off', 
                    top='off', labelbottom='off') 
    ax2.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax2.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4', '0.5'])
    ax2.axhline(y=0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax2.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
    ax2.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
    ax2.grid(axis='x')
    #ax3.set_title('SF '+choice+', '+str(N_objects)+' objects, '+str(nbins)+' bin means')
 
    
    
    # ##############################################
    # Plot the SF from approximate_mu_sigma  vs tau  
    # ##############################################
    
    print ' Plotting SF vs Delta t... ' 
    ax3 = fig1.add_subplot(413)
    ax3.scatter(np.log10(mean_tau), SF, s=p_size, alpha=p_al)
    ax3.errorbar(np.log10(mean_tau), SF, err_SF, linestyle='None')
    ax3.set_ylim(bottom=y_bott, top=y_top)
    ax3.set_xlim(left=x_left, right=x_right)
    ax3.set_ylabel(r'$SF $',fontsize=20)
    ax3.tick_params( axis='x', which='both',  bottom='off', 
                    top='off', labelbottom='off')
    ax3.grid(axis='x')
    ax3.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax3.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4', '0.5']) 
    ax3.axhline(y=0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax3.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax3.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al) 
        
    # ###########################################
    # Plot mu approximate ... 
    # ###########################################
    
    ax4 = fig1.add_subplot(414)
    ax4.scatter(np.log10(mean_tau), mu_approx,s=p_size, alpha=p_al )
    ax4.errorbar(np.log10(mean_tau), mu_approx, err_mu_approx, linestyle='None')
    ax4.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax4.set_ylim(top=0.6, bottom=-0.2)
    ax4.set_xlim(left=x_left, right=x_right)
    ax4.set_yticks([-0.1,0,0.1,0.2,0.3,0.4,0.5])
    ax4.set_yticklabels(['-0.1','0.0','0.1', '0.2', '0.3', '0.4', '0.5'])  
    ax4.set_ylabel(r'$\mu_{approx}$', fontsize=20)
    ax4.grid(axis='x')
    ax4.set_xlabel(r'$log_{10} (\Delta _{t})$ [days]',fontsize=20)
    
    title2 = 'SF_'+choice+'_'+str(nbins)+'_bins_'+sample+'.png'
    fig1.subplots_adjust(hspace=0)
    plt.savefig(title2)
    plt.show()
  
    return SF
 

def plot_single_SF(inDir,  choice, good_ids):
    # Load the 'master' files... 
    inDir = inDir
    good_ids = good_ids
    investigate_bins = 'off'  # set 'on' to play with number of bins : 
               # MAKE SURE  that then the range in the loop below is not large!  
    
    masterFiles = os.listdir(inDir)
    delflx = np.empty(0,dtype=float)
    tau    = np.empty(0,dtype=float)
    err    = np.empty(0,dtype=float)
    master_acc_list = np.empty(0, dtype=str)
        
    for i in range(3): # len(masterFiles)
        master = np.genfromtxt(inDir+masterFiles[i], dtype=str)
        master_names = master[:,3]
        unique_names = np.unique(master_names)

        # choose which rows are from quasars that we want... 
        mask_unique = np.in1d(unique_names,good_ids)
        unique_acc = unique_names[mask_unique]
        master_mask = np.in1d(master_names, unique_acc)
        
        # accepted quasars from the master files:
        master_acc = master_names[master_mask]
        print '\n We accepted', len(master_acc), ' out of ', len(master_names),\
        ' rows of master file', i
        
        # read in tau,  del_mag,  del_mag_err for quasars on the list 
        delflx = np.append(delflx, master[:,0][master_mask].astype(float))
        tau = np.append(tau, master[:,1][master_mask].astype(float))
        err = np.append(err, master[:,2][master_mask].astype(float))
        master_acc_list  = np.append(master_acc_list, master_acc)
        
        # Calculate appropriate bin size from Knuths rule...
        bin_size = 2.7 * np.std(delflx) / np.power(float(len(delflx)), 0.25)
        bin_knuth = int( (max(tau) - min(tau)) / bin_size)
        print 'According to Knuths rule, bin number should be  ', bin_knuth
        
        # plotting(tau, delflx)
                          
        if investigate_bins == 'on' : 
            for n_bins in np.arange(1,11)*100 : 
                print ' Plotting SF for  ', len(tau), 'points, with ', n_bins, ' bins' 
                outcome = sf_plotting(tau,delflx, err, master_names = master_acc_list, nbins=n_bins, 
                          sample=str(i), choice=choice )
        else:
            n_bins = 200
            print ' Plotting SF for  ', len(tau), 'points, with ', n_bins, ' bins' 
            outcome = sf_plotting(tau,delflx, err, master_names = master_acc_list, nbins=n_bins, 
                          sample=str(i), choice=choice )  
                          
    return outcome 
       # sf_plotting(tau,delflx, err, master_names = master_acc_list, nbins=300, 
        #            choice='qso-knuth', sample=str(i) )

def plot_qso_stats(qso_cat):
    plt.clf()
    nbins = 20  
    fig2 = plt.figure(figsize=(12,12))
    ax1 = fig2.add_subplot(311)
    ax1.hist(qso_cat['CRTS_avg_m'], bins=nbins)
    ax1.set_title('CRTS Quasars Magnitude')
    
    ax2 = fig2.add_subplot(312)
    ax2.hist(qso_cat['CRTS_avg_e'], bins=nbins)
    ax2.set_title('CRTS Quasars Error')
    
    ax3 = fig2.add_subplot(313)
    ax3.hist(qso_cat['redshift'], bins=nbins)
    ax3.set_title('CRTS Quasars Redshift')
    
    
    fig2.subplots_adjust(hspace=0.5)
    
    plt.savefig('QSO_CRTS_histograms.png')
    
def old_way_get_mu_sigma(delflx, delflxerr, approx, bin_number, y_34, nbins, bin_count, bin_stdev):
    if approx == 'else' : 
        p_mu =[]
        p_sigma=[]
        for N in np.unique(bin_number):
            xi = delflx[bin_number == N]
            ei = delflxerr[bin_number == N]
            p_mu_N, p_sigma_N = p_sigma_mu(xi,ei)
            p_mu.append(p_mu_N)
            p_sigma.append(p_sigma_N)
        
        # arrays of p_mu and p_sigma : a distribution of p_mu and p_sigma per bin,
        # 100 values for sigma, 40 values for mu      
        p_mu = np.array(p_mu)
        p_sigma = np.array(p_sigma)
        
        
        # from p_sigma distribution extract mode (sig1), expectation value (sig2),
        # median (sig3)
        
        # retrieve sigma and mu arrays along which p_sigma and p_mu are evaluated...
        sigma, mu = p_sigma_mu(xi,ei, get_sigma=True)    
        
    
        # find mu_3
        if y_34 == 'mode' : 
            # find  sigma_1 finding sigma for which probability is maximized 
            sig_1 = []
            mu_1  = []
            for N in range(0,nbins):
                #max_ind = argrelextrema(p_sigma[N,:],np.greater)[0] 
                #sig_1.append(sigma[max_ind])
                
                sig_max = sigma[p_sigma[N,:] == max(p_sigma[N,:])][0]
                sig_1.append(sig_max)
                
                mu_max = mu[p_mu[N,:] == max(p_mu[N,:])][0]
                mu_1.append(mu_max)
            
            SF = np.array(sig_1)
            mu_plot = np.array(mu_1)
            
        if y_34 == 'exp' :
            # find sigma_2     
            sig_2 = []
            mu_2 = []
            for N in range(0,nbins):
                exp_val_sig = (sigma[1]-sigma[0])*np.sum(p_sigma[N,:]*sigma)
                sig_2.append(exp_val_sig)   
            
                exp_val_mu = (mu[1]-mu[0])*np.sum(p_mu[N,:]*mu)
                mu_2.append(exp_val_mu)
            
            SF = np.array(sig_2)
            mu_plot = np.array(mu_2)
        
        if y_34 == 'median': 
            # find sigma_3
            sig_3 = []
            mu_3 = []
            for N in range(0,nbins):
                sig_3.append(np.median(p_sigma[N,:]))
                mu_3.append(np.median(p_mu[N,:]))
                
            SF = np.array(sig_3)
            mu_plot = np.array(mu_3)
            
        SF = SF.reshape(nbins,)
        mu_plot = mu_plot.reshape(nbins,)
        
        err_SF = SF * 1.06 / np.sqrt(bin_count)
        err_mu_plot = bin_stdev / np.sqrt(bin_count)
        
def old_histograms():
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax1.plot(qso_plot['bin_delflx_cen'],qso_plot['hist_delflx'], ls='steps', 
                      label='QSO', color='black', lw=lw_h)
    # blue stars
    ax1.plot(star_plot['bin_delflx_cen'],star_plot['hist_delflx'], ls='steps',
             label='Blue *', color='blue',lw=lw_h)
    
    # red stars 
    ax1.plot(star_plot1['bin_delflx_cen'], star_plot1['hist_delflx'], ls='steps',
             label='Red *', color='red',lw=lw_h)
    
    # gaussian  with sigma from sigma_std
    # Gaussian 
    bins = qso_plot['bin_delflx_cen']
    #hist = qso_plot['hist_delflx']
    st = qso_plot['st']
    mu = st['mu_mean']
    sigma = st['sigma_stdev']
    area = qso_plot['area_delflx']
    y = area*gaussian(bins,mu, sigma)
    ax1.plot(bins, y, 'g--', lw=lw_h, label='Gauss QSO')
 
    # add a 'best fit' line
    #y = mlab.normpdf( qso_plot['bin_delflx_cen'], mu_gauss, sigma_gauss)
    #l = ax1.plot(bins, y, 'r--', linewidth=2)
    
    #ax1.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu_gauss, sigma_gauss))
    #ax1.grid(True)
         
    
    
    ax1.set_xlabel(r'$\Delta m$')   
    ax1.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$') 
    #ax1.set_yticklabels([' '])
    ax1.set_xlim(left=-2, right=2)
    ax1.legend(framealpha=0.7)

    
    # Subplot 2 : plot Delta_Mag / delflxerr 
    
    ax2 = fig3.add_subplot(222)
    # Quasars
    ax2.plot(qso_plot['bin_flx_err'], qso_plot['hist_flx_err'], ls='steps', 
                      label='Quasars', color='black', lw=lw_h)
                                      
   
    
    # red stars 
  
    
    ax2.set_xlabel(r'$\Delta m /  err(\Delta m)$')   
    ax2.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$') 
    ax2.set_xlim(left=-5, right=5)
    
    
    # Overplot a Gaussian fitted to QSO distribution of delflx / delflxerr
    bins = qso_plot['bin_flx_err']
    st = qso_plot['st']
    mu = st['mu_gauss_2']
    sigma = st['sigma_gauss_2']
    area = qso_plot['area_flx_err']
    y = area*gaussian(bins,mu, sigma)
    ax2.plot(bins, y, 'g--', lw=lw_h)
    ax2.legend(framealpha=0.7)
    
    # Subplot3 : Plot red stars error distr.
    ax3 = fig3.add_subplot(223)
    ax3.set_xlabel(r'$\Delta m$') 
    ax3.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$')
    ax2.plot(star_plot1['bin_flx_err'], star_plot1['hist_flx_err'], ls='steps',
             label='Red stars', color='red',lw=lw_h)
    ax3.set_xlim(left=-2, right=2)
    ax3.legend(framealpha=0.7)
    
    # Subplot4 : plot blue stars error distr. 
    ax4 = fig3.add_subplot(224)
    ax4.set_xlabel(r'$\Delta m$') 
    ax4.set_ylabel(r'$n_{bin} / (N * \Delta_{bin})$')
    ax4.plot(star_plot['bin_flx_err'], star_plot['hist_flx_err'], ls='steps',
             label='Blue stars', color='blue',lw=lw_h)
    
    fig3.subplots_adjust(wspace=0.3)    
    
    title = 'Hist_delflx_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'     
    plt.savefig(title)
    plt.show()
    