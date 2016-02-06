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

def get_qso_catalog(catalog):
    if catalog == 's82drw':
        File = 'CRTS_SDSS_cross_matched_qso_s82drw_catalog.txt'
    if catalog == 'DB_QSO':
        File = 'CRTS_SDSS_cross_matched_qso_DB_QSO_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    qso_catalog = {}
    print 'Zipping CRTS-SDSS quasars catalog from ', File, ' ...'
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

cols1, qso_cat, qso_names = get_qso_catalog(catalog='DB_QSO') 
cols2 , star_cat= get_stars_catalog()


# Perform cuts 
def cut_qso(qso_cat=qso_cat, qso_names=qso_names, mMin=-9, mMax=19, 
            mErrMin = -9, mErrMax = 0.3):

    mask_mag = (qso_cat['r'] > mMin) * (qso_cat['r'] < mMax) 
    mask_err = (qso_cat['CRTS_avg_e'] > mErrMin) * (qso_cat['CRTS_avg_e'] < mErrMax)
    mask = mask_mag * mask_err 
    qso_id = qso_names[mask]
    print '\n These cuts reduced the number of qso  in the sample from', \
          len(qso_cat['redshift']), ' to ', len(qso_id)
    return  qso_id, mask 

def cut_stars(star_cat=star_cat, mMin=-9, mMax=19, mErrMin = -9, 
              mErrMax = 0.3, gi_Min = -1, gi_Max=1  ):

    mask_mag = (star_cat['r_mMed'] > mMin) * (star_cat['r_mMed'] < mMax) 
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
        
 
  
def model_sf(t, sf_inf=0.25, tau = 1.0):
    br = 1.0-np.exp(-t/tau)
    sf = sf_inf * np.power(br,0.5)
    return sf


def get_plotted_quantities(data, nbins, mu_sig_sample, mu_sig_generic, err_factor=1.0 ):
    '''
    Create input for sf_plot_panels()
    data : the input data for a given bin, consisting of all del_mag, tau, err per bin
    nbins : number of bins for the entire panel 
    '''    
    delflx = data[0]
    tau = data[1]
    delflxerr = data[2]*err_factor #  for error experiment 
    
    # UNPACK PARAMETERS 
         
    y_34      = mu_sig_generic['y_34']
    approx    = mu_sig_generic['approx']
    
    # Save PARAMETERS of this function, eg. for reliable testing ... 
    pars = {}    
    colnames = ['nbins', 'mu_sig_sample', 'mu_sig_generic']
    datatable = [nbins,  mu_sig_sample, mu_sig_generic]
    for label, column in zip(colnames, datatable):
        pars[label] = column  
    
       
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

    plot_data = {}
    print ' passing on the  plot data...'
    colnames = ['mean_tau', 'bin_stdev', 'err_stdev', 'bin_sigma_G', 'err_sigma_G',
                'mu_plot', 'SF', 'err_SF', 'err_mu_plot', 'err_median']
    datatable = [mean_tau, bin_stdev, err_stdev, bin_sigma_G, err_sigma_G, 
                 mu_plot, SF, err_SF, err_mu_plot, err_median]
                 
    for label, column in zip(colnames, datatable):
        plot_data[label] = column    
    
    return plot_data


def sf_plot_panels(qso_data,star_data_blue, star_data_red, sample, choice, 
                   nbins,  err_factor, approx, y_34):  
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

    mu_sig_generic={}   
    colnames = ['y_34', 'approx']
    datatable = [y_34, approx]
    for label, column in zip(colnames, datatable):
        mu_sig_generic[label] = column 
    
  
    # DEFINE HOW MU AND SIGMA IS CALCULATED FOR QUASARS
    
    mu_sig_sample={}
    colnames = [ 'mu_s', 'sig_s','mu_lim', 'sig_lim']
    datatable = [ 2*40, 2*100, [-0.007,0.007 ], [0.17, 0.19]]
    for label, column in zip(colnames, datatable):
        mu_sig_sample[label] = column  
        
    ##############################################
    # CALCULATE PLOTTED QUANTITIES  : QUASARS 
    ############################################## 
        
    qso_plot  = get_plotted_quantities(qso_data, nbins, mu_sig_sample, mu_sig_generic, err_factor )   
   
    # DEFINE HOW MU AND SIGMA IS CALCULATED FOR RED STARS
   
    mu_sig_sample={}
    colnames = [ 'mu_s', 'sig_s','mu_lim', 'sig_lim']
    datatable = [ 2*40, 2*100, [-0.01,0.01 ], [0.08, 0.10]]
    for label, column in zip(colnames, datatable):
        mu_sig_sample[label] = column 
    ##############################################
    # CALCULATE PLOTTED QUANTITIES  : RED STARS  
    ##############################################
              
    star_plot = get_plotted_quantities(star_data_blue, nbins, mu_sig_sample, mu_sig_generic, err_factor)
    #return star_plot
    
    # DEFINE HOW MU AND SIGMA IS CALCULATED FOR RED STARS
    
    mu_sig_sample={}
    colnames = [ 'mu_s', 'sig_s','mu_lim', 'sig_lim']
    datatable = [ 2*40, 2*100, [-0.005,0.005 ], [0.02, 0.05]]
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
    ' rows of master file', i, masterFiles[i]
    
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
    
    for i in range(len(masterFiles_Q)): 
        qso_data = add_tau_delflx(masterFiles_Q,inDir_Q, good_ids_Q, i, 
                                  qso_data)

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

good_ids_S_blue  = cut_stars(mMax=20, mErrMax = 0.3, gi_Min = -1, gi_Max=1)
good_ids_S_red = cut_stars(mMax=20, mErrMax = 0.3, gi_Min = 1, gi_Max=3)
good_ids_QSO, mask_qso = cut_qso(mErrMax = 0.3 , mMax = 20)

out, qso, star_b, star_r = plot_both_SF(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO, choice='NEW_1.0Eboth0.3', nbins=200, bins_hist=200,
                  err_factor=1.0, approx=False, y_34 = 'mode')
     

    