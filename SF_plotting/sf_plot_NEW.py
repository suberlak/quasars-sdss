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
import sys
from astroML.plotting import scatter_contour
from astroML.stats import median_sigmaG

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
    """Equation 5.63: gaussian likelihood with gaussian errors"""
    ndim = len(np.broadcast(sigma, mu).shape)

    xi = xi.reshape(xi.shape + tuple(ndim * [1]))
    ei = ei.reshape(ei.shape + tuple(ndim * [1]))

    s2_e2 = sigma ** 2 + ei ** 2
    return -0.5 * np.sum(np.log(s2_e2) + (xi - mu) ** 2 / s2_e2,
                         -1 - ndim)
              
def p_sigma_mu(xi, ei):
    '''
    Instead of calculating approximate mu and sigma with approximate_mu_sigma(),
    I calculate their distribution p_sigma, p_mu, and choose the best value...
    The exact method, more time-consuming than the approximate_mu_sigma() but 
    better 
    '''
    # I assume sigma and mu range as I think they are for my distribution 
    sigma = np.linspace(0.00, 0.5, 100)
    mu = np.linspace(-0.2, 0.2, 40)
    
    logL = gaussgauss_logL(xi, ei, mu, sigma[:, np.newaxis])
    logL -= logL.max()
    L = np.exp(logL)
    
    p_sigma = L.sum(1)
    p_sigma /= (sigma[1] - sigma[0]) * p_sigma.sum()
    
    p_mu = L.sum(0)
    p_mu /= (mu[1] - mu[0]) * p_mu.sum()
    
    return p_mu, p_sigma
#######################
# PLOTTING FUNCTIONS  #
#######################


plt.clf()    
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

def get_histogram(xdata, nbins,N_hist_ttl):
    # Calculate unnormalized histogram, divided by the N_sample at the end 

    hist, bin_edges = np.histogram(xdata, bins=nbins, density=False) 
    # density = False instead of density = True 
    # ensures I do the normalisation my way : then I need to calculate the area 
    # under the curve, and scale up the gaussian so that they are normalised with 
    # respect to each other , but not to unit area  
    
    bin_cen = (bin_edges[:-1] + bin_edges[1:])/2 
    hist  = hist / N_hist_ttl
    print 'N_hist_ttl=', N_hist_ttl
    
    bin_width = (bin_edges[-1] - bin_edges[0]) / float(nbins)
    area = np.sum(bin_width * hist)
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


def get_plotted_quantities(data,  nbins,bins_hist, err_factor=1.0):
    '''
    Create input for sf_plot_panels()
    '''    
    delflx = data[0]
    tau = data[1]
    delflxerr = data[2]*err_factor #  for error experiment 
    
    # Pick delflx and delflxerr values for histograms
    mask =np.log10(tau)<1.7
    tau_sm = tau[mask]
    N_hist_ttl = float(len(tau_sm))
    delflx_sm = delflx[mask]
    delflxerr_sm = delflxerr[mask]
    # For delta_Mag / err : need to limit the range fo histogram, etc. 
    flx_err = delflx_sm / delflxerr_sm
    flx_err = flx_err[(flx_err<5.0) * (flx_err>-5.0)]
    N_flx_err = float(len(flx_err))
    # Calculate histograms... (values of counts are normalized : 
    #                          n_bin / N_in_sample)
    print 'Calculating histograms content...'
    hist_tau , hist_tau_n, bin_tau_cen , area_tau = get_histogram(xdata=tau_sm, nbins=bins_hist,
                                           N_hist_ttl=N_hist_ttl)
    hist_delflx , hist_delflx_n, bin_delflx_cen, area_delflx= get_histogram(xdata=delflx_sm, nbins=bins_hist,
                                                N_hist_ttl=N_hist_ttl)
    hist_err, hist_err_n, bin_err_cen, area_err = get_histogram(xdata=delflxerr_sm, 
                                                      nbins=bins_hist,
                                                      N_hist_ttl=N_hist_ttl)
    hist_flx_err, hist_flx_n, bin_flx_err, area_flx_err = get_histogram(xdata= flx_err,
                                                 nbins=bins_hist, N_hist_ttl = N_flx_err)                  
    
    # Calculate all the statistics for the sample... 
    # sigma_std,   sigma_G,  mean,  SF , etc.. 
    
    st = gauss_stats(tau_sm, delflx_sm, delflxerr_sm, flx_err)
    
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
    # Approximate way 
    # I loop over bins: each mu_i, sigma_i, is an approximate value for that 
    # calculated for delflx and delflxerr in a given bin  
    
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
    
    # Exact way : calculating log-Likelihood for a range of values and 
    # Choosing the best one 
    # this code comes verbatim from fig. 5.8 on AstroML, but I make 
    # arrays of sigma and mu values more relevant to our problem 

  
    p_mu =[]
    p_sigma=[]
    for N in np.unique(bin_number):
        xi = delflx[bin_number == N]
        ei = delflxerr[bin_number == N]
        mu_i_p, sigma_i_p = p_sigma_mu(xi,ei)
        p_mu.append(mu_i_p)
        p_sigma.append(sigma_i_p)
        
    p_mu = np.array(p_mu)
    p_sigma = np.array(p_sigma)
    
    plot_data = {}
    print ' passing on the  plot data...'
    colnames = ['mean_tau', 'bin_stdev', 'err_stdev', 'bin_sigma_G', 'err_sigma_G',
                'mu_approx', 'SF', 'err_SF', 'err_mu_approx', 'err_median', 
                'hist_tau', 'hist_tau_n', 'bin_tau_cen', 'hist_delflx' ,
                'hist_delflx_n',  'bin_delflx_cen', 'hist_err' , 'hist_err_n', 
                'bin_err_cen', 'st', 'area_tau', 'area_delflx','area_err' ,
                'hist_flx_err', 'hist_flx_n', 'bin_flx_err', 'area_flx_err',
                'p_mu', 'p_sigma']
    datatable = [mean_tau, bin_stdev, err_stdev, bin_sigma_G, err_sigma_G, 
                 mu_approx, SF, err_SF, err_mu_approx, err_median, 
                 hist_tau , hist_tau_n, bin_tau_cen, hist_delflx , 
                 hist_delflx_n, bin_delflx_cen,  hist_err , hist_err_n, 
                 bin_err_cen, st, area_tau, area_delflx, area_err,
                 hist_flx_err, hist_flx_n, bin_flx_err, area_flx_err,
                 p_mu, p_sigma]
                 
    for label, column in zip(colnames, datatable):
        plot_data[label] = column    
    
    return plot_data

def gauss_stats(tau,y, y_err, flx_err):
    # Define functions for bin statistics 
    rms_robust = lambda x : 0.7414 *(np.percentile(x,75) - np.percentile(x,25))
    rms_std = lambda x : np.std(x)
        
    # Calculate statistics on histograms... 
    sigma_stdev = binned_statistic(tau,y, statistic=rms_std, bins=1)[0][0]
    sigma_G = binned_statistic(tau,y, statistic=rms_robust, bins=1)[0][0]
    
    # Calculate sigma_approx and mu_approx using Fig.5.7 code
    xi = y   # delflx[mask]
    ei = y_err  # delflxerr[mask]
    mu_app, sigma_app = approximate_mu_sigma(xi, ei)
    SF = sigma_app
    mu_mean = np.mean(y)
    
    # fit a gaussian
    # http://stackoverflow.com/questions/7805552/fitting-a-histogram-with-python
    # best fit of data
    (mu_gauss, sigma_gauss) = norm.fit(y)
    (mu_gauss_2, sigma_gauss_2) = norm.fit(flx_err)
    
    print ' Returning statistics on the chosen subsample of the  plot data...'
    stat_data = {}    
    colnames = ['sigma_stdev', 'sigma_G', 'SF', 'mu_app','mu_gauss', 
                'sigma_gauss', 'mu_mean','mu_gauss_2', 'sigma_gauss_2']
    datatable = [sigma_stdev, sigma_G, SF, mu_app, mu_gauss, 
                 sigma_gauss, mu_mean , mu_gauss_2, sigma_gauss_2]
    for label, column in zip(colnames, datatable):
        stat_data[label] = column    
    
    return stat_data  
           
def sf_plot_panels(qso_data,star_data_blue, star_data_red, sample, choice, 
                   nbins, bins_hist, err_factor):  
    '''
    NEW : instead of sf_plotting, this routine is more versatile, as it 
    calls external function  get_plotted_quantities()  to actually 
    calculate the things to be plotted : means of delmag  per bin, etc. 
    
    It plots the four panels, first getting quantities to plot for stars, 
    and then for quasars, and then plotting them altogether . It also 
    creates a separate figure that plots mu_approx on a linear scale 
    '''               
    
    # First calculate plotted quantities 
    
    qso_plot  = get_plotted_quantities(qso_data, nbins, bins_hist,err_factor)
    star_plot = get_plotted_quantities(star_data_blue, nbins, bins_hist,err_factor)
    star_plot1 = get_plotted_quantities(star_data_red, nbins, bins_hist,err_factor)
         

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
    ax3.scatter(np.log10(qso_plot['mean_tau']), qso_plot['SF'], s=p_size, 
                alpha=p_al,c = col1)
    ax3.errorbar(np.log10(qso_plot['mean_tau']), qso_plot['SF'], 
                 qso_plot['err_SF'], linestyle='None',c = col1)
                 
    # blue stars 
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
    ax4.scatter(np.log10(qso_plot['mean_tau']), qso_plot['mu_approx'],
                s=p_size, alpha=p_al,c = col1 )
    ax4.errorbar(np.log10(qso_plot['mean_tau']), qso_plot['mu_approx'], 
                 qso_plot['err_mu_approx'], linestyle='None',c = col1)
                
    # blue stars 
    ax4.scatter(np.log10(star_plot['mean_tau']), star_plot['mu_approx'],
                s=p_size, alpha=p_al,c = col2 )
    ax4.errorbar(np.log10(star_plot['mean_tau']), star_plot['mu_approx'], 
                 star_plot['err_mu_approx'], linestyle='None',c = col2)
                 
    # red stars 
    ax4.scatter(np.log10(star_plot1['mean_tau']), star_plot1['mu_approx'],
                s=p_size, alpha=p_al,c = col3 )
    ax4.errorbar(np.log10(star_plot1['mean_tau']), star_plot1['mu_approx'], 
                 star_plot1['err_mu_approx'], linestyle='None',c = col3)
                 
                 
    ax4.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax4.set_ylim(top=y_mu_top, bottom=y_mu_bott)
    ax4.set_xlim(left=x_left, right=x_right)
    ax4.set_yticks([-0.05,0,0.05])
    ax4.set_yticklabels(['-0.05','0.0', '0.05'])  
    ax4.set_ylabel(r'$\mu_{approx}$', fontsize=20)
    ax4.grid(axis='x')
    ax4.set_xlabel(r'$log_{10} (\Delta _{t})$ [days]',fontsize=20)
    
    title2 = 'SF_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'
    fig1.subplots_adjust(hspace=0)
    plt.savefig(title2)
    plt.show()
    
    # Make a separate  figure for the sigma_G vs delta(t) on a linear scale
    # plotting EXACTLY THE SAME  as for ax4 above, but on a linear x-scale
    # thus I am NOT TAKING np.log10  of x-quantities 
    
    plt.clf()
    fig2 = plt.figure(figsize=(12,4))
    ax1 = fig2.add_subplot(111)
    ax1.scatter(qso_plot['mean_tau'], qso_plot['mu_approx'],
                s=p_size, alpha=p_al,c = col1 )
    ax1.errorbar(qso_plot['mean_tau'], qso_plot['mu_approx'], 
                 qso_plot['err_mu_approx'], linestyle='None',c = col1)
                 
    ax1.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.set_ylim(top=y_mu_top, bottom=y_mu_bott)
    #ax1.set_xlim(left=x_left, right=x_right)
    ax1.set_yticks([-0.05,0,0.05])
    ax1.set_yticklabels(['-0.05','0.0','0.05'])  
    ax1.set_ylabel(r'$\mu_{approx}$', fontsize=20)
    ax1.grid(axis='x')
    ax1.set_xlabel(r'$\Delta _{t}$ [days]',fontsize=20)
    
    title3 = 'Sigma_del_t_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'
    plt.savefig(title3)
    plt.show()
    
    # Make histograms : plot on one subplot histograms of delflx for stars, qso
    # on another subplot, make hist of delflxerr for stars, qso 
    
    plt.clf()
    fig3 = plt.figure(figsize=(12,12))
    
       
    
    
    # Subplot1 : Delta_Mag  (delflx)
    ax1  = fig3.add_subplot(221)  # rows, columns, subplot_id
    # quasars
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
         
    
    
    ax1.set_xlabel(r'$\Delta M$')   
    ax1.set_ylabel(r'$n_{bin} / N_{sample}$') 
    #ax1.set_yticklabels([' '])
    ax1.set_xlim(left=-2, right=2)
    ax1.legend(framealpha=0.7)
    
    # Subplot2 : Delta_Mag_err (delflxerr)
    #ax2 = fig3.add_subplot(222)
    # Quasars
    #ax2.plot(qso_plot['bin_err_cen'], qso_plot['hist_err'], ls='steps', 
    #                  label='Quasars', color='black', lw=lw_h)
                                      
      # blue stars
    #ax2.plot(star_plot['bin_err_cen'], star_plot['hist_err'], ls='steps',
    #          label='Blue stars', color='blue',lw=lw_h)
    
    # red stars 
    #ax2.plot(star_plot1['bin_err_cen'], star_plot1['hist_err'], ls='steps',
    #          label='Red stars', color='red',lw=lw_h)
    
   # ax2.set_xlabel(r'$\sigma(\Delta M): err$')   
    #ax2.set_ylabel(r'$n_{bin} / N_{sample}$') 
    #ax2.set_xlim(left=-2, right=2)
    
    
    # Subplot 2 : plot Delta_Mag / delflxerr 
    
    ax2 = fig3.add_subplot(222)
    # Quasars
    ax2.plot(qso_plot['bin_flx_err'], qso_plot['hist_flx_err'], ls='steps', 
                      label='Quasars', color='black', lw=lw_h)
                                      
    # blue stars
    ax2.plot(star_plot['bin_flx_err'], star_plot['hist_flx_err'], ls='steps',
             label='Blue stars', color='blue',lw=lw_h)
    
    # red stars 
    ax2.plot(star_plot1['bin_flx_err'], star_plot1['hist_flx_err'], ls='steps',
             label='Red stars', color='red',lw=lw_h)
    
    ax2.set_xlabel(r'$\Delta M /  err(\Delta M)$')   
    ax2.set_ylabel(r'$n_{bin} / N_{sample}$') 
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
    
    # Subplot3 : plot only Quasars and gaussian 
    ax3 = fig3.add_subplot(223)
        
    # quasars
    ax3.plot(qso_plot['bin_delflx_cen'],qso_plot['hist_delflx'], ls='steps', 
                      label='Quasars', color='black', lw=lw_h)
    
    # Gaussian 
    bins = qso_plot['bin_delflx_cen']
    st = qso_plot['st']
    mu = st['mu_mean']
    sigma = st['sigma_stdev']
    area = qso_plot['area_delflx']
    y = area*gaussian(bins,mu, sigma)
    #ax3.set_title('Normalized distribution'+r'$  \sigma=%.3f$' %(sigma))
    ax3.set_xlabel(r'$\Delta M$') 
    ax3.plot(bins, y, 'b--', lw=lw_h, label ='Gaussian')     
    ax3.legend(framealpha=0.7)
    ax3.set_ylabel(r'$n_{bin} / N_{sample}$')
    ax3.set_xlim(left=-2, right=2)
    
    # Subplot 4 : plot the CDF of Quasars and Gaussian CDF 
    # http://stackoverflow.com/questions/9378420/how-to-plot-cdf-in-matplotlib-in-python
    
    ax4 = fig3.add_subplot(224)
    c_gauss = np.cumsum(y)
    c_qso = np.cumsum(qso_plot['hist_delflx'])
    ax4.plot(bins, c_gauss, 'b--', lw=lw_h, label='Gauss')
    ax4.plot(bins, c_qso, 'black', lw=lw_h, label='QSO')
    ax4.set_ylabel('CDF')
    ax4.set_xlabel(r'$\Delta M$') 
    ax4.legend(loc=4,framealpha=0.7)
    ax4.set_xlim(left=-2, right=2)
    
    fig3.subplots_adjust(wspace=0.3)    
    
    title = 'Hist_delflx_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'     
    plt.savefig(title)
    plt.show()
    
    
    
    return qso_plot
    
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
                 good_ids_QSO, choice, nbins=200, bins_hist=50, err_factor=1.0):
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
    
    for i in range(1): # len(masterFiles_Q)
        #delflx_Q, tau_Q, err_Q, master_acc_list_Q
        qso_data = add_tau_delflx(masterFiles_Q,inDir_Q, good_ids_Q, i, 
                                  qso_data)
        #delflx_S, tau_S, err_S, master_acc_list_S
        star_data_blue = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_blue, i, 
                                   star_data_blue)
        
        star_data_red = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_red, i, 
                                   star_data_red)                            
                                   
        out = sf_plot_panels(qso_data, star_data_blue, star_data_red,  i, 
                             choice, nbins, bins_hist, err_factor)
    return out, qso_data
    
inDirStars   = 'sf_TRY/sf_stars/'
inDirQSO = 'sf_TRY/sf_qso/'



# Experiment on how different CRTS errors affect SF of both quasars and stars 
# A step forward from plotting just quasars (as below)

# Require  Merr < 0.2 mag for quasars and stars ( big error) : same as was done 
# before 

good_ids_S_blue  = cut_stars(mMax=20, mErrMax = 0.2, gi_Min = -1, gi_Max=1)
good_ids_S_red = cut_stars(mMax=20, mErrMax = 0.2, gi_Min = 1, gi_Max=3)
good_ids_QSO, mask_qso = cut_qso(mErrMax = 0.2 , mMax = 20)

out, qso = plot_both_SF(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO, choice='1.0Eboth0.2', nbins=200, bins_hist=100,
                  err_factor=1.0)
                  

# CHECK CHECK CHECK 
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