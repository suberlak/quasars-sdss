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
from matplotlib import rcParams
from matplotlib.patches import Rectangle
rcParams['ytick.labelsize'] = 25
rcParams['xtick.labelsize'] = 25
rcParams['axes.labelsize'] = 35
rcParams['axes.linewidth'] = 3
rcParams['font.size'] = 25
#rcParams.update({'figure.subplot.hspace' : 0})

rcParams.update({'figure.autolayout': False})

from scipy.optimize import curve_fit

from scipy.stats import binned_statistic
from astroML.stats import median_sigmaG
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})
#import seaborn as sns
#sns.set_context("poster")
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

   
def model_sf(t, sf_inf=0.25, tau = 1.0):
    '''
    A fiducial Structure Function model 
    based on the DRW 
    
    Parameters:
    -----------
    t : delta_t  points at which the model should be evaluated
        (along x axis )
    tau      : tau (in days)
    sf_inf :sf_infinity , no units
    
    Returns:
    --------
    sf : structure function 
    
    '''
    
    br = 1.0-np.exp(-t/tau)
    sf = sf_inf * np.power(br,0.5)
    return sf

def split_bin(array, length):
    ''' Split an input array into chunks of desired length
    The chunks will all have the same length, the part that
    is left (non-divisible by length) will form the last chunk'''
    
    length = max(1, length)
    return [array[i:i + length] for i in range(0, len(array), length)]

def get_plotted_quantities(data, nbins, mu_sig_sample, mu_sig_generic, err_factor=1.0,save_bin = False ):
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
    
    #####
    ##### Panel 3 : SF  , Panel 4 : mu_approx   
    #####
    
    N_bin = [] 
     
     
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
            
            # store the number of points per bin
            N_bin.append(len(xi)) 
            
            if N == 25 and save_bin == True : 
                DATA = np.column_stack((xi,ei))
                outfile = 'QSO_bin_'+str(N)+'_xi_ei.txt'
                np.savetxt(outfile, DATA, delimiter= ' ', fmt = "%s")            
                return xi  
                
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
            
            # store the number of points per bin
            N_bin.append(len(xi)) 
            
            
            if save_bin == True : 
                DATA = np.column_stack((xi,ei))
                outfile = 'QSO_bin_'+str(N)+'_xi_ei.txt'
                np.savetxt(outfile, DATA, delimiter= ' ', fmt = "%s")            
                return xi           
                
            # Split into chunks, calculate mu, sigma per chunk and average 
            chunk_size = 40000 
            if len(xi) > chunk_size : 
                split_xi = split_bin(xi,chunk_size)
                split_ei = split_bin(ei, chunk_size)
                mu_chunks, sigma_chunks = [],[]
                
                for i in range(len(split_xi)):
                    mu, sigma = get_sigma_mu(split_xi[i],split_ei[i],approx, y_34)
                    print('Indiv : mu_chunk =%f and sigma_chunk=%f'%(mu, sigma))
                    mu_chunks.append(mu) , sigma_chunks.append(sigma)
                mu_bin, sigma_bin = np.average(mu_chunks), np.average(sigma_chunks)    
  
            else:
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
                'mu_plot', 'SF', 'err_SF', 'err_mu_plot', 'err_median', 'N_bin']
    datatable = [mean_tau, bin_stdev, err_stdev, bin_sigma_G, err_sigma_G, 
                 mu_plot, SF, err_SF, err_mu_plot, err_median, N_bin]
                 
    for label, column in zip(colnames, datatable):
        plot_data[label] = column    
    
    return plot_data


def sf_plot_panels(qso_data,star_data_blue, star_data_red, sample, choice, 
                   nbins,  err_factor, approx, y_34,sf_panel_only,save_bin,
                   multipanel, bin_hist_info):  
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
        
    qso_plot  = get_plotted_quantities(qso_data, nbins, mu_sig_sample, mu_sig_generic, err_factor, save_bin )  
    #print qso_plot.keys()
    #return qso_plot
    if save_bin == True : 
        return qso_plot
    
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
     
    ##############################################
    #    PLOT SETTINGS COMMON FOR ALL PLOTS 
    ##############################################
       
        
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

    # y limits for mu  
    y_mu_top = 0.1
    y_mu_bott = -0.1

    # x limits for ALL PANELS 
    x_left = 0.5
    x_right = 3.7
    
    # colors for quasars, blue and red stars 
    col1 = 'black'
    col2 = 'blue'
    col3   = 'red'
    
    # print some info about what is plotted below 
    print('The general panel info ... :') 
    print('SF '+str(N_qso)+' qso, '+str(N_star_b)+' blue, and '+ 
      str(N_star_r)+' red  stars, '+ str(nbins)+  ' bin means')
      
    ##########################################
    ######### BIN HIST INFO ##################
    ##########################################
    if bin_hist_info == True :  
        plt.clf()
        fig, ax = plt.subplots(3,1,  figsize=(12,12))
        x = np.arange(nbins)
        ax[0].plot(x, qso_plot['N_bin'], ls='steps', color='black')
        ax[1].plot(x, star_plot['N_bin'], ls='steps', color='blue')
        ax[2].plot(x, star_plot1['N_bin'], ls='steps', color='red')  
        title = 'N_points_per_bin_qso_stars_SF_panel_'+str(nbins)+'_'+str(sample)+'.png'    
        plt.savefig(title)
        plt.show()
        return qso_plot
    ##########################################
    ############ SF PANEL ONLY ###############
    ##########################################   
   
    if sf_panel_only == True : 
        
        plt.clf()
        fig2 = plt.figure(figsize=(12,5))
        
        print ' Plotting SF vs Delta t... ' 
        ax = fig2.add_subplot(111)
            
        # Plot fiducial DRW  
        xdata = qso_plot['mean_tau']
        sf = qso_plot['SF']    
        popt, pcov = curve_fit(model_sf, xdata, sf)
        y = model_sf(xdata, sf_inf=popt[0], tau = popt[1]) # tau 1 year in days 
        
        ax.plot(np.log10(xdata), y , lw=3, c = 'orange', ls='--')
        text = r'$ \mathrm{Model:}\ \tau=%.3f \, \mathrm{days} \, , \ SF_{\infty}=%.3f \, \mathrm{mags}$'% \
                 (popt[1],popt[0])
        ax.text(x=0.75, y=0.3,s = text )
        
        # quasars
        ax.scatter(np.log10(qso_plot['mean_tau']), qso_plot['SF'], s=p_size, 
                    alpha=p_al,c = col1)
        ax.errorbar(np.log10(qso_plot['mean_tau']), qso_plot['SF'], 
                     qso_plot['err_SF'], linestyle='None',c = col1)
                     
        # blue stars 
        ax.scatter(np.log10(star_plot['mean_tau']), star_plot['SF'], s=p_size, 
                    alpha=p_al,c = col2)
        ax.errorbar(np.log10(star_plot['mean_tau']), star_plot['SF'], 
                     star_plot['err_SF'], linestyle='None',c = col2)
                     
        # red stars
        ax.scatter(np.log10(star_plot1['mean_tau']), star_plot1['SF'], s=p_size, 
                    alpha=p_al,c = col3)
        ax.errorbar(np.log10(star_plot1['mean_tau']), star_plot1['SF'], 
                     star_plot1['err_SF'], linestyle='None',c = col3)
                     
        ax.set_ylim(bottom=y_bott, top=y_top)
        ax.set_xlim(left=x_left, right=x_right)
        ax.set_ylabel(r'$SF $')
        ax.tick_params( axis='x', which='both',  bottom='on', 
                        top='off', labelbottom='on')
        ax.grid(axis='x')
        ax.set_yticks([0,0.1,0.2,0.3,0.4])
        ax.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
        ax.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)    
        ax.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
        ax.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al) 
        ax.set_xlabel(r'$log_{10} (\Delta _{t})$ [days]')
        
        # draw a rectangle ... 
#        someX, someY = 0.52, 0.02
#        width, height = 1.2, 0.1
#        ax.add_patch(Rectangle((someX, someY),width , height, edgecolor='red', alpha=0.9, fill=False, lw=4, ls='solid' ))
#        
        title = 'SF_panel'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'    
        plt.tight_layout()
        plt.savefig(title)
        plt.show()
    
        
        return qso_plot
   
   
   
    ##########################################
    ############ MULTIPANEL ##################
    ##########################################   
   
    if multipanel == True : 
        plt.clf()
        fig1 = plt.figure(figsize=(12,16)) 
        
        ##########################################
        ############ Panel 1 #####################
        ##########################################
        
        print 'Plotting Standard Deviation vs Delta t ... '    
        ax1 = fig1.add_subplot(411)
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
        
             
        ax1.set_ylabel(r'$\sigma$')  
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
        # ax1.spines['top'].set_linewidth(0.5)
     
        
        ##########################################
        ############ Panel 2 #####################
        ##########################################
        
        print 'Plotting sigma_Gaussian  vs Delta t ... '
        ax2 = fig1.add_subplot(412)
        
            
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
        ax2.set_ylabel(r'$\sigma_{G} $')
        ax2.tick_params( axis='x', which='both',  bottom='off', 
                        top='off', labelbottom='off') 
        ax2.set_yticks([0,0.1,0.2,0.3,0.4])
        ax2.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
        ax2.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
        ax2.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
        ax2.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
        ax2.grid(axis='x')
          
        
        ##########################################
        ############ Panel 3 #####################
        ##########################################
        
        print ' Plotting SF vs Delta t... ' 
        ax3 = fig1.add_subplot(413)
            
        # Plot fiducial DRW
        # Fitting to QSOs , and y-data from panel 3 (SF)
        # http://stackoverflow.com/questions/26058792/correct-fitting-with-scipy-curve-fit-including-errors-in-x
            
        xdata = qso_plot['mean_tau']
        sf = qso_plot['SF']    
        popt, pcov = curve_fit(model_sf, xdata, sf)
        y = model_sf(xdata, sf_inf=popt[0], tau = popt[1]) # tau 1 year in days 
        
        # -->  plot the model sigma with error... perhaps use error envelope ?  
        
        #err_sig = qso_plot['err_SF']
        #sf_folded = np.sqrt((y ** 2.0)+ (err_sig ** 2.0) )
        
        ax3.plot(np.log10(xdata), y , lw=3, c = 'orange', ls='--')
        text = r'$ \mathrm{Model:}\ \tau=%.3f \, \mathrm{days} \, , \ SF_{\infty}=%.3f \, \mathrm{mags}$'\
            %(popt[1],popt[0])
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
        ax3.set_ylabel(r'$SF $')
        ax3.tick_params( axis='x', which='both',  bottom='off', 
                        top='off', labelbottom='off')
        ax3.grid(axis='x')
        ax3.set_yticks([0,0.1,0.2,0.3,0.4])
        ax3.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
        ax3.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)    
        ax3.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
        ax3.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al) 
        
         
        # draw a rectangle ... 
#        someX, someY = 0.52, 0.02
#        width, height = 1.2, 0.1
#        ax3.add_patch(Rectangle((someX, someY),width , height, edgecolor='red', alpha=0.9, fill=False, lw=4, ls='solid' ))
#            
            
        ##########################################
        ############ Panel 4 #####################
        ##########################################
        
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
        ax4.set_ylabel(r'$\mu$')
        ax4.grid(axis='x')
        ax4.set_xlabel(r'$log_{10} (\Delta _{t})$ [days]')
        
        fig1.subplots_adjust(hspace=0)
        
        title2 = 'SF_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'    
        #plt.tight_layout()
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
    
def plot_panels(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                 good_ids_QSO, choice, nbins=200,  err_factor=1.0,
                 approx=True, y_34 = 'mode',sf_panel_only=False, save_bin=False,
                 multipanel=True, bin_hist_info=False):
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
    
    for i in range(len(masterFiles_Q)): #  
        qso_data = add_tau_delflx(masterFiles_Q,inDir_Q, good_ids_Q, i, 
                                  qso_data)

        star_data_blue = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_blue, i, 
                                   star_data_blue)
        
        star_data_red = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_red, i, 
                                   star_data_red)                            
                                   
        out = sf_plot_panels(qso_data, star_data_blue, star_data_red,  i, 
                             choice, nbins, err_factor, approx, y_34,sf_panel_only, 
                             save_bin,multipanel, bin_hist_info)
    
    return out, qso_data, star_data_blue, star_data_red
    
inDirStars   = 'sf_TRY/sf_stars/'
inDirQSO = 'sf_TRY/sf_qso/'

# Experiment on how different CRTS errors affect SF of both quasars and stars 
# A step forward from plotting just quasars (as below)

# Require  Merr < 0.2 mag for quasars and stars ( big error) : same as was done 
# before 

#
#  Standard run, selecting all objects with SDSS r_mMag < 20  
#
good_ids_S_blue  = cut_stars(mMax=19, mErrMax = 0.3, gi_Min = -1, gi_Max=1)
good_ids_S_red = cut_stars(mMax=19, mErrMax = 0.3, gi_Min = 1, gi_Max=3)
good_ids_QSO, mask_qso = cut_qso(mErrMax = 0.3 , mMax = 19)

out, qso, star_b, star_r = plot_panels(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO, choice='NEW_chunks_1.0E_err_0.3_mode_', nbins=200, 
                  err_factor=1.0, approx=False, y_34 = 'mode',sf_panel_only=True, save_bin=False,
                  multipanel=False, bin_hist_info=False)
     
#
# do an experiment choosing only certain ranges for the plot, i.e. 
# making narrower SDSS mag cuts      
#   
#
     
mMin = [17,18,18.5]
mMax = [18,18.5,19]
#
#for i in range(len(mMin)):
#    Min = mMin[i]
#    Max = mMax[i]
#    print('Using now only lightcurves with SDSS  %f< r_mMed < %f' % (Min, Max))
#    
#    good_ids_S_blue  = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1, gi_Max=1)
#    good_ids_S_red = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, gi_Max=3)
#    good_ids_QSO, mask_qso = cut_qso(mMin = Min, mMax=Max, mErrMax = 0.3)
#    
#    ch = '_1.0E_err_0.3_approx_mag_'+str(Min)+'-'+str(Max)+'_'
#    
#    out, qso, star_b, star_r = plot_panels(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
#                  good_ids_QSO, choice=ch, nbins=200, 
#                  err_factor=1.0, approx=True, y_34 = 'mode',sf_panel_only=True, save_bin=False,
#                  multipanel=False, bin_hist_info=False)