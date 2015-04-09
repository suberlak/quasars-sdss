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
    return  qso_id

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


def get_plotted_quantities(data,  nbins):
    '''
    Create input for plot_panels()
    '''    
    delflx = data[0]
    tau = data[1]
    delflxerr = data[2]
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
    
    plot_data = {}
    print ' passing on the  plot data...'
    colnames = ['mean_tau', 'bin_stdev', 'err_stdev', 'bin_sigma_G', 'err_sigma_G',
                'mu_approx', 'SF', 'err_SF', 'err_mu_approx', 'err_median']
    datatable = [mean_tau, bin_stdev, err_stdev, bin_sigma_G, err_sigma_G, 
                 mu_approx, SF, err_SF, err_mu_approx, err_median]
    for label, column in zip(colnames, datatable):
        plot_data[label] = column    
    
    return plot_data
           
def sf_plot_panels(qso_data,star_data_blue, star_data_red, sample, choice, nbins):  
    '''
    Plot the four panels, first getting quantities to plot for stars, and then 
    for quasars, and then plotting them altogether 
    '''               
    
    # First calculate plotted quantities 
    qso_plot  = get_plotted_quantities(qso_data, nbins)
    star_plot = get_plotted_quantities(star_data_blue, nbins)
    star_plot1 = get_plotted_quantities(star_data_red, nbins)
         

    # ####################################
    # Plot  the standard deviation  vs tau 
    # ####################################
    
    # nbins above, check also how many stars are in the sample :
    N_qso = len(np.unique(qso_data[3]))
    N_star_b = len(np.unique(star_data_blue[3]))
    N_star_r = len(np.unique(star_data_red[3]))
    lh_w   = 1.0  # horizontal line thickness 
    lh_st  = '--' # horizontal line style 
    lh_al  = 0.5  # horizontal line alpha parameter 
    p_size = 10
    p_al   = 0.5 
    y_top  = 0.6
    y_bott = -0.05
    x_left = 0.5
    x_right = 3.7
    col1 = 'black'
    col2 = 'blue'
    col3   = 'red'
    
    plt.clf()
    fig1 = plt.figure(figsize=(12,16))

    print 'Plotting Standard Deviation vs Delta t ... '    
    ax1 = fig1.add_subplot(411)
    
    ax1.set_title('SF '+str(N_qso)+' qso, '+str(N_star_b)+' blue, and '+ 
      str(N_star_r)+' red  stars, '+ str(nbins)+  ' bin means', fontsize=20)    
    
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
    ax1.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax1.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4', '0.5'])
    ax1.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.grid(axis='x')
 
    
    # ####################################
    # Plot the Gaussian sigma  vs tau  
    # ####################################
    
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
    ax2.set_ylabel(r'$\sigma_{G} = 0.7413 \cdot (q_{75} - q_{25}) $',fontsize=20)
    ax2.tick_params( axis='x', which='both',  bottom='off', 
                    top='off', labelbottom='off') 
    ax2.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax2.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4', '0.5'])
    ax2.axhline(y=0.1, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
    ax2.axhline(y=0.2, color='black', lw=lh_w, ls=lh_st, alpha=lh_al)
    ax2.grid(axis='x')
    #ax3.set_title('SF '+choice+', '+str(N_objects)+' objects, '+str(nbins)+' bin means')
 
    
    
    # ##############################################
    # Plot the SF from approximate_mu_sigma  vs tau  
    # ##############################################
    
    print ' Plotting SF vs Delta t... ' 
    ax3 = fig1.add_subplot(413)
    
    # qusasars
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
    ax3.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax3.set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4', '0.5'])    
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
    ax4.set_ylim(top=0.6, bottom=-0.2)
    ax4.set_xlim(left=x_left, right=x_right)
    ax4.set_yticks([-0.1,0,0.1,0.2,0.3,0.4,0.5])
    ax4.set_yticklabels(['-0.1','0.0','0.1', '0.2', '0.3', '0.4', '0.5'])  
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
    print 'Hello Im here! '
    fig2 = plt.figure(figsize=(12,4))
    ax1 = fig2.add_subplot(111)
    ax1.scatter(qso_plot['mean_tau'], qso_plot['mu_approx'],
                s=p_size, alpha=p_al,c = col1 )
    ax1.errorbar(qso_plot['mean_tau'], qso_plot['mu_approx'], 
                 qso_plot['err_mu_approx'], linestyle='None',c = col1)
                 
    ax1.axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
    ax1.set_ylim(top=0.6, bottom=-0.2)
    #ax1.set_xlim(left=x_left, right=x_right)
    ax1.set_yticks([-0.1,0,0.1,0.2,0.3,0.4,0.5])
    ax1.set_yticklabels(['-0.1','0.0','0.1', '0.2', '0.3', '0.4', '0.5'])  
    ax1.set_ylabel(r'$\mu_{approx}$', fontsize=20)
    ax1.grid(axis='x')
    ax1.set_xlabel(r'$\Delta _{t}$ [days]',fontsize=20)
    
    title3 = 'Sigma_del_t_'+choice+'_'+str(nbins)+'_bins_'+str(sample)+'.png'
    plt.savefig(title3)
    plt.show()
    
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
        
    for i in range(len(masterFiles)): # len(masterFiles)
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
                 good_ids_QSO):
    inDir_S       = inDirStars
    good_ids_S_blue    = good_ids_S_blue
    good_ids_S_red    = good_ids_S_red
    inDir_Q       = inDirQSO
    good_ids_Q    = good_ids_QSO
    
    nbins = 200 
    
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
        #delflx_Q, tau_Q, err_Q, master_acc_list_Q
        qso_data = add_tau_delflx(masterFiles_Q,inDir_Q, good_ids_Q, i, 
                                  qso_data)
        #delflx_S, tau_S, err_S, master_acc_list_S
        star_data_blue = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_blue, i, 
                                   star_data_blue)
        
        star_data_red = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_red, i, 
                                   star_data_red)                            
                                   
        out = sf_plot_panels(qso_data, star_data_blue, star_data_red,  i, 
                             'both', nbins)
    return out
    
inDirStars   = 'sf_TRY/sf_stars/'
inDirQSO = 'sf_TRY/sf_qso/'


# Plotting stars and quasars on one plot .... 

good_ids_S_blue  = cut_stars(mMax=19, mErrMax = 0.2, gi_Min = -1, gi_Max=1)
good_ids_S_red = cut_stars(mMax=19, mErrMax = 0.2, gi_Min = 1, gi_Max=3)
good_ids_QSO = cut_qso()

out = plot_both_SF(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO)

#plot_qso_stats(qso_cat)

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