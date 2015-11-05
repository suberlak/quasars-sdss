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
import seaborn as sns 
sns.set_context("poster")
from matplotlib import rcParams
rcParams['ytick.labelsize'] = 25
rcParams['xtick.labelsize'] = 25
rcParams['axes.labelsize'] = 25
rcParams['axes.linewidth'] = 3
rcParams['font.size'] = 25
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
     
def get_sigma_mu(xi,ei, return_sigma=False, return_p=False, y_34 = 'mode',  
                 mu_s=100, sig_s=40, sig_lim=[0.00,0.5], mu_lim=[-0.2,0.2],
                 approx=False):
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
        mu_i, sigma_i = approximate_mu_sigma(xi, ei)
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
    """ Calculate statistics for the y, y_err 1D arrays
    
    This function calls get_sigma_mu to calculate mu, sigma, based on the 
    p_sigma and p_mu distributions based on Fig. 5.7 AstroML.  
    Size of sigma  and mu space over which 
    p_sigma and p_mu are calculated is variable, depends on  mu_s, sig_s which 
    define the size of linear spaces, and  mu_lim, sig_lim , which define their
    extent.  From probability distributions both sigma and mu are calculated 
    as mean, mode or expectation value of the probability distribution, 
    depending on the value of y_34 parameter.  
    
    One can also use the approximate way of finding mu, sigma, but by default 
    it is disabled. 
    """     
    y_34    = mu_sig_sample['y_34']
    mu_s    = mu_sig_sample['mu_s']
    sig_s   = mu_sig_sample['sig_s']
    sig_lim = mu_sig_sample['sig_lim']
    mu_lim  = mu_sig_sample['mu_lim']
     
    mu_bin, sigma_bin = get_sigma_mu(y, y_err, y_34, sig_s=sig_s, mu_s=mu_s,
                                 sig_lim = sig_lim, mu_lim = mu_lim)
    #mu_app, sigma_app = approximate_mu_sigma(xi, ei)
    SF = sigma_bin
    #mu_mean = np.mean(y)    
        
    return SF, mu_bin 

       
def get_histogram(xdata, nbins):
    """ Calculate simple histogram given the 1D data array, and number of bins. 
    
    It first calculates unnormalized histogram (because numpy normalization 
    works in  a way so that the integral over all space=1 , which is not what
    I want). Then it divides it by the number of elements in the input data array,
    and the bin width. 
    In case a normalisation was needed, I calculate the area under the histogram,
    to scale up the gaussian by the area. 
    In case a normalised gaussian was needed, I also calculate it, but bins 
    are the same as for not-normalised gaussian . 
    """
    
    N_hist_ttl = float(len(xdata))
    hist, bin_edges = np.histogram(xdata, bins=nbins, density=False) 
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
    return hist, hist_n, bin_cen, area        
    
    
def combine_gaussians(hist_xlim, sigma, mu, xi, ei, err_f=1.0 ):
    """ A function to convolve gaussians according to chosen model
    
    """
    
    xgrid = np.linspace(hist_xlim[0], hist_xlim[1], 10000)
    N = float(len(xi))
           
    
    model2 = np.zeros_like(xgrid)
    for i in range(len(xi)):
        sig_com = np.sqrt(sigma ** 2.0 + (err_f*ei[i])**2.0)
        Gi =  gaussian(xgrid,mu, sig_com)
        model2 += Gi  
    model2= (1.0 / N) * model2   
    
    return xgrid, model2
    
# PLOT THINGS...
# Cols are # delflx  delflxerr  obj_ID  SDSS_r_mMed
    
#starB_gMag = read_file('StarB_g')
#starB_rMag = read_file('StarB')
#starR_gMag = read_file('StarR_g') 
#starR_rMag = read_file('StarR')
 

#fig1 = plt.figure(figsize=(12,12))
#ax1 = fig1.add_subplot(221)
#
## Use SDSS g magnitudes...
##ax1.scatter( starB_gMag['SDSS_g_mMed'], starB_gMag['delflx'] / starB_gMag['delflxerr']  )
#scatter_contour(starB_gMag['SDSS_g_mMed'], starB_gMag['delflx'] / starB_gMag['delflxerr'],
#                threshold=400, log_counts=True, ax=ax1,histogram2d_args=dict(bins=40),
#                plot_args=dict(marker=',', linestyle='none', color='black'))
#ax1.set_xlabel('SDSS g magnitudes')
#ax1.set_ylabel(r'$\chi= \Delta (CRTS \, mag) \, / \,(CRTS \, err)$')
#ax1.set_title('Blue Stars')

# Use SDSS r magnitudes... 

def print_chi(data_dict, axis):
    """ Make the chi plot 
    
    Some more doc if need be ... 
    """
    #plt.clf()
    #fig = plt.figure()
    #axis = fig.add_subplot(111)
    
    mags = data_dict['SDSS_g_mMed'] 
    chi = data_dict['delflx'] / data_dict['delflxerr']
    
    # remove silly values     
    m = ( mags > 17 ) * (mags< 20)
    mags = mags[m];  chi = chi[m]
    
    print('We throw away points with magnitude < 15, so that we are left with %d'%len(chi)+" points")    
    nlev = 15
    scatter_contour(mags, chi, levels=nlev, threshold=50, log_counts=True, ax=axis,
                    histogram2d_args=dict(bins=40), 
                    plot_args=dict(marker=' ', linestyle='none', color='black'),
                    contour_args=dict(cmap=plt.cm.jet))
                    
    step = (max(chi)-min(chi)) / nlev
    print('Colourscale ranges from %f to %f, with 15 levels, i.e. steps every %f '%\
        (min(chi), max(chi), step) ) 
        
    rms_robust = lambda x : 0.7414 *(np.percentile(x,75) - np.percentile(x,25))
    rms_std = lambda x : np.std(x)
   
   # The code below starts at 40 bins and calculates standard deviation sigma, 
   # and robust sigma G for each bin, and finds the smallest number of bins 
   # where there are no NaNs in any bin.
   
    
    nbins1 = 40
    diagnostic = 10
    while  diagnostic > 0 : 
        rms1_all = binned_statistic(mags, chi, statistic=rms_std, bins=nbins1)
        nbins1 = nbins1 - 2
        diagnostic = np.isnan(rms1_all[0]).sum() 
    print('For standard deviation we use %d bins' %nbins1)
    nbins2 = 40
    diagnostic = 10
    while  diagnostic > 0 : 
        rms2_all = binned_statistic(mags, chi, statistic=rms_robust, bins=nbins2)
        nbins2 = nbins2 - 2
        diagnostic = np.isnan(rms1_all[0]).sum()
    print('For robust sigma G we use %d bins' %nbins2)
    
    # Evaluated at the same bins as respective standard deviations... 
    bin_means1 = binned_statistic(mags, chi, statistic = 'mean', bins=nbins1+2)[0]
    bin_means2 = binned_statistic(mags, chi, statistic = 'mean', bins=nbins2+2)[0]
        
        
    chi_rms1 = rms1_all[0]
    mags1 = rms1_all[1][:-1] # number of bins is one more than that of values
    
    chi_rms2 = rms2_all[0]
    mags2 = rms2_all[1][:-1]
    
    #return mags1, mags2, bin_means1, bin_means2, chi_rms1, chi_rms2
    
    sym = 56 # symbol size 
    axis.scatter(mags1, bin_means1 + 2*chi_rms1,s=sym, alpha=1.0, color='yellow', edgecolors='black')
    axis.scatter(mags1, bin_means1 - 2*chi_rms1,s=sym, alpha=1.0, color='yellow', edgecolors='black')
    
    axis.scatter(mags2, bin_means2 + 2*chi_rms2,s=sym, alpha=1.0, color='orange', edgecolors='black' )
    axis.scatter(mags2, bin_means2 - 2*chi_rms2,s=sym, alpha=1.0, color='orange', edgecolors='black' )
          
    axis.set_xlim(xmin=17)
    axis.set_ylim(ymin=-6, ymax=6)
    axis.xaxis.set_label_position('top')
    axis.xaxis.tick_top()
    axis.set_xlabel('SDSS g mag' )
    axis.set_ylabel(r'$\chi $')
    axis.set_xticks([17,18,19,20])
    axis.set_xticklabels(['17','18','19','20'])
 
            
    #plt.savefig('test_test_test.png')
    #plt.show()
    #return mags1, mags2, bin_means1, bin_means2, chi_rms1, chi_rms2

def load_xi_ei(data_dict, mag_min, mag_max):
    """ Load delflx, delflxerr from sample files (star or QSO)
    
    It assumes that the data dictionary has fields delflx, delflxerr,
    SDSS_r_mMed. So it'll not work if we use SDSS_g_mMed, or anything else...
    """
    
    mask = (data_dict['SDSS_g_mMed'] < mag_max) * (data_dict['SDSS_g_mMed'] >mag_min)
    delflx = data_dict['delflx'][mask]
    delflxerr = data_dict['delflxerr'][mask]
    return delflx, delflxerr
    
    
# PLOT  HISTOGRAMS and GAUSSIANS for r magnitudes.... 

# SETTINGS THAT CONTROL THE APPEARANCE OF HISTOGRAM, AND SIGMA, MU LINSPACES 
#mu_sig_sample={}
#colnames = ['bins_hist', 'hist_xlim', 'mu_s', 'sig_s','mu_lim', 'sig_lim','y_34', 'approx']
#datatable = [200, [-1.5,1.5], 40, 100, [-0.007,0.007 ], [0.17, 0.19],'mode', True]

#for label, column in zip(colnames, datatable):
#    mu_sig_sample[label] = column  
    
    
# BLUE Stars   mag: 17-18 
#delflx, delflxerr = load_xi_ei(starB_rMag, 17, 18)

# QSO


def plot_p_distr(delflx, delflxerr):
    ''' A function to find SF (sigma) and mu for subsample of stars or QSO, 
    from p-distributions. One essentially manually adjusts mu_lim, sig_lim, 
    until the distribution contains the peak with good sampling (increasing
    the sigma or mu space does not really work because the system easily 
    runs out of memory when mu has 200 points and sig 80...
    
    '''

    mu_lim = [-0.004,-0.001]
    sig_lim = [0.149, 0.153]
    
    mu,sigma = get_sigma_mu(xi=delflx, ei=delflxerr, sig_lim=sig_lim, 
                            mu_lim=mu_lim, return_sigma=True)
    mu_bin, sig_bin, p_mu, p_sigma = get_sigma_mu(xi=delflx, ei=delflxerr, 
                            sig_lim=sig_lim, mu_lim=mu_lim, return_p=True)
    fig = plt.figure()
    plt.clf()
   
    ax1 = fig.add_subplot(211)
    ax1.plot(mu, p_mu, '-o')
    ax1.set_title('p(mu)')
    ax1.scatter(mu_bin, p_mu[np.abs(mu-mu_bin) < 0.00001], marker='o', c='y')
    
    ax2 = fig.add_subplot(212)
    ax2.plot(sigma, p_sigma, '-o')
    ax2.set_title('p(sigma)')
    fig.tight_layout()
    plt.savefig('p_sigma_p_mu_QSO_18.5-19.png')
    plt.show()    

    print 'Calculating mode of p_sigma , p_mu, we get sigma=',sig_bin, \
    ' mu=', mu_bin
    
    return mu, sigma, p_mu, p_sigma, mu_bin, sig_bin

#out= plot_p_distr(delflx, delflxerr)
 #SF, mu = get_stats(delflx, delflxerr, mu_sig_sample) 

def do_panel_plot(data, mag_min, mag_max, sig, mu, hist_lim, filename, err_factor):
    ''' A short loop that loops over the four panels for QSO / Blue Stars: 
    in the upper left corner it makes chi plot, and in the following 
    three panels (clockwise) it plots histograms of QSOs in the 
    sample, and overplots two models: one with error_factor = 1.0, 
    and another with error_factor = 1.3 , where 
    model = SUM(Gaussians(mu, sqrt(SF^2+(err*err_f)^2)), and is more
    defined in combined_gaussians() function.
    
    Accepts:
    - data, which is  a dictionary with keys ['SDSS_r_mMed', 'delflx', 
    'delflxerr', 'obj_ID'], 
    - mag_min , mag_max : being lists of lower and upper magnitude limits in each
    subsample of the three plots 
    - sig, mu, which are lists of previously-computed mu and sigma 
    from p-distributions for each subsample using  plot_p_distr()
    - hist: limits on the histogram 
    '''
    # do the QSO plot 
    plt.clf()
    fig, axs = plt.subplots(2,2, figsize=(12, 12), facecolor='w', edgecolor='k')
    #fig.subplots_adjust(hspace = .5, wspace=.5)
    axs = axs.ravel()
    # from http://stackoverflow.com/questions/17210646/python-subplot-within-a-loop-first-panel-appears-in-wrong-position
    
    # First plot the chi distr in the upper left corner 
   
    print_chi(data_dict=data, axis=axs[0])
    
    # Then loop over different ranges and plot the models... 
    for i in range(len(mag_min)):
        
        delflx, delflxerr = load_xi_ei(data, mag_min[i], mag_max[i])
        len_original = float(len(delflx))
        mask = (delflx < hist_lim[1] )*(delflx>hist_lim[0])
        delflx = delflx[mask]
        delflxerr = delflxerr[mask]
        percent = (100.0*float(len(delflx))) / len_original
        print 'The hist limits criteria are satisfied by ', percent, ' percent'
        
        hist, hist_n, bins, area = get_histogram(xdata=delflx, nbins=150)
        
        # first model : SF^2 + ( 1.0 * err) ^2
        #xgrid, model1 = combine_gaussians(xi=delflx, ei=delflxerr, 
         #                            hist_xlim=hist_lim, sigma=sig[i], 
         #                            mu=mu[i], err_f=0 )
                                     
        # second model : SF^2 + ( 1.3 * err) ^2
        xgrid, model2 = combine_gaussians(xi=delflx, ei=delflxerr, 
                                     hist_xlim=hist_lim, sigma=sig[i], 
                                     mu=mu[i], err_f=err_factor[i] )                             
                                     
        # Plot the histogram    
        axs[i+1].plot(bins, hist,ls='steps', lw=2, label='data', color='black')
        axs[i+1].set_xlabel(r'$\Delta m$')
        axs[i+1].set_ylabel(r'$\mathcal{N}$')
        # Plot  model1
        #model1_lab = r'$ f_{c}=1.0$' # % (sig[i], mu[i])
        #axs[i+1].plot(xgrid, model1, lw=2, label=model1_lab )
        
        # Plot model2
        model2_lab = r'$ f_{c}=%.2f$' %(err_factor[i]) #% (sig[i], mu[i])
        axs[i+1].plot(xgrid, model2, lw=3, label=model2_lab, ls='--', color='red' )
    
        # Set x limits
        axs[i+1].set_xlim(xmin=hist_lim[0], xmax=hist_lim[1])
        axs[i+1].set_ylim(ymin=-0.05)
        
        # Set x ticks
        axs[i+1].set_xticks([-0.5,0,0.5])
        axs[i+1].set_xticklabels(['-0.5', '0', '0.5'])
   
        
        # Get y limits  
        ymin, ymax = axs[i+1].get_ylim()
        # Add text saying the limits ... 
        axs[i+1].text(0.1, 0.75*ymax, r'$%.1f - %.1f$'%(mag_min[i], mag_max[i]))     
        
        # For the 17 < m < 18 plot (upper-right), and 
        # For the 18.5 < m < 19 plot (lower-righ), 
        # move ticklabels to the right
        
        if i==0 or i==2 : 
            axs[i+1].yaxis.tick_right()             
            axs[i+1].yaxis.set_label_position('right')
            axs[i+1].yaxis.set_ticks_position('both')
            
        # For 17 - 18 plot move the xlabels to the top
            
        if i==0 : 
            axs[i+1].xaxis.set_label_position('top')
            axs[i+1].xaxis.tick_top()
            
    #lgd =  plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))   
    fig.subplots_adjust(hspace=0.1, wspace=0.1)       
    plt.savefig(filename+'_sample_panels.png') #,
#                bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()
    
    




# What are the ranges of magnitude on each subplot 
# from mag_min[i] to mag_max[i]
mag_min = [17, 18, 18.5]
mag_max = [18, 18.5, 19]

# What are sigma (SF), mu, histogram limits (x-axis) for QSO and stars ? 
# they are different for QSO's and stars 

qso_sig = [0.0997, 0.10456, 0.15105]
qso_mu = [-0.00097, -0.00033, -0.00254]
starB_sig = [0,0,0, 0]
starB_mu = [0,0, 0]
hist_qso = [-0.6,0.6]
hist_starB = [-0.6,0.6]

# Set the error factor for the plots ...
#err_inc = 0.72

err_factor = [0.72, 0.91, 1.07]
star_sigma_list = [0,0,0]
qso_sigma_list = [0.045, 0.052, 0.07]


# Read files 
QSO = read_file('QSO')
starB_rMag = read_file('StarB_g')


#do_panel_plot(QSO, mag_min, mag_max, qso_sig, qso_mu, hist_qso, 'new_QSO', err_inc)
#do_panel_plot(QSO,mag_min, mag_max, qso_sigma_list, qso_mu, 
#              hist_starB, 'NEW_QSO', err_factor )
do_panel_plot(starB_rMag,mag_min, mag_max, star_sigma_list, starB_mu, 
              hist_starB, 'NEW_StarB', err_factor )