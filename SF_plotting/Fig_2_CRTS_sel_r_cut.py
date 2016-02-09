# -*- coding: iso-8859-1 -*-
import os
import numpy as np 
import matplotlib.pyplot as plt
from astroML.stats import median_sigmaG
from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
import CRTS_paper_modules as mod

''' A program to make Fig.2 for the CRTS paper describing the error problems
with CRTS. A notebook explaining this code is Fig_2_CRTS_sel_r_cut.ipynb,  
and any new development is first checked there.  

To change the magnitude range for the four-panels made here without any correction, 
in an approximate way, change Min and Max parameters in line 77 

'''

# Set a directory to save the sample...
outDir = os.path.join(os.getcwd(), 'Fig_2_data', datetime.datetime.now().strftime('%Y-%m-%d')+ '/')

cols1, qso_cat = mod.get_qso_catalog(catalog='DB_QSO') 
cols2 , star_cat= mod.get_stars_catalog()




Min = 18.5
Max = 19
magnitudes = ['r']  # don't need g magnitude , decided to use only r magnitude 

objects_in_cut = {}

for mag in magnitudes : 
    cut_mag = mag
    report_mag = mag
    
    print('\nUsing now only lightcurves with SDSS  %f< %s < %f' % (Min, cut_mag, Max))
    print('\n Reporting SDSS %s  '% report_mag)

    good_ids_S_blue = mod.cut_stars(star_cat = star_cat,mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1,
                                              gi_Max=1, cut_mag=cut_mag + '_mMed', 
                                              report_mag=report_mag + '_mMed')
    
    good_ids_S_red = mod.cut_stars(star_cat = star_cat,mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, 
                                           gi_Max=3, cut_mag=cut_mag + '_mMed', 
                                           report_mag=report_mag + '_mMed')
    
    good_ids_QSO = mod.cut_qso(qso_cat=qso_cat, mMin = Min, mMax=Max, mErrMax = 0.3, 
                                               cut_mag=cut_mag,report_mag=report_mag)
    objects_in_cut[mag] = {'starsB':good_ids_S_blue, 'starsR':good_ids_S_red, 
                           'qso':good_ids_QSO}
    

bins = {}
bin_types = ['r_cut']  # removed 'g_cut'

objects = objects_in_cut['r'].keys()

# first need to explicitly initialize the dictionaries 
for b in bin_types:
    bins[b] = {}
    
for obj in objects : 
    bins['r_cut'][obj] =  objects_in_cut['r'][obj]


# inside the main loop : get tau, delflx from a master file, either qso or star


inDirStars   = 'sf_file_per_LC/star/'
inDirQSO = 'sf_file_per_LC/qso/'

out_dic = {}

#for b in bins.keys():
# read in only r_cut 

b = 'r_cut'
print 'Reading in xi, ei for bin ', b
out_dic[b] = {}   # initialize the dic 

good_ids_S_blue = bins[b]['starsB']
good_ids_S_red = bins[b]['starsR']
good_ids_QSO = bins[b]['qso']

qso, starB, starR = mod.read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
              good_ids_QSO,xi_ei_data=None, fc=None)

# put into a dictionary : makes it more explicit 
out_dic[b] = {'starsB': starB, 'starsR': starR, 'qso':qso}


# Straight after reading-in xi, ei,   one can proceed directly to part 9) (one bin) or 10 : all bins sigma comparison 
# or to Saving just the log(tau) samples of xi, tau, ei. 


# Run the calculation over all bins : 
# - keep the option to plot a bin if one wants to compare the three methods : need to provide object N
#   
# - as it is, calculate three values for sigma, etc., and plot it  

# Set correction factor
fc = 1.0
# Set the bin  : here we use only r_cut 
b = 'r_cut'


# Loop over objects in the dictionary  
for obj in out_dic[b].keys() : 
    print '\n For obj = ', obj
    # obj = 'qso'  # or starsB,  starsR 
    m_ij = out_dic[b][obj][0]
    tau =  out_dic[b][obj][1]
    e_ij = fc * out_dic[b][obj][2]

    #n =  out_dic[b][obj][3]

    #def calculate_mu_sig_single_bin(m_ij,e_ij,tau, N=1, plotLnL= True):

    nbins = 200 

    # Pull out some tau to plot means : common to all panels 
    binned_tau = binned_statistic(tau, tau, statistic='mean', bins=nbins)
    mean_tau = binned_tau[0]
    # Take N from each bin... 'count' function works like a regular histogram
    binned_count = binned_statistic(tau, tau, statistic='count', bins=nbins)
    bin_count = binned_count[0]
    #bin_names = np.arange(1,len(binned_count[2]))

     # Calculate median preprocessed photometric error per bin 
    binned_err_median = binned_statistic(tau, e_ij, statistic='median', bins=nbins) 
    err_median = binned_err_median[0]

    # checking for empty bins : either mean or some custom function, but not
    # count! If statistic='count', then check for 0's , and not for nan's/ 
    non_empty_bins = np.bitwise_not(np.isnan(mean_tau))

    # reassign number of points in a bin and  tau position 

    bin_count = bin_count[non_empty_bins]
    mean_tau = mean_tau[non_empty_bins]
    err_median = err_median[non_empty_bins]

    # Which point belongs to which bin
    bin_number  = binned_tau[2]
    
    ####
    ####  Panel 1 : Standard Deviation 
    ####

    rms_std = lambda x : np.std(x)
    stdev_binned = binned_statistic(tau, m_ij, statistic = rms_std, 
                                              bins=nbins)


    bin_stdev = stdev_binned[0][non_empty_bins]  
    #bin_number = stdev_binned[2]  
     # since each point belongs to some bin : len(bin_number) =len(delflx)


    # error on standard deviation in the bin     
    err_stdev = bin_stdev / np.sqrt(2.0*(bin_count - 1.0))

    #####
    ##### Panel 2  : Gaussian rms  
    #####
    rms_robust = lambda x : 0.7414 *(np.percentile(x,75) - np.percentile(x,25))
    bin_sigma_G = binned_statistic(tau, m_ij, statistic = rms_robust, 
                                      bins=nbins)[0][non_empty_bins]

    # error on Gaussian estimate of rms in the bin 
    err_sigma_G = bin_sigma_G* 1.06 / np.sqrt(bin_count)


    #####
    ##### Panel 3 (SF)   and Panel 4   (mu)
    #####


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

    # Loop over all bins  calculating approximate mu and sigma 

    mu_bins = {}
    sig_bins = {}

    sig_bins['approx'] = np.zeros(200)
    mu_bins['approx'] = np.zeros(200)

    for N in np.unique(bin_number):
        print('\r --- Calculating mu, sigma for bin %d' % N),
        xi = m_ij[bin_number == N]
        ei = e_ij[bin_number == N]

        # 1) Calculate in an approximate way 
        mu_approx, sig_approx = approximate_mu_sigma(xi, ei)
        sig_bins['approx'][N-1] = sig_approx
        mu_bins['approx'][N-1] = mu_approx

    # Calculate error of points 
    sig_bins['approx_err'] = sig_bins['approx']* 1.06 / np.sqrt(bin_count)
    mu_bins['approx_err'] = bin_stdev / np.sqrt(bin_count)


    print '\n In this calculation, fc=', fc 

    
    # Save full results for plotting 4 panels 
    # Save the results of calculation  : Panel 1,2,3,4 

    fname = outDir+ b+'_'+str(Min)+'-'+str(Max)+'_'+obj+'_fc-'+str(fc)+'_mean_tau_sig_sigG_SF_mu.txt'
    print 'Results saved as ', fname

    data = np.column_stack((mean_tau, bin_stdev, err_stdev, bin_sigma_G, err_sigma_G, 
                            sig_bins['approx'], sig_bins['approx_err'],
                            mu_bins['approx'],mu_bins['approx_err'] ))

    header = 'mean_tau   bin_stdev   err_stdev  bin_sigma_G  err_sigma_G   SF  SF_err  mu  mu_err  '
    np.savetxt(fname, data, fmt = '%s', delimiter = ' ' , header=header )



# Read the results 

# Initialize the dictionary 

plot_data = {}
for obj in out_dic[b].keys() : 
    plot_data[obj] = {}
    
    fname = outDir+ b+'_'+str(Min)+'-'+str(Max)+'_'+obj+'_fc-'+str(fc)+'_mean_tau_sig_sigG_SF_mu.txt'
    print 'Reading from', fname

    colnames = open(fname,'r').read().splitlines()[0][1:].split()
    d = np.genfromtxt(fname, dtype=float)

    for label, column in zip(colnames, d.T):
        plot_data[obj][label] = column


import matplotlib 
from matplotlib import rcParams
matplotlib.rcParams['font.size'] = 17

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
colors= {}
colors['qso'] = 'black'
colors['starsB'] = 'blue'
colors['starsR']   = 'red'

fig,ax = plt.subplots(4,1, figsize=(8,12), sharex=True)
fig.subplots_adjust(hspace=0)

# Panel 1 
for obj in plot_data.keys() : 
    ax[0].scatter(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['bin_stdev'], s=p_size, 
                    alpha=p_al, c = colors[obj])
    ax[0].errorbar(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['bin_stdev'],plot_data[obj]['err_stdev'], 
                   linestyle='None', c = colors[obj]  )

ax[0].set_ylabel(r'$\sigma_{stdev}$',fontsize=20)  
ax[0].tick_params( axis='x', which='both',  bottom='off', 
                top='off', labelbottom='off') 
ax[0].set_ylim(bottom=y_bott, top=y_top)
ax[0].set_xlim(left=x_left, right=x_right)
ax[0].set_yticks([0,0.1,0.2,0.3,0.4])
ax[0].set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
ax[0].axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[0].axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[0].axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[0].grid(axis='x')


# Panel 2
for obj in plot_data.keys() : 
    ax[1].scatter(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['bin_sigma_G'], s=p_size, 
                    alpha=p_al, c = colors[obj])
    ax[1].errorbar(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['bin_sigma_G'],plot_data[obj]['err_sigma_G'], 
                   linestyle='None', c = colors[obj] )


ax[1].set_ylabel(r'$\sigma_{G}$',fontsize=20)  
ax[1].tick_params( axis='x', which='both',  bottom='off', 
                top='off', labelbottom='off') 
ax[1].set_ylim(bottom=y_bott, top=y_top)
ax[1].set_xlim(left=x_left, right=x_right)
ax[1].set_yticks([0,0.1,0.2,0.3,0.4])
ax[1].set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
ax[1].axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[1].axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[1].axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[1].grid(axis='x')


# Panel 3
 
def model_sf(t, sf_inf=0.25, tau = 1.0):
    br = 1.0-np.exp(-t/tau)
    sf = sf_inf * np.power(br,0.5)
    return sf

# Calculate the model DRW fit for QSO
obj = 'qso'
xdata = plot_data[obj]['mean_tau']
sf = plot_data[obj]['SF']
popt, pcov = curve_fit(model_sf, xdata, sf)
y = model_sf(xdata, sf_inf=popt[0], tau = popt[1]) # tau 1 year in days 

# Fold-in the error to the model SF , plot 
# both folded and not-folded version \
err_sig = plot_data[obj]['SF_err']
y_fold = np.sqrt((y ** 2.0)+ (err_sig ** 2.0) )
ax[2].plot(np.log10(xdata), y_fold , lw=3, c = 'green', ls='--')
ax[2].plot(np.log10(xdata), y , lw=3, c = 'orange', ls='--')

text = r'$ \mathrm{Model:}\ \tau=%.3f \,\mathrm{days} , \ SF_{\infty}=%.3f \,\mathrm{mag}$'%(popt[1],popt[0])
ax[2].text(x=1.0, y=0.3,s = text )

    
for obj in plot_data.keys() : 
    ax[2].scatter(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['SF'], s=p_size, 
                    alpha=p_al, c = colors[obj])
    ax[2].errorbar(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['SF'],plot_data[obj]['SF_err'], 
                   linestyle='None', c = colors[obj]  )

ax[2].set_ylim(bottom=y_bott, top=y_top)
ax[2].set_xlim(left=x_left, right=x_right)
ax[2].set_ylabel(r'$SF $',fontsize=20)
ax[2].tick_params( axis='x', which='both',  bottom='off', 
                top='off', labelbottom='off')
ax[2].grid(axis='x')
ax[2].set_yticks([0,0.1,0.2,0.3,0.4])
ax[2].set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
ax[2].axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)    
ax[2].axhline(y=0.1, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[2].axhline(y=0.2, color='black', lw=lh_w, ls=lh_st,alpha=lh_al) 


# Panel 4 
for obj in plot_data.keys() : 
    ax[3].scatter(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['mu'], s=p_size, 
                    alpha=p_al, c = colors[obj])
    ax[3].errorbar(np.log10(plot_data[obj]['mean_tau']), plot_data[obj]['mu'],plot_data[obj]['mu_err'], 
                   linestyle='None', c = colors[obj]  )


ax[3].axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[3].set_ylim(top=y_mu_top, bottom=y_mu_bott)
ax[3].set_xlim(left=x_left, right=x_right)
ax[3].set_yticks([-0.05,0,0.05])
ax[3].set_yticklabels(['-0.05','0.0', '0.05'])  
ax[3].set_ylabel(r'$\mu$', fontsize=20)
ax[3].grid(axis='x')
ax[3].set_xlabel(r'$log_{10} (\Delta _{t})$ [days]',fontsize=20)

figname = 'Fig2_'+str(Min)+'-'+str(Max)+'_panels_fc-'+str(fc)+'_NEW.png'
plt.savefig(figname)
print 'Figure is saved as ', figname 

# Note : selecting 17-19 , err < 0.3 ,  

