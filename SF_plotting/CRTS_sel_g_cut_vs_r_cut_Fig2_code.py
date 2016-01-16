# -*- coding: iso-8859-1 -*-
import os
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import binned_statistic

def get_qso_catalog(catalog):
    if catalog == 's82drw':
        File = 'CRTS_SDSS_cross_matched_qso_s82drw_catalog.txt'
    if catalog == 'DB_QSO':
        File = 'CRTS_SDSS_cross_matched_qso_DB_QSO_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File, dtype=str)
    qso_catalog = {}
    print 'Zipping CRTS-SDSS quasars catalog from ', File, ' ...'
    for label, column in zip(colnames, datatable.T):
        qso_catalog[label] = column
    
    print 'Read in ', len(qso_catalog['redshift']), ', quasars from CRTS'
    return  colnames, qso_catalog
    
def get_stars_catalog():
    File = 'CRTS_SDSS_cross_matched_stars_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    stars_catalog = {}
    print 'zipping CRTS-SDSS stars catalog...'
    for label, column in zip(colnames, datatable.T):
        stars_catalog[label] = column
        
    return  colnames, stars_catalog

cols1, qso_cat = get_qso_catalog(catalog='DB_QSO') 
cols2 , star_cat= get_stars_catalog()

# Perform cuts 
def cut_qso(qso_cat=qso_cat, mMin=-9, mMax=19, 
            mErrMin = -9, mErrMax = 0.3,cut_mag='r', report_mag = 'r'):

    mask_mag = (qso_cat[cut_mag].astype(float) > mMin) * (qso_cat[cut_mag].astype(float) < mMax) 
    mask_err = (qso_cat['CRTS_avg_e'].astype(float) > mErrMin) * (qso_cat['CRTS_avg_e'].astype(float) < mErrMax)
    mask = mask_mag * mask_err 
    qso_id = qso_cat['CRTS_id'][mask]
    qso_mags = qso_cat[report_mag][mask]
    print '\n These cuts reduced the number of qso  in the sample from', \
          len(qso_cat['redshift']), ' to ', len(qso_id)
    return  qso_id

def cut_stars(star_cat=star_cat, mMin=-9, mMax=19, mErrMin = -9, 
              mErrMax = 0.3, gi_Min = -1, gi_Max=1 , cut_mag='r_mMed',
              report_mag = 'r_mMed'):

    mask_mag = (star_cat[cut_mag] > mMin) * (star_cat[cut_mag] < mMax) 
    mask_err = (star_cat['CRTS_Merr'] > mErrMin) * (star_cat['CRTS_Merr'] < mErrMax)
    SDSS_gi = star_cat['g_mMed'] - star_cat['i_mMed']
    mask_color = (SDSS_gi > gi_Min ) * (SDSS_gi < gi_Max)
    mask = mask_mag * mask_err * mask_color
    star_id_f = star_cat['crts_id'][mask]
    star_mags = star_cat[report_mag][mask]
    # convert floats to strings without comma and zeros
    star_id = np.array(["{:.0f}".format(name) for name in star_id_f])
    print '\n These cuts reduced the number of stars  in the sample from', \
          len(star_cat['CRTS_M']), ' to ', len(star_id)
    return  star_id


Min = 17
Max = 19
magnitudes = ['g','r']

objects_in_cut = {}

for mag in magnitudes : 
    cut_mag = mag
    report_mag = mag
    
    print('\nUsing now only lightcurves with SDSS  %f< %s < %f' % (Min, cut_mag, Max))
    print('\n Reporting SDSS %s  '% report_mag)

    good_ids_S_blue = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1,
                                              gi_Max=1, cut_mag=cut_mag + '_mMed', 
                                              report_mag=report_mag + '_mMed')
    
    good_ids_S_red = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, 
                                           gi_Max=3, cut_mag=cut_mag + '_mMed', 
                                           report_mag=report_mag + '_mMed')
    
    good_ids_QSO = cut_qso(mMin = Min, mMax=Max, mErrMax = 0.3, 
                                               cut_mag=cut_mag,report_mag=report_mag)
    objects_in_cut[mag] = {'starsB':good_ids_S_blue, 'starsR':good_ids_S_red, 
                           'qso':good_ids_QSO}
    
##################      3


bins = {}
bin_types = ['g_cut','r_cut']

objects = objects_in_cut['g'].keys()

# first need to explicitly initialize the dictionaries 
for b in bin_types:
    bins[b] = {}
    
for obj in objects : 
    bins['g_cut'][obj] = objects_in_cut['g'][obj]
    bins['r_cut'][obj] =  objects_in_cut['r'][obj]

    
##################      4

# inside the main loop : get tau, delflx from a master file, either qso or star
def add_tau_delflx(File, inDir, data, fc):
    # read in storage arrays
    delflx = data[0]  
    tau = data[1]
    err = data[2]
    master_acc_list = data[3]   
    
    # grab the object name 
    master_name = File[3:-4]
    
    # read in the i-th master file 
    master =  np.genfromtxt(inDir+File, dtype=str)
    
    # read in tau,  del_mag,  del_mag_err for quasars on the list 
    delflx = np.append(delflx, master[:,0].astype(float))
    tau = np.append(tau, master[:,1].astype(float))
    
    if fc is not None :  # correct new master rows only if  asked for 
        err = np.append(err, master[:,2].astype(float)*fc)
    else:                # otherwise read in without any correction
        err = np.append(err, master[:,2].astype(float))
    master_names  = np.append(master_acc_list, np.array(len(master[:,0])*[master_name]))
    
    return delflx, tau, err, master_names
    
def read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                 good_ids_QSO, xi_ei_data=None, fc=None):
                     
    inDir_S       = inDirStars
    good_ids_S_blue    = good_ids_S_blue
    good_ids_S_red    = good_ids_S_red
    inDir_Q       = inDirQSO
      
    
    # Read the Stellar Master file names 
    masterFiles_S = os.listdir(inDir_S)
    masterFilesS1 = [name[3:-4] for name in masterFiles_S]
    
    good_masterSB = np.array(masterFiles_S)[np.in1d(masterFilesS1, good_ids_S_blue)]
    good_masterSR = np.array(masterFiles_S)[np.in1d(masterFilesS1, good_ids_S_red)]
    
    # Read the QSO Master file names 
    masterFiles_Q = os.listdir(inDir_Q)
    masterFilesQ1 = [name[3:-4] for name in masterFiles_Q]
    good_masterQ = np.array(masterFiles_Q)[np.in1d(masterFilesQ1, good_ids_QSO)]
    

  
    # If no previous read-in xi, ei exists, initialize arrays    
    if xi_ei_data is None : 
        print 'making new delflx, tau, xi arrays'
        delflx_S      = np.empty(0,dtype=float)
        tau_S         = np.empty(0,dtype=float)
        err_S         = np.empty(0,dtype=float)
        master_acc_list_S = np.empty(0, dtype=str)
    
        
       
        delflx_Q      = np.empty(0,dtype=float)
        tau_Q         = np.empty(0,dtype=float)
        err_Q         = np.empty(0,dtype=float)
        master_acc_list_Q = np.empty(0, dtype=str)
        
        # Initialize the data structures to which more and more delta_t and delta_mag
        # are addded from each consecutive master file 
        qso_data = [delflx_Q, tau_Q, err_Q, master_acc_list_Q] 
        star_data_blue = [delflx_S, tau_S, err_S, master_acc_list_S]
        star_data_red  = [delflx_S, tau_S, err_S, master_acc_list_S]
        
    else:
        print 'using existing xi ei arrays'
        qso_data = xi_ei_data[0]
        star_data_blue = xi_ei_data[1]
        star_data_red = xi_ei_data[2]
        
    print('\n')
    c = 0
    for File in good_masterQ: #  len(masterFiles_Q)
        qso_data = add_tau_delflx(File,inDir_Q, qso_data, fc)
        c += 1 
        if c % 5 == 0:
            pers = (100.0*c) / float(len(good_masterQ))
            print('\r----- Already read %d%% of qso'%pers),
    
    print('\n')
    c = 0                   
    for File in good_masterSB:    # [:len(good_masterQ)]
        star_data_blue = add_tau_delflx(File, inDir_S,star_data_blue, fc)
        c += 1 
        if c % 5 == 0:
            pers = (100.0*c) / float(len(good_masterSB))
            print('\r----- Already read %d%% of Blue Stars'%pers),  
    print('\n')
    c = 0                         
    for File in good_masterSR:   # [:len(good_masterQ)]
        star_data_red = add_tau_delflx(File, inDir_S, star_data_red, fc)      
        c += 1               
        if c % 5 == 0:
            pers = (100.0*c) / float(len(good_masterSR))
            print('\r----- Already read %d%% of Red Stars'%pers),          
                     
    print('returning xi, ei for ... %d objects'%len(good_masterQ))
                            
    return  qso_data, star_data_blue, star_data_red


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

qso, starB, starR = read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
              good_ids_QSO,xi_ei_data=None, fc=None)

# put into a dictionary : makes it more explicit 
out_dic[b] = {'starsB': starB, 'starsR': starR, 'qso':qso}


# Straight after reading-in xi, ei,   one can proceed directly to part 9) (one bin) or 10 : all bins sigma comparison 
# or to Saving just the log(tau) samples of xi, tau, ei. 


# ################   11 aa 

# Run the calculation over all bins : 
# - keep the option to plot a bin if one wants to compare the three methods : need to provide object N
#   
# - as it is, calculate three values for sigma, etc., and plot it  

# Set correction factor
fc = 1.0

b = 'r_cut'   # or g_cut
obj = 'qso'  # or starsB,  starsR 
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
    
    


from astroML.stats import median_sigmaG
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

#####
##### Panel 3 (SF)   and Panel 4   (mu)
#####

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

fname =  b+'_'+str(Min)+'-'+str(Max)+'_'+obj+'_fc-'+str(fc)+'_mean_tau_sig_sigG_SF_mu.txt'
print 'saved as ', fname

data = np.column_stack((mean_tau, bin_stdev, err_stdev, bin_sigma_G, err_sigma_G, 
                        sig_bins['approx'], sig_bins['approx_err'],
                        mu_bins['approx'],mu_bins['approx_err'] ))

header = 'mean_tau   bin_stdev   err_stdev  bin_sigma_G  err_sigma_G   SF  SF_err  mu  mu_err  '
np.savetxt(fname, data, fmt = '%s', delimiter = ' ' , header=header )


################     11 bb

Min=17
Max=19
obj='qso'
fc=1.0

fname =  b+'_'+str(Min)+'-'+str(Max)+'_'+obj+'_fc-'+str(fc)+'_mean_tau_sig_sigG_SF_mu.txt'
print fname

colnames = open(fname,'r').read().splitlines()[0][1:].split()
d = np.genfromtxt(fname, dtype=float)
plot_data = {}

for label, column in zip(colnames, d.T):
    plot_data[label] = column


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

fig,ax = plt.subplots(4,1, figsize=(8,12), sharex=True)
fig.subplots_adjust(hspace=0)

# Panel 1 

ax[0].scatter(np.log10(plot_data['mean_tau']), plot_data['bin_stdev'], s=p_size, 
                alpha=p_al, c = col1)
ax[0].errorbar(np.log10(plot_data['mean_tau']), plot_data['bin_stdev'],plot_data['err_stdev'], 
               linestyle='None', c = col1  )

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

ax[1].scatter(np.log10(plot_data['mean_tau']), plot_data['bin_sigma_G'], s=p_size, 
                alpha=p_al, c = col1)
ax[1].errorbar(np.log10(plot_data['mean_tau']), plot_data['bin_sigma_G'],plot_data['err_sigma_G'], 
               linestyle='None', c = col1  )


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

ax[2].scatter(np.log10(plot_data['mean_tau']), plot_data['SF'], s=p_size, 
                alpha=p_al, c = col1)
ax[2].errorbar(np.log10(plot_data['mean_tau']), plot_data['SF'],plot_data['SF_err'], 
               linestyle='None', c = col1  )

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

ax[3].scatter(np.log10(plot_data['mean_tau']), plot_data['mu'], s=p_size, 
                alpha=p_al, c = col1)
ax[3].errorbar(np.log10(plot_data['mean_tau']), plot_data['mu'],plot_data['mu_err'], 
               linestyle='None', c = col1  )


ax[3].axhline(y=0.0, color='black', lw=lh_w, ls=lh_st,alpha=lh_al)
ax[3].set_ylim(top=y_mu_top, bottom=y_mu_bott)
ax[3].set_xlim(left=x_left, right=x_right)
ax[3].set_yticks([-0.05,0,0.05])
ax[3].set_yticklabels(['-0.05','0.0', '0.05'])  
ax[3].set_ylabel(r'$\mu$', fontsize=20)
ax[3].grid(axis='x')
ax[3].set_xlabel(r'$log_{10} (\Delta _{t})$ [days]',fontsize=20)

figname = 'Fig2_'+str(Min)+'-'+str(Max)+'_'+obj+'_fc-'+str(fc)+'.png'
plt.savefig(figname)

# Note : selecting 17-19 , err < 0.3 ,  


