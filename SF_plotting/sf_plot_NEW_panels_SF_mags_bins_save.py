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
from matplotlib import rcParams

rcParams['ytick.labelsize'] = 25
rcParams['xtick.labelsize'] = 25
rcParams['axes.labelsize'] = 35
rcParams['axes.linewidth'] = 3
rcParams['font.size'] = 25
rcParams.update({'figure.autolayout': False})

from scipy.stats import binned_statistic


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


       
        
#######################
# PLOTTING FUNCTIONS  #
#######################

   
def get_plotted_quantities(data, nbins, pre, save_bin , bin_range,correction ):
    '''
    Create input for sf_plot_panels()
    data : the input data for a given bin, consisting of all del_mag, tau, err per bin
    nbins : number of bins for the entire panel 
    '''    
    
    delflx = data[0]
    tau = data[1]
    delflxerr = data[2] 
          
    # Define functions for bin statistics 
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
    ##### Panel 3 : SF  , Panel 4 : mu_approx   
    #####
    
    if  save_bin == True :
        
        # give correct name depending on whether there is or isnt
        # any correction 
        if correction is not None:
            d = pre+'_bins_'+ bin_range+'_'+correction 
        else: 
            d = pre+'_bins_'+ bin_range
        
        if not os.path.exists(d): #make the dir if it does not exist 
            os.makedirs(d)
        N_bin = [] 
        
        for N in np.unique(bin_number):
            xi = delflx[bin_number == N]
            ei = delflxerr[bin_number == N]
            
            # store the number of points per bin
            N_bin.append(len(xi)) 
        
            if N % 25 == 0:
                print 'Saving bin N=', N
            DATA = np.column_stack((xi,ei))
            outfile = d+'/'+pre+'_bin_'+str(N).zfill(3)+'_xi_ei'+bin_range+'.txt'
            np.savetxt(outfile, DATA, delimiter= ' ', fmt = "%s")            
               
          

    plot_data = {}
    print ' passing on the  plot data...'
    colnames = ['mean_tau', 'bin_stdev', 'err_stdev','err_median', 'N_bin']
    datatable = [mean_tau, bin_stdev, err_stdev,  err_median, N_bin]
                 
    for label, column in zip(colnames, datatable):
        plot_data[label] = column    
    
    return plot_data


def sf_plot_panels(qso_data,star_data_blue, star_data_red, nbins, correction, 
                   save_bin, bin_range):  
    '''
    NEW : instead of sf_plotting, this routine is more versatile, as it 
    calls external function  get_plotted_quantities()  to actually 
    calculate the things to be plotted : means of delmag  per bin, etc. 
    
    It plots the four panels, first getting quantities to plot for stars, 
    and then for quasars, and then plotting them altogether . It also 
    creates a separate figure that plots mu_approx on a linear scale 
    '''               
      
    ##############################################
    # CALCULATE PLOTTED QUANTITIES  : QUASARS 
    ############################################## 
        
    qso_plot  = get_plotted_quantities(qso_data, nbins, 'QSO', save_bin, 
         bin_range,correction )  
    
    ##############################################
    # CALCULATE PLOTTED QUANTITIES  : RED STARS  
    ##############################################
              
    #star_plot = get_plotted_quantities(star_data_blue, nbins,'StarB',save_bin, 
    #                                   bin_range, correction)

    ##############################################
    # CALCULATE PLOTTED QUANTITIES  : BLUE STARS 
    ##############################################
    
   # star_plot1 = get_plotted_quantities(star_data_red, nbins,'StarR',save_bin, 
    #                                    bin_range, correction)
     
    return qso_plot



# inside the main loop : get tau, delflx from a master file, either qso or star
def add_tau_delflx(masterFiles, inDir, good_ids, i, data, fc):
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
    
    if fc is not None :  # correct new master rows only if  asked for 
        err = np.append(err, master[:,2][master_mask].astype(float)*fc)
    else:                # otherwise read in without any correction
        err = np.append(err, master[:,2][master_mask].astype(float))
    master_acc_list  = np.append(master_acc_list, master_acc)
    
    return delflx, tau, err, master_acc_list
    
def read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                 good_ids_QSO, xi_ei_data=None, fc=None):
    inDir_S       = inDirStars
    good_ids_S_blue    = good_ids_S_blue
    good_ids_S_red    = good_ids_S_red
    inDir_Q       = inDirQSO
    good_ids_Q    = good_ids_QSO
    
    
    # Read the Stellar Master file names 
    masterFiles_S = os.listdir(inDir_S)
    
     # Read the QSO Master file names 
    masterFiles_Q = os.listdir(inDir_Q)

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
    
    for i in range(len(masterFiles_Q)): #  len(masterFiles_Q)
        qso_data = add_tau_delflx(masterFiles_Q,inDir_Q, good_ids_Q, i, 
                                  qso_data, fc)
        print np.shape(qso_data)
        star_data_blue = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_blue, i, 
                                   star_data_blue, fc)
        
        star_data_red = add_tau_delflx(masterFiles_S, inDir_S, good_ids_S_red, i, 
                                   star_data_red, fc)                            
                                   
    
    return  qso_data, star_data_blue, star_data_red
    
inDirStars   = 'sf_TRY/sf_stars/'
inDirQSO = 'sf_TRY/sf_qso/'


##############################################################################
##############################################################################
##############################################################################

no_corr= False
w_corr = True

#
#  with correction
#

if w_corr == True:

     
    #
    # Part 1 :  17-19 
    #
    
    mMin = [17,18,18.5]
    mMax = [18,18.5,19]
   #fc = [1.0,1.0,1.0]
    fc = [0.72, 0.91, 1.07]
    #fc=[1.3,1.3,1.3]
    names = ['qso', 'starB', 'starR']
  
   
    out = None
    
    for i in range(len(mMin)): # len(mMin) 
        Min = mMin[i]
        Max = mMax[i]
        print('\nUsing now only lightcurves with SDSS  %f< r_mMed < %f' % (Min, Max))
        print('\n Using fc=%f' % fc[i])
        
        # Select the magnitude range input 
        good_ids_S_blue  = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1, gi_Max=1)
        good_ids_S_red = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, gi_Max=3)
        good_ids_QSO, mask_qso = cut_qso(mMin = Min, mMax=Max, mErrMax = 0.3)
        
        
        # Read in the corresponding xi, ei  lines 
        # from qso, starB, starR
        # From all master files 
        
        
        qso, starB, starR = read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                      good_ids_QSO,xi_ei_data=out, fc=fc[i])
                      
      
        out = [qso, starB, starR]
        
    #  Just save bin data using all master files and all chunks together     
    pl_1 = sf_plot_panels(qso,starB,starR, nbins=200, 
                          correction='corr', save_bin=True, bin_range = '17-19')
    #
    # Part 2 :  18.5-19 
    #
   
    
    names = ['qso', 'starB', 'starR']
    Min = 18.5
    Max = 19
    fc = 1.07 #1.3 # 1.07
    print('Using now only lightcurves with SDSS  %f< r_mMed < %f' % (Min, Max))
    
    # Select the magnitude range input 
    good_ids_S_blue  = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1, gi_Max=1)
    good_ids_S_red = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, gi_Max=3)
    good_ids_QSO, mask_qso = cut_qso(mMin = Min, mMax=Max, mErrMax = 0.3)
    
    # Read in the corresponding xi, ei  lines 
    # from qso, starB, starR
    # From all master files 
    qso, starB, starR = read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO, xi_ei_data=None, fc=fc)
                  
  
    
    # Just save bin data using all master files and all chunks together     
    pl_2 = sf_plot_panels(qso, starB, starR, nbins=200,
                          correction='corr', save_bin=True, bin_range = '18.5-19')


#
# no correction
#

# Initialize Figure 

if no_corr == True:
   
    #
    # Part 1 :  17-19 
    # initialize dict to store all chunks - stitch them together 
    #
    
    Min = 17
    Max = 19
    print('Using now only lightcurves with SDSS  %f< r_mMed < %f' % (Min, Max))
    
    # Select the magnitude range input 
    good_ids_S_blue  = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1, gi_Max=1)
    good_ids_S_red = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, gi_Max=3)
    good_ids_QSO, mask_qso = cut_qso(mMin = Min, mMax=Max, mErrMax = 0.3)
    
    # Read in the corresponding xi, ei  lines from all master files 
    qso, starB, starR = read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO, xi_ei_data=None, fc=None)
                      
    #  Just save bin data using all master files and all chunks together     
    pl_3 = sf_plot_panels(qso, starB, starR, nbins=200, 
                          correction=None, save_bin=True, bin_range = '17-19')
    #
    # Part 2 :  18.5-19 
    #
    
    Min = 18.5
    Max = 19
    
    print('Using now only lightcurves with SDSS  %f< r_mMed < %f' % (Min, Max))
    
    # Select the magnitude range input 
    good_ids_S_blue  = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1, gi_Max=1)
    good_ids_S_red = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, gi_Max=3)
    good_ids_QSO, mask_qso = cut_qso(mMin = Min, mMax=Max, mErrMax = 0.3)
    
    # Read in the corresponding xi, ei  lines 
    # from qso, starB, starR
    # From all master files 
    qso, starB, starR = read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
                  good_ids_QSO, xi_ei_data=None, fc=None)
                 
    #  Just save bin data using all master files and all chunks together     
    pl_4 = sf_plot_panels(qso, starB, starR, nbins=200,
                          correction=None, save_bin=True, bin_range = '18.5-19')
             
  