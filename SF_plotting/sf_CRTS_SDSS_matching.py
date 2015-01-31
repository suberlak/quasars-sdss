# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 16:10:39 2015

@author: suberlak

Match CRTS stars/ quasars to SDSS stars / quasars

Big files with lots of data : 

SDSS stars : ~/Desktop/SDSS/stars_SDSS/stripe82calibStars_v2.6.dat  (many cols)
SDSS quasars : ~/Desktop/SDSS/s82drw/s82drw_g.dat (Chelsea's fits have redshift....)

Small files with just  names of each object for which we have a lightcurve:

CRTS stars: ~/Desktop/SDSS/SF_plotting/stars_bright_names.txt   stars_faint_names.txt
CRTS quasars : ~/Desktop/SDSS/SF_plotting/quasars_bright_names.txt   quasars_faint_names.txt


"""

import numpy as np


matching_type = 'star'  # star or qso   ? 

############################
####  SDSS STARS   #########
############################

if matching_type == 'star' : 
    File = '../stars_SDSS/stripe82calibStars_v2.6.dat' 
    datatable = np.genfromtxt(File)
    colnames = ['calib_fla', 'ra', 'dec', 'raRMS', 'decRMS', 'nEpochs', 'AR_val', 
                'u_Nobs', 'u_mMed', 'u_mMean', 'u_mErr', 'u_rms_scatt', 'u_chi2',
                'g_Nobs', 'g_mMed', 'g_mMean', 'g_mErr', 'g_rms_scatt', 'g_chi2',
                'r_Nobs', 'r_mMed', 'r_mMean', 'r_mErr', 'r_rms_scatt', 'r_chi2',
                'i_Nobs', 'i_mMed', 'i_mMean', 'i_mErr', 'i_rms_scatt', 'i_chi2',
                'z_Nobs', 'z_mMed', 'z_mMean', 'z_mErr', 'z_rms_scatt', 'z_chi2']
    
    data_stars = {}
    print 'Zipping the stars...'
    
    for label, column in zip(colnames, datatable.T):
        data_stars[label] = column
        
    print 'I read in data for ', len(data_stars['ra']), ' SDSS stars'

############################
#### SDSS QUASARS  #########
############################

if matching_type == 'qso' : 
    File = '../s82drw/s82drw_g.dat'
    datatable = np.genfromtxt(File)
    colnames = ['SDR5ID', 'ra', 'dec', 'redshift', 'M_i', 'mass_BH', 
                'chi^2_pdf', 'log10(tau[days])', 'log10(sigma_hat)',
                'log10(tau_lim_lo)', 'log10(tau_lim_hi)',  'log10(sig_lim_lo)',
                'log10(sig_lim_hi)', 'edge_flag', 'Plike',  'Pnoise', 'Pinf',
                'mu', 'npts' ]
    
    data_quasars = {}
    print 'zipping quasars...'
    for label, column in zip(colnames, datatable.T):
        data_quasars[label] = column
        
    print 'I read in data for ', len(data_quasars['ra']), ' SDSS quasars'
############################
####  CRTS STARS #######\
############################
# name scheme  out_1000647.dat.txt

if matching_type == 'stars':
    
    # load names of stars 
    File = 'stars_bright_names.txt'
    CRTS_stars = np.loadtxt(File, dtype= str)
    File = 'stars_faint_names.txt'
    CRTS_stars2 = np.loadtxt(File, dtype= str)
    CRTS_stars = np.append(CRTS_stars, CRTS_stars2)
    
    print 'I loaded names of ', len(CRTS_stars), ' CRTS stars'
    
    # load ra dec info for matching...
    File = 'radec.00'
    CRTS_strs_radec_table= np.genfromtxt(File)
    colnames = ['CRTS_ID', 'ra', 'dec']
    
    data_CRTS_stars = {}
    for label, column in zip(colnames, CRTS_strs_radec_table.T):
        data_CRTS_stars[label] = column

# merge the two, and at the  end save as separate files...

############################
#### CRTS QUASARS ######
############################
#  name scheme : out_000456.17+000645.5.txt

if matching_type == 'qso':
    File = 'qso_bright_names.txt'
    CRTS_qso = np.loadtxt(File, dtype= str)
    
    File = 'qso_faint_names.txt'
    CRTS_qso1 = np.loadtxt(File, dtype= str)
    
    CRTS_qso = np.append(CRTS_qso, CRTS_qso1)
    
    
    print 'I loaded names of ', len(CRTS_qso), ' CRTS quasars'
    # merge the two, and at the end save as separate files...




