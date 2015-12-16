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

# rename the star key so that I could use the same name throughout  
#star_cat['r'] = star_cat.pop('r_mMed')

# Perform cuts 
def cut_qso(qso_cat=qso_cat, qso_names=qso_names, mMin=-9, mMax=19, 
            mErrMin = -9, mErrMax = 0.3):

    mask_mag = (qso_cat['r'] > mMin) * (qso_cat['r'] < mMax) 
    mask_err = (qso_cat['CRTS_avg_e'] > mErrMin) * (qso_cat['CRTS_avg_e'] < mErrMax)
    mask = mask_mag * mask_err 
    qso_id = qso_names[mask]
    qso_rmags = qso_cat['r'][mask]
    print '\n These cuts reduced the number of qso  in the sample from', \
          len(qso_cat['redshift']), ' to ', len(qso_id)
    return  qso_id, mask , qso_rmags

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
    return  star_id, star_mags


#
#  Since in this program I'm not doing any plotting
#  I simply removed those lines, which are unchanged in the 
#  version  sf_plot_NEW_panels_SF_mags_split.py


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
        #print 'Reading in ', File
        
        qso_data = add_tau_delflx(File,inDir_Q, qso_data, fc)
        c += 1 
        if c % 5 == 0:
            pers = (100.0*c) / float(len(good_masterQ))
            print('\r----- Already read %d%% of qso'%pers),
    
    print('\n')
    c = 0                   
    for File in good_masterSB:    # [:len(good_masterQ)]
        #print 'Reading in ', File
        star_data_blue = add_tau_delflx(File, inDir_S,star_data_blue, fc)
        c += 1 
        if c % 5 == 0:
            pers = (100.0*c) / float(len(good_masterQ))
            print('\r----- Already read %d%% of Blue Stars'%pers),  
    print('\n')
    c = 0                         
    for File in good_masterSR:   # [:len(good_masterQ)]
        #print 'Reading in ', File
        star_data_red = add_tau_delflx(File, inDir_S, star_data_red, fc)      
        c += 1               
        if c % 5 == 0:
            pers = (100.0*c) / float(len(good_masterQ))
            print('\r----- Already read %d%% of Red Stars'%pers),          
                     
    print('returning xi, ei for ... %d objects'%len(good_masterQ))
                            
    return  qso_data, star_data_blue, star_data_red
    
inDirStars   = 'sf_file_per_LC/star/'
inDirQSO = 'sf_file_per_LC/qso/'


##############################################################################
##############################################################################
##############################################################################



# Select the magnitude range input 
cut_mags = False
if cut_mags == True:
    Min = 15
    Max = 16
    cut_mag = 'g_mMed'
    report_mag = 'g_mMed'
    
    print('\nUsing now only lightcurves with SDSS  %f< %s < %f' % (Min, cut_mag, Max))
    print('\n Reporting SDSS %s  '% report_mag)

    good_ids_S_blue, S_blue_mags  = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = -1, gi_Max=1,
                                               cut_mag=cut_mag, report_mag=report_mag)
    good_ids_S_red, S_red_mags = cut_stars(mMin = Min, mMax=Max, mErrMax = 0.3, gi_Min = 1, gi_Max=3,
                                            cut_mag=cut_mag, report_mag=report_mag)
    good_ids_QSO, mask_qso, qso_rmags = cut_qso(mMin = Min, mMax=Max, mErrMax = 0.3)

# read in objects exactly the same as Sample_4
sample_4_redo = True
if sample_4_redo == True:
    
    S_blue_str_zeros = np.loadtxt('Sample_4_stars_b_unique_names.txt', dtype=str)
    good_ids_S_blue = np.array([s.rstrip(".0") for s in S_blue_str_zeros])
    S_blue_float =  S_blue_str_zeros.astype(float)
    S_blue_mags     = np.array([star_cat['g_mMed'][star_cat['crts_id'] == star][0]  for star in S_blue_float])
    
    S_red_str_zeros  = np.loadtxt('Sample_4_stars_r_unique_names.txt',dtype=str)
    good_ids_S_red = np.array([s.rstrip(".0") for s in S_red_str_zeros])
    S_red_float =  S_red_str_zeros.astype(float)
    S_red_mags      = np.array([star_cat['g_mMed'][star_cat['crts_id'] == star][0]  for star in S_red_float])
    
    good_ids_QSO    = np.loadtxt('Sample_4_qso_unique_names.txt', dtype=str)
    qso_rmags       = np.array([qso_cat['r'][qso_names == qso][0] for qso in good_ids_QSO])

obj_type = ['qso', 'starB', 'starR' ]

# Read in the master files  if necessary 




    
read_from_master = True
read_from_prev_xi = False

#if not os.path.exists(re) : 
if read_from_master == True : 
    qso, starB, starR = read_xi_ei(inDirStars, good_ids_S_blue, good_ids_S_red, inDirQSO,
              good_ids_QSO,xi_ei_data=None, fc=None)
          
    # put into a list to loop over... 
    out = [qso, starB, starR]
   
# Save the output of reading-in the master files...

    for i in range(3):   
        re = 'All_xi_tau_ei_name_Sample_4_redo_'+obj_type[0]+'.txt'
        print('\n We are saving delflx, tau, err, master_names')
        print('\n We save read-in master files according to this cut as %s' % re )
        DATA = np.column_stack((out[i][0],out[i][1],out[i][2],out[i][3] ))
        np.savetxt(re, DATA, fmt = '%s', delimiter = ' ')

if read_from_prev_xi == True : 
    
    print('\n This exact cut was already used before to read-in master files.')
    out = [[],[],[]]
    for i in range(3):   
        re = 'All_xi_tau_ei_name_Sample_4_redo_'+obj_type[i]+'.txt'
        print('\n Reading in xi, tau, ei, name from %s' % re)
        data = np.genfromtxt(re, dtype=str)
        out[i] = [data[:,0].astype(float), data[:,1].astype(float), data[:,2].astype(float),
                 data[:,3]] 

            

ids = [good_ids_QSO, good_ids_S_blue, good_ids_S_red]
rep_mags = [qso_rmags, S_blue_mags, S_red_mags]
   # NOTE : rmags and ids are ordered in the same fashion, and have 
   # identical lengths 


# Select only those that satisfy cuts...

tauMin = [0, 2.3, 3]
tauMax= [1.7, 2.5, 3.2]
hist_min = [-1.5, -1.0, -1.0 ]
hist_max = [1.5, 1.0, 1.0]

for i in range(len(out)): # loop over qso, StarB, StarR 
    print '\nFor ', obj_type[i]
    for j in range(1): # len(tauMin) :  loop over selection bins for tau 
        tau_min = tauMin[j]
        tau_max = tauMax[j]
        print('\nWe are selecting a sample for %.1f <log10(tau)< %.1f'%(tau_min, tau_max))    
        
        # take vectors from the given object type
        xi  = out[i][0]
        tau = out[i][1]
        ei  = out[i][2]
        n   = out[i][3]
        
        print 'The median error before any selection is ', np.median(ei)
        print('Before any cut we have  %d objects'% len(np.unique(n))) 
        # make selection
        xmin = hist_min[i]
        xmax = hist_max[i]
        
        mask_delflx = (xi > xmin) * (xi < xmax)  
        
        mask_tau =  (np.log10(tau) > tau_min) * (np.log10(tau)<tau_max)

        print('After log(tau) cut we have %d objects'% len(np.unique(n[mask_tau])))        
        
        mask = mask_tau * mask_delflx
        
        print 'The median error after tau an outliers selection is ', np.median(ei[mask])        
        print('After log(tau) cut and cutting outliers we have %d objects'% len(np.unique(n[mask])))
        # grab names to get magnitudes 
        names_sample_uniq = np.unique(n[mask])
        
        # grab r mags corresponding to names 
        mask_mags = np.in1d(ids[i], names_sample_uniq)
        mags_sample = rep_mags[i][mask_mags]
        
        keys = names_sample_uniq   #  unique names in the sample
        values = mags_sample     # and corresponding rmags 
        
        d = dict(zip(keys,values))
        
        # make an array of rmags corresponding to the original names 
        # using the dictionary : we know what mag corresponds to 
        # which object name, just need to make a new array, 
        # based on the full list of names 
        
        # map mags onto names 
        mags = [d[key] for key in n[mask]]
        
        # save the sample ... 
        
        data = np.column_stack((xi[mask], ei[mask], tau[mask], n[mask], mags))
        dir = '/astro/store/scratch/tmp/suberlak/SF_Samples/Sample_4_redo_try/'
        fname = 'Sample_'+obj_type[i]+'_tau_'+str(tau_min)+'-'+str(tau_max)+'_'+ \
            str(len(mags))+'_lines_'+str(len(names_sample_uniq))+'_objects.txt'
        print 'In directory', dir
        if not os.path.exists(dir): #make the dir if it does not exist 
            os.makedirs(dir)
        print 'Saving the Sample as ', fname    
        np.savetxt(dir+fname, data, delimiter=' ', fmt = '%s')
        
   

 