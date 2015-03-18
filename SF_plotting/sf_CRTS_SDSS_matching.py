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
import os

############################
####  SDSS STARS   #########
############################

def load_sdss_stars():

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
    
    return data_stars

############################
#### SDSS QUASARS  #########
############################

def load_sdss_qso():
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
    
    return data_quasars
    
############################
####  CRTS STARS #######\
############################
# name scheme  out_1000647.dat.txt

def load_crts_stars(inDir):
    # load names of stars 
    inFiles = os.listdir(inDir)
    CRTS_star_names = np.loadtxt(inFiles, dtype= str)
    
    print 'I loaded names of ', len(CRTS_star_names), ' CRTS stars'
    
    # load ra dec info for matching...
    File = 'radec.00'
    CRTS_star_radec_table= np.genfromtxt(File)
    colnames = ['CRTS_ID', 'ra', 'dec']
    
    CRTS_star_radec = {}
    for label, column in zip(colnames, CRTS_star_radec_table.T):
        CRTS_star_radec[label] = column
    
    return CRTS_star_names, CRTS_star_radec

# merge the two, and at the  end save as separate files...

############################
#### CRTS QUASARS ######
############################
#  name scheme : out_000456.17+000645.5.txt

def load_crts_qso(inDir):
    # load names of quasars
    inFiles = 'out.list'
    CRTS_qso = np.loadtxt(inDir+inFiles, dtype= str)  
    print 'I loaded names of ', len(CRTS_qso), ' CRTS quasars'
    # merge the two, and at the end save as separate files...
    return CRTS_qso

crts_dirs = ['../QSO_CRTS_processed_err_w/','../stars_CRTS_proc_err_w_good/']


############################
###  MATCHING FUNCTIONS  ###
############################

# from  ../s82drw/chelsea_results_load_crts.py

##############################
# CONVERT RA DEC TO DEGREES  #
##############################

def get_ra_dec_CRTS(ra_hms, dec_hms):
    """
    Extracting RA, DEC information from the QSO  name for the CRTS case 
    """
    
    dec_hms_split = np.empty(0, dtype=str)
    for i in range(0,len(dec_hms)):
        dec_hms_split = np.append(dec_hms_split,  dec_hms[i][0:3]+' '+dec_hms[i][3:5]+' '+dec_hms[i][5:9])
    
    
    ra_hms_split = np.empty(0, dtype=str)
    for i in range(0,len(ra_hms)):
        ra_hms_split = np.append(ra_hms_split, ra_hms[i][0:2]+' '+ ra_hms[i][2:4] + ' ' + ra_hms[i][4:9])
    return ra_hms_split, dec_hms_split
     
def HMS2deg(ra='', dec=''):
    """
    From http://www.bdnyc.org/2012/10/15/decimal-deg-to-hms/  
    Converting  ra and dec from h:m:s   and deg:m:s  to  degrees.decimal 
    I assume they are using ICRS coordinates 
    """
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
      
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)
      
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC
    
    
def convert_to_deg(ra_split, dec_split):
    '''
    Converts ra and dec from h:m:s  extracted from the quasar name with 
    get_ra_dec_CRTS()  to  degrees.decimal , using HMS2deg() function
    '''
    ra_deg  =   np.empty(0, dtype=float)
    dec_deg =   np.empty(0, dtype=float)
    
    for i in range(0,len(dec_split)):
        dec_deg = np.append(dec_deg,float(HMS2deg(dec=dec_split[i])))
    
    
    for i in range(0,len(ra_split)):
        ra_deg = np.append(ra_deg,float(HMS2deg(ra=ra_split[i])))
    return ra_deg , dec_deg
    



################################ 
# MATCH CRTS AND SDSS DATASETS # 
################################


DEG_PER_HR = 360. / 24.             # degrees per hour
DEG_PER_MIN = DEG_PER_HR / 60.      # degrees per min
DEG_PER_S = DEG_PER_MIN / 60.       # degrees per sec
DEG_PER_AMIN = 1./60.               # degrees per arcmin
DEG_PER_ASEC = DEG_PER_AMIN / 60.   # degrees per arcsec
RAD_PER_DEG = np.pi / 180.             # radians per degree
def ang_sep(ra1, dec1, ra2, dec2):
    """
    John Ruan's  matching code 
    """
    ra1 = np.asarray(ra1);  ra2 = np.asarray(ra2)
    dec1 = np.asarray(dec1);  dec2 = np.asarray(dec2)

    ra1 = ra1 * RAD_PER_DEG           # convert to radians
    ra2 = ra2 * RAD_PER_DEG
    dec1 = dec1 * RAD_PER_DEG
    dec2 = dec2 * RAD_PER_DEG

    sra1 = np.sin(ra1);  sra2 = np.sin(ra2)
    cra1 = np.cos(ra1);  cra2 = np.cos(ra2)
    sdec1 = np.sin(dec1);  sdec2 = np.sin(dec2)
    cdec1 = np.cos(dec1);  cdec2 = np.cos(dec2)

    csep = cdec1*cdec2*(cra1*cra2 + sra1*sra2) + sdec1*sdec2

    # An ugly work-around for floating point issues.
    #if np.any(csep > 1):  print csep
    csep = np.where(csep > 1., 1., csep)

    degsep = np.arccos(csep) / RAD_PER_DEG
    # only works for separations > 0.1 of an arcsec or  >~2.7e-5 dec
    degsep = np.where(degsep < 1e-5, 0, degsep)
    return degsep
    
    
#########################
# Matching  function    #
#########################

def match_CRTS_to_SDSS(ra_deg_CRTS, dec_deg_CRTS, sdss_data, archive_file):
    
    dec_deg = dec_deg_CRTS
    ra_deg = ra_deg_CRTS
    SDSS_matching_rows = np.zeros_like(ra_deg, dtype=int)
    mask_mismatched = np.zeros_like(SDSS_matching_rows,dtype=bool)
    matched_radius = np.zeros_like(ra_deg,dtype=float)
    
    
    ra_sdss_deg = sdss_data['ra']
    dec_sdss_deg = sdss_data['dec']
    not_yet_matched = np.ones_like(ra_sdss_deg, dtype=bool) # starts all True, 
     # and is set to false when an sdss object is definitely matched to CRTS one.
    
    # Looping over CRTS data (shorter than SDSS)
    for i in range(0,len(ra_deg)):
        angular_separation = ang_sep(ra_deg[i],dec_deg[i], ra_sdss_deg[not_yet_matched], 
                                    dec_sdss_deg[not_yet_matched])
        indices = np.where(angular_separation <= 0.001)
        N_sdss_rows = len(np.where(not_yet_matched == True)[0])

        print 'CRTS ra_deg row ', i, 'using  ', N_sdss_rows, 'SDSS rows'
        check = np.array(indices)
        dims = check.shape
        if (dims[1] == 1.0):
            matched_radius[i] = 0.001
            print 'Obj from CRTS coords: ', ra_deg[i], dec_deg[i]
            print 'Obj from SDSS coords: ', ra_sdss_deg[indices], dec_sdss_deg[indices]
            SDSS_matching_rows[i] = int(indices[0])
            not_yet_matched[indices] = False
        else :  
            #chelsea_matching_rows[i] = 999999    # 1e6-1
            mask_mismatched[i] = True 
            print 'For CRTS object ra_deg',ra_deg[i] ,'dec_deg',\
            dec_deg[i],'There is / are',dims[1],' matching entries from SDSS results'
       
         
    # secondary matching : allow bigger margin for those that were not matched : case by case
    print '\n Secondary matching'
    
    for i in range(0,len(ra_deg[mask_mismatched])):
        ttl_index = np.where(mask_mismatched == True)[0][i]
        N_sdss_rows = len(np.where(not_yet_matched == True)[0])
        print '> Matching ra_jav_deg index number', ttl_index, 'using  ', N_sdss_rows, 'SDSS rows'
        angular_separation = ang_sep(ra_deg[mask_mismatched][i],dec_deg[mask_mismatched][i], 
                                    ra_sdss_deg[not_yet_matched], dec_sdss_deg[not_yet_matched])
        indices = np.where(angular_separation <= 0.001)
        check = np.array(indices)
        dims = check.shape
        radius = 0.001
        if(dims[1] == 0.0 ): 
            print 'No matched object for QSO at javelin coords ra:', \
            ra_deg[mask_mismatched][i], ' dec: ',  dec_deg[mask_mismatched][i]  
            while True:
                radius = radius + 0.001    
                # print 'Increase matching distance to ',radius
                indices = np.where(angular_separation <= radius)
                check = np.array(indices)
                dims = check.shape
                if(dims[1] ==1):
                    print 'Matched to ' , ra_sdss_deg[indices], dec_sdss_deg[indices]
                    print 'Matched radius is', radius 
                    matched_radius[ttl_index] = radius
                    SDSS_matching_rows[ttl_index] = int(indices[0])
                    break
        if(dims[1] > 1.0):
            print 'Too many objects matched for QSO at javelin coords ra:',\
            ra_deg[mask_mismatched][i], ' dec: ',  dec_deg[mask_mismatched][i]
            while True:
                radius = radius - 0.0001    
                print 'Decrease matching distance to ',radius
                indices = np.where(angular_separation <= radius)
                check = np.array(indices)
                dims = check.shape
                if(dims[1] ==1.0):
                    print 'Matched to ' , ra_sdss_deg[indices], dec_sdss_deg[indices]
                    matched_radius[ttl_index] = radius
                    SDSS_matching_rows[ttl_index] = int(indices[0])
                    break
    
    
    
    # checking whether they were all increased 
    mk=np.array(np.where(mask_mismatched == True))
    dims = mk.shape
    
    if( len(np.where(mask_mismatched == True)[0]) == len(np.where(matched_radius > 0.001)[0])) :
        print '\n All mismatched objects had the matching radius increased to make a match'
        for i in range(0,dims[1]):
            if(np.where(matched_radius > 0.001)[0][i] == np.where(mask_mismatched == True)[0][i]):
                print 'The radius was increased from 0.001 to ', matched_radius[mask_mismatched][i]

    np.savez(archive_file, SDSS_matching_rows = SDSS_matching_rows, matched_radius=matched_radius ) 
    return SDSS_matching_rows, matched_radius 


##############  ACTION  : MATCHING QUASARS  ###############


def match_quasars():

    # load names from CRTS \
    DIR = crts_dirs[0]
    crts_qso_names = load_crts_qso(DIR)
    
    # load data from SDSS
    sdss_qso_data =  load_sdss_qso()
        
    # LOOP OVER QUASARS 
    archive_file_qso='CRTS_qso_avg_mag_err_ra_dec.npz'
    # Check whether this has not been done already :
    if not os.path.exists(archive_file_qso) :
        length= len(crts_qso_names)
        print '- Computing average mag, err , extracting ra, dec for %i points' % length
        
        avg_mag=[]
        avg_err=[]
        ra_ls =[]
        dec_ls=[]
        
        for i in range(length):
            file = str(crts_qso_names[i])
            print '\nFile ',i, 'out of',  length
            mjd,flx4,err = np.loadtxt(DIR+'%s' % (file),usecols=(0,1,2),unpack=True)
            avg_mag.append(np.mean(flx4))
            avg_err.append(np.mean(err))
            ra_ls.append(file[4:13])
            dec_ls.append(file[13:-4])
        np.savez(archive_file_qso, avg_mag=avg_mag, avg_err=avg_err, ra_ls=ra_ls, 
                 dec_ls=dec_ls )   
                 
    else: 
        print '- Using precomputed CRTS qso average values for mag, err, and ra, dec results'
        archive = np.load(archive_file_qso)
        avg_mag = archive['avg_mag']
        avg_err  = archive['avg_err']
        ra_ls = archive['ra_ls']
        dec_ls  = archive['dec_ls']
       
    # Split CRTS  ra, dec from hms to h m s 
    ra_hms_split, dec_hms_split = get_ra_dec_CRTS(ra_ls, dec_ls)
    # Convert CRTS  ra, dec from hms to deg  
    ra_deg_CRTS, dec_deg_CRTS = convert_to_deg(ra_hms_split, dec_hms_split)
    
    # Matching CRTS to SDSS  : which SDSS row corresponds to which CRTS row... 
    archive_file_matching = 'CRTS_SDSS_qso_matched_rows_radii.npz'
    
    if not os.path.exists(archive_file_matching) :
        print '- Computing the SDSS matching rows to CRTS quasars'
        SDSS_matching_rows , matched_radius= match_CRTS_to_SDSS(ra_deg_CRTS, dec_deg_CRTS, sdss_data = sdss_qso_data, archive_file=archive_file_matching) 
    else:
        print '- Using precomputed SDSS rows matched to CRTS quasars'
        archive =np.load(archive_file_matching)
        SDSS_matching_rows = archive['SDSS_matching_rows']
        matched_radius = archive['matched_radius']
        
        
        
    # Saving a combined cross-matched SDSS-CRTS quasars dataset 
        
    ind = SDSS_matching_rows
    
    sdss_qso_data.keys()
    datatable=np.array([avg_mag, avg_err, sdss_qso_data['M_i'][ind], sdss_qso_data['redshift'][ind], 
               sdss_qso_data['ra'][ind], sdss_qso_data['dec'][ind], ra_deg_CRTS, dec_deg_CRTS, matched_radius])
    colnames = ['CRTS_avg_mag','CRTS_avg_err','M_i', 'redshift', 'dec_CRTS', 'ra_CRTS', 'dec_SDSS','ra_SDSS', 'R_match']
    # colnames is read from the right....
    
    data_qso_SDSS_CRTS= {}
    print 'Zipping the stars...'
    
    for label, column in zip(colnames, datatable):
        data_qso_SDSS_CRTS[label] = column
    print 'I made a dictionary with data for ', len(data_qso_SDSS_CRTS['R_match']), ' SDSS-CRTS cross-matched quasars'
    
    print 'Saving the SDSS-CRTS cross-matched QSO catalog...' 
    
    archive_SDSS_CRTS_qso = 'CRTS_SDSS_cross_matched_qso_catalog.txt' 
    keys = colnames
    DATA = np.column_stack((data_qso_SDSS_CRTS[keys[8]], 
                            data_qso_SDSS_CRTS[keys[7]], 
                            data_qso_SDSS_CRTS[keys[6]],data_qso_SDSS_CRTS[keys[5]],
                            data_qso_SDSS_CRTS[keys[4]], data_qso_SDSS_CRTS[keys[3]], 
                            data_qso_SDSS_CRTS[keys[2]], data_qso_SDSS_CRTS[keys[1]],
                            data_qso_SDSS_CRTS[keys[0]]))    
    
    header=''
    for key in keys[::-1] : 
        header= header+'{:<10}'.format(key[:10])
    
    
    
    fmt = ['%s', '%.4e', '%10.5f']
    
    np.savetxt(archive_SDSS_CRTS_qso, DATA, delimiter =' ', fmt=fmt[2], header=header)

#match_quasars()       


##############  ACTION  : MATCHING STARS  ###############

 

def match_stars():
    # load names from CRTS
    DIR = crts_dirs[1]
    
    crts_star_names, crts_star_radec = load_crts_stars(DIR)
        
    # LOOP OVER STARS 
    archive_file_qso='CRTS_stars_avg_mag_err_ra_dec.npz'
    # Check whether this has not been done already :
    if not os.path.exists(archive_file_qso) :
        length= len(crts_star_names)
        print '- Computing average mag, err , extracting ra, dec for %i points' % length
        
        avg_mag=[]
        avg_err=[]
        ra_ls  =[]
        dec_ls =[]
    
        
        for i in range(length):
            file = str(crts_star_names[i])
            print '\nFile ',i, 'out of',  length
            mjd,flx4,err = np.loadtxt(DIR+'%s' % (file),usecols=(0,1,2),unpack=True)
            avg_mag.append(np.mean(flx4))
            avg_err.append(np.mean(err))
            
            name_mask = crts_star_radec['CRTS_ID'] == float(crts_star_names[i][4:-8])
            ra_ls.append(crts_star_radec['ra'][name_mask][0])
            dec_ls.append(crts_star_radec['dec'][name_mask][0])
        np.savez(archive_file_qso, avg_mag=avg_mag, avg_err=avg_err, ra_ls=ra_ls, 
                 dec_ls=dec_ls )   
                 
    else: 
        print '- Using precomputed CRTS qso average values for mag and  err results'
        archive = np.load(archive_file_qso)
        avg_mag = archive['avg_mag']
        avg_err  = archive['avg_err']
        ra_ls = archive['ra_ls']
        dec_ls  = archive['dec_ls']
        
    # My CRTS coordinates are already in  degrees...
        
    ra_deg_CRTS = ra_ls
    dec_deg_CRTS = dec_ls  
    # Matching CRTS to SDSS  : which SDSS row corresponds to which CRTS row... 
    archive_file_matching = 'CRTS_SDSS_stars_matched_rows_radii.npz'
    
    return ra_ls, dec_ls 
    
    # Load data from SDSS
#    sdss_star_data =  load_sdss_stars()
#    
#    if not os.path.exists(archive_file_matching) :
#        print '- Computing the SDSS matching rows to CRTS stars'
#        SDSS_matching_rows , matched_radius= match_CRTS_to_SDSS(ra_deg_CRTS, dec_deg_CRTS, sdss_data = sdss_star_data, archive_file=archive_file_matching) 
#    else:
#        print '- Using precomputed SDSS rows matched to CRTS stars'
#        archive =np.load(archive_file_matching)
#        SDSS_matching_rows= archive['SDSS_matching_rows']
#        matched_radius = archive['matched_radius']
#        
#        
        
   

ra, dec = match_stars()