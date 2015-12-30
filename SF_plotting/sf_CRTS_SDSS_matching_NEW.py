# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 16:52:41 2015

@author: suberlak

New matching program, using astropy  instead of my custom-writen functions
--> I have done it just like that guy (writing my own angular -separation routines, etc.)
    https://www.sites.google.com/site/mrpaulhancock/blog/theage-oldproblemofcross-matchingastronomicalsources

--> But now I can do it  using astropy built-in modules, which should really save all 
    computational time! 

UPDATE: 2015/12/16
--> changed the was in which CRTS QSO is matched to DB_QSO, adding another column
     to the saved CRTS_SDSS_cross_matched_qso_DB_QSO_catalog.txt   : crts_id, 
    which is the 000453+123443  coordinate that also serves as a  name for CRTS 
    quasars  .  Now the DB_QSO saving format is %s,  which is easier to code 

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

def load_sdss_qso(catalog = 'DB_QSO'):
    '''
    Choose where the catalog is coming : is it 
    DB_QSO  : DB_QSO_S82.dat, with SDSS ugriz, for 9258 objects
    master_QSO : master database with SDSS-2MASS cross-matched photometry &
                 errors for 8974 objects
    s82_drw : s82_drw*.dat  files, which are very similar to DB_QSO, except 
              that they quote DRW fitting results for a given SDSS band 
              photometry (hence 5 s82_drw*  files where '*' = u,g,r,i,z )
              
    NOTE:  all catalogs taken from
    http://www.astro.washington.edu/users/ivezic/cmacleod/qso_dr7/Southern.html
    
    '''
    if catalog == 's82drw' : 
        # note : ra and dec in degrees
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
        print 'From catalog ', File
        return data_quasars
        
    if catalog == 'DB_QSO' : 
        # note : ra and dec in degrees
        File = 'DB_QSO_S82.dat'
        datatable = np.genfromtxt(File)
        colnames = ['dbID', 'ra', 'dec', 'SDR5ID', 'M_i', 'M_i_corr', 'redshift',
                   'mass_BH', 'Lbol', 'u', 'g', 'r', 'i', 'z', 'Au']
        data_quasars = {}
        print 'zipping quasars...'
        for label, column in zip(colnames, datatable.T):
            data_quasars[label] = column
        
        print 'I read in data for ', len(data_quasars['ra']), ' SDSS quasars'
        print 'From catalog ', File
        return data_quasars
        
    if catalog == 'master_QSO' :
        File = 'master_QSO_S82.dat'
        datatable = np.genfromtxt(File, usecols=np.arange(14))
        colnames = ['dbID', 'ra', 'dec','redshift', 'u', 'err_u', 'g', 'err_g',
                    'r', 'err_r', 'i', 'err_i', 'z','err_z'] 
        data_quasars = {}
        print 'zipping quasars...'
        for label, column in zip(colnames, datatable.T):
            data_quasars[label] = column
        
        print 'I read in data for ', len(data_quasars['ra']), ' SDSS quasars'
        print 'From catalog ', File
        return data_quasars            
   
    
############################
####  CRTS STARS #######\
############################
# name scheme  out_1000647.dat.txt

def load_crts_stars(inDir):
    # load names of stars 
    inFiles = os.listdir(inDir)
    CRTS_star_names = np.loadtxt(inFiles, dtype= str)
    
    print '\nI loaded names of ', len(CRTS_star_names), ' CRTS stars'
    
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
    # load names of quasars, which already contain ra and dec infor 
    CRTS_qso  = os.listdir(inDir)
    print '\nI loaded names of ', len(CRTS_qso), ' CRTS quasars'
    return CRTS_qso

crts_dirs = ['../QSO_CRTS_processed_err_w/','../stars_CRTS_proc_err_w_good/']


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
    

def match_catalogs(cat1_ra, cat1_dec, cat2_ra, cat2_dec, archive_file):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    cat1 = SkyCoord(ra=cat1_ra*u.degree, dec=cat1_dec*u.degree)
    cat2 = SkyCoord(ra=cat2_ra*u.degree, dec=cat2_dec*u.degree)
    idx, sep2d, dist3d = cat1.match_to_catalog_sky(cat2)
    #np.savez(archive_file, SDSS_matching_rows = idx, matched_radius=sep2d ) 
    return idx, sep2d 
    
##############  ACTION  : MATCHING STARS  ############### 
    
    
def match_stars(purge_pickle=True):
    # load names from CRTS
    DIR = crts_dirs[1]
    crts_star_names, crts_star_radec = load_crts_stars(DIR)
     
    ########################################################## 
    # EXTRACT LC PARAMETERS FOR EACH STAR OBSERVED WITH CRTS #
    ##########################################################
     
    archive_file='CRTS_stars_LC_params.npz'
    # Check whether this has not been done already :
    if not os.path.exists(archive_file) or  purge_pickle:
        length= len(crts_star_names)
        print '\n- Computing average mag, err , extracting ra, dec for %i points' % length
        
        avg_mag=[]
        avg_err=[]
        ra_ls  =[]
        dec_ls =[]
        crts_id=[]
        mjd_span=[]
        mjd_uniq_N=[]
        N_rows= []
        c=0
        for i in range(length):
            file = str(crts_star_names[i])
            #print '\nCRTS stars file ',i, 'out of',  length
            mjd,flx4,err = np.loadtxt(DIR+'%s' % (file),usecols=(0,1,2),unpack=True)
            # 1) Average brightness per LC
            avg_mag.append(np.mean(flx4))
            
            # 2) Average error per LC
            avg_err.append(np.mean(err))
            
            # 3) Delta_MJD : span in days between final and first day
            mjd_span.append(int(mjd.max()-mjd.min()))  
            
            # 4) N_MJD : number of days (i.e. taking the integer part of mjd 
            # date and picking only those that are unique)
            unique_mjds = np.unique([int(day) for day in mjd])
            mjd_uniq_N.append(len(unique_mjds))  # number of days observed 
            
            # 5) N_rows : number of rows per LC 
            N_rows.append(len(mjd))
            
            # 6) CRTS ID  
            crts_id_i= float(crts_star_names[i][4:-8])
            crts_id.append(crts_id_i)
            
            # 7) CRTS  ra and dec for that object 
            name_mask = crts_star_radec['CRTS_ID'] == crts_id_i
            ra_ls.append(crts_star_radec['ra'][name_mask][0])
            dec_ls.append(crts_star_radec['dec'][name_mask][0])#
            
            c += 1 
            if c % 5 == 0:
                pers = (100.0*c) / float(length)
                print('\r----- Already read %d%% of stars '%pers), 

        print '\nSaving the results of all LC parameters for CRTS stars to...'
        print archive_file
        np.savez(archive_file, avg_mag=avg_mag, avg_err=avg_err, ra_ls=ra_ls, 
                 dec_ls=dec_ls, crts_id=crts_id, mjd_span=mjd_span, 
                 mjd_uniq_N=mjd_uniq_N, N_rows=N_rows )   
                 
    else: 
        print '\n- Using precomputed LC parameters (avg_mag, err, ra,dec, crts_id, ',\
        ' mjd_span, mjd_uniq_N, N_rows) for CRTS stars from ...'
        print archive_file
        archive = np.load(archive_file)
        avg_mag    = archive['avg_mag']
        avg_err    = archive['avg_err']
        ra_ls      = archive['ra_ls']
        dec_ls     = archive['dec_ls']
        crts_id    = archive['crts_id']
        mjd_span   = archive['mjd_span']
        mjd_uniq_N = archive['mjd_uniq_N']
        N_rows     = archive['N_rows']
    # My CRTS coordinates are already in  degrees...    
    ra_deg_CRTS = ra_ls
    dec_deg_CRTS = dec_ls  


    ########################################################## 
    # MATCH CRTS TO SDSS (ROW BY ROW)                        #
    ##########################################################
  
    archive_file_matching = 'CRTS_SDSS_stars_matched_rows_radii.npz'
    if not os.path.exists(archive_file_matching) or purge_pickle :
        print '\n- Computing the SDSS matching rows to CRTS stars'
          #     Load data from SDSS
        sdss_star_data =  load_sdss_stars()
        SDSS_matching_rows , matched_radius= match_catalogs(cat1_ra=ra_deg_CRTS, 
                                                            cat1_dec=dec_deg_CRTS, 
                                                            cat2_ra= sdss_star_data['ra'], 
                                                            cat2_dec=sdss_star_data['dec'], 
                                                            archive_file=archive_file_matching) 
        np.savez(archive_file_matching, SDSS_matching_rows=SDSS_matching_rows)
    else:
        print '\n- Using precomputed SDSS rows matched to CRTS stars from'
        print archive_file_matching
        archive =np.load(archive_file_matching)
        SDSS_matching_rows= archive['SDSS_matching_rows']
        sdss_star_data =  load_sdss_stars()
    
    ########## SAVE ################
    # Saving a combined cross-matched SDSS-CRTS stars dataset 
    
    ind = SDSS_matching_rows
    datatable=np.array([avg_mag, avg_err,  sdss_star_data['dec'][ind], sdss_star_data['ra'][ind],
                         dec_deg_CRTS, ra_deg_CRTS, sdss_star_data['g_Nobs'][ind], 
                        sdss_star_data['g_mMed'][ind], sdss_star_data['r_mMed'][ind],
                        sdss_star_data['i_mMed'][ind],
                        crts_id, mjd_span, mjd_uniq_N, N_rows ])
    colnames = ['CRTS_M','CRTS_Merr', 'dec_SDSS', 'ra_SDSS', 'dec_CRTS',
                'ra_CRTS', 'g_Nobs', 'g_mMed','r_mMed', 'i_mMed', 'crts_id', 'mjd_span', 
                'mjd_N', 'N_rows']
    
#    colnames = ['calib_fla', 'ra', 'dec', 'raRMS', 'decRMS', 'nEpochs', 'AR_val', 
#                'u_Nobs', 'u_mMed', 'u_mMean', 'u_mErr', 'u_rms_scatt', 'u_chi2',
#                'g_Nobs', 'g_mMed', 'g_mMean', 'g_mErr', 'g_rms_scatt', 'g_chi2',
#                'r_Nobs', 'r_mMed', 'r_mMean', 'r_mErr', 'r_rms_scatt', 'r_chi2',
#                'i_Nobs', 'i_mMed', 'i_mMean', 'i_mErr', 'i_rms_scatt', 'i_chi2',
#                'z_Nobs', 'z_mMed', 'z_mMean', 'z_mErr', 'z_rms_scatt', 'z_chi2']
#    
    
    data_SDSS_CRTS= {}
    print 'Zipping the stars...'
    
    for label, column in zip(colnames, datatable):
        data_SDSS_CRTS[label] = column
    print 'I made a dictionary with data for ', len(data_SDSS_CRTS['dec_SDSS']), ' SDSS-CRTS cross-matched stars'
    
    print 'Saving the SDSS-CRTS cross-matched stars catalog...' 
    
    archive_SDSS_CRTS = 'CRTS_SDSS_cross_matched_stars_catalog.txt' 
    keys = colnames
    DATA = np.column_stack((datatable))    
    
    header=''
    for key in keys: 
        header= header+'{:<10}'.format(key[:10])+' '
    
    fmt = ['%s', '%.4e', '%10.5f']   # formatters to choose from...  
    
    np.savetxt(archive_SDSS_CRTS, DATA, delimiter =' ', fmt=fmt[2], header=header)    
    print 'All done with star catalogs, please see: ' , archive_SDSS_CRTS
    
    return SDSS_matching_rows, ra_deg_CRTS, dec_deg_CRTS, avg_mag, avg_err, sdss_star_data 
       
    #  fmt='%11.5f'*8+'%6.i'+'%5.i'*2
    
    
    
    
##############  ACTION  : MATCHING QUASARS  ############### 

    
def match_quasars(catalog, purge_pickle=True):
    '''
    For catalog matching need to choose SDSS catalog: 
    
    
    'DB_QSO'  : DB_QSO_S82.dat, with SDSS ugriz, for 9258 objects
    'master_QSO' : master database with SDSS-2MASS cross-matched photometry &
                 errors for 8974 objects
    's82_drw' : s82_drw*.dat  files, which are very similar to DB_QSO, except 
              that they quote DRW fitting results for a given SDSS band 
              photometry (hence 5 s82_drw*  files where '*' = u,g,r,i,z )
              
    '''
    # load names from CRTS 
    DIR = crts_dirs[0]
    crts_qso_names = load_crts_qso(DIR)
    
    # load data from SDSS
    sdss_qso_data =  load_sdss_qso(catalog=catalog)
        
    # LOOP OVER CRTS QUASARS 
    archive_file='CRTS_qso_LC_params.npz'
    # Check whether this has not been done already :
    if not os.path.exists(archive_file) or purge_pickle :
        length= len(crts_qso_names)
        print '- Computing average mag, err , extracting ra, dec for %i points' % length
        
        avg_mag=[]
        avg_err=[]
        ra_ls =[]
        dec_ls=[]
        mjd_span = []
        mjd_uniq_N = []
        N_rows = []        
        
        c=0
        for i in range(length):
            file = str(crts_qso_names[i])
            #print '\nCRTS quasars file ',i, 'out of',  length
            mjd,flx4,err = np.loadtxt(DIR+'%s' % (file),usecols=(0,1,2),unpack=True)
            
            # 1) Average brightness per LC
            avg_mag.append(np.mean(flx4))
            
            # 2) Average error per LC
            avg_err.append(np.mean(err))
            
            # 3) Delta_MJD : span in days between final and first day
            mjd_span.append(int(mjd.max()-mjd.min()))  
            
            # 4) N_MJD : number of days (i.e. taking the integer part of mjd 
            # date and picking only those that are unique)
            unique_mjds = np.unique([int(day) for day in mjd])
            mjd_uniq_N.append(len(unique_mjds))  # number of days observed 
            
            # 5) N_rows : number of rows per LC 
            N_rows.append(len(mjd))
            
            # 6) CRTS  ra and dec for that object ( no  need to pull ra, dec 
            # from a separate file, matching by name, because the qso name
            # already includes that... )
            ra_ls.append(file[4:13])
            dec_ls.append(file[13:-4])
            
            c += 1 
            if c % 5 == 0:
                pers = (100.0*c) / float(length)
                print('\r----- Already read %d%% of QSO '%pers), 

        print '\nSaving the results of all LC parameters for CRTS quasars to...'
        print archive_file   
        np.savez(archive_file, avg_mag=avg_mag, avg_err=avg_err, ra_ls=ra_ls, 
                 dec_ls=dec_ls, mjd_span=mjd_span, mjd_uniq_N=mjd_uniq_N, 
                 N_rows=N_rows )   
                 
    else: 
        print '\n - Using precomputed LC parameters (avg_mag, err, ra,dec, crts_id, ',\
        ' mjd_span, mjd_uniq_N, N_rows) for CRTS quasars  from ...'
    
        archive = np.load(archive_file)
        avg_mag    = archive['avg_mag']
        avg_err    = archive['avg_err']
        ra_ls      = archive['ra_ls']
        dec_ls     = archive['dec_ls']
        mjd_span   = archive['mjd_span']
        mjd_uniq_N = archive['mjd_uniq_N']
        N_rows     = archive['N_rows']
       
    # Split CRTS  ra, dec from hms to h m s 
    ra_hms_split, dec_hms_split = get_ra_dec_CRTS(ra_ls, dec_ls)
    # Convert CRTS  ra, dec from hms to deg  
    ra_deg_CRTS, dec_deg_CRTS = convert_to_deg(ra_hms_split, dec_hms_split)
    
    # Matching CRTS to SDSS  : which SDSS row corresponds to which CRTS row... 
    archive_file_matching = 'CRTS_SDSS_qso_'+catalog+'_matched_rows.npz'
    
    if not os.path.exists(archive_file_matching) or purge_pickle :
        print '\n- Computing the SDSS matching rows to CRTS quasars'
        SDSS_matching_rows , matched_radius= match_catalogs(cat1_ra=ra_deg_CRTS, 
                                                            cat1_dec=dec_deg_CRTS, 
                                                            cat2_ra= sdss_qso_data['ra'], 
                                                            cat2_dec=sdss_qso_data['dec'], 
                                                            archive_file=archive_file_matching) 
        np.savez(archive_file_matching, SDSS_matching_rows=SDSS_matching_rows)
        print '\n- Saved the SDSS-CRTS quasars matched rows to ', archive_file_matching
   
    else:
        print '\n- Using precomputed SDSS rows matched to CRTS quasars from ', archive_file_matching
        archive =np.load(archive_file_matching)
        SDSS_matching_rows = archive['SDSS_matching_rows']
        
        
    # Saving a combined cross-matched SDSS-CRTS quasars dataset 
        
    ind = SDSS_matching_rows
    
    sdss_qso_data.keys()
    
    
    # Save the list of names of CRTS Quasars 
    # independent of SDSS catalog choice, because 
    # saving here only CRTS names to which SDSS is matched 
    
    qso_names_file = 'CRTS_SDSS_cross_matched_qso_names.txt'
    np.savetxt(qso_names_file, crts_qso_names, fmt='%s')
    print '\nSaving the CRTS quasar file names to',   qso_names_file
    
    # Make a list of just CRTS QSO id for catalog identification
    crts_qso_names_radec = [name[4:-4] for name in crts_qso_names]    
    
    # Save all other measurable quantities for CRTS - SDSS quasars 
    if catalog == 's82drw' : 
        datatable=np.array([avg_mag, avg_err, sdss_qso_data['M_i'][ind], 
                          sdss_qso_data['redshift'][ind], sdss_qso_data['ra'][ind], 
                          sdss_qso_data['dec'][ind], ra_deg_CRTS, dec_deg_CRTS,
                          mjd_span, mjd_uniq_N, N_rows])
        colnames = ['CRTS_avg_mag','CRTS_avg_err','M_i', 'redshift', 'dec_CRTS',
                    'ra_CRTS',  'dec_SDSS','ra_SDSS', 'mjd_span', 'mjd_uniq_N',
                    'N_rows']
        # NOTE: colnames is read from the right....
    if catalog == 'DB_QSO' :
        datatable=np.array([crts_qso_names_radec, avg_mag, avg_err, sdss_qso_data['z'][ind], 
                          sdss_qso_data['i'][ind], sdss_qso_data['r'][ind],
                          sdss_qso_data['g'][ind], sdss_qso_data['u'][ind], 
                          sdss_qso_data['redshift'][ind], sdss_qso_data['ra'][ind], 
                          sdss_qso_data['dec'][ind], ra_deg_CRTS, dec_deg_CRTS,
                          mjd_span, mjd_uniq_N, N_rows])
        colnames = ['CRTS_id','CRTS_avg_mag','CRTS_avg_err', 'z', 'i', 'r', 'g', 'u', 
                    'redshift', 'ra_SDSS', 'dec_SDSS', 'ra_CRTS', 'dec_CRTS', 
                    'mjd_span', 'mjd_uniq_N', 'N_rows']
    if catalog == 'master_QSO' :
        datatable=np.array([avg_mag, avg_err, sdss_qso_data['u'][ind], 
                          sdss_qso_data['err_u'][ind], sdss_qso_data['g'][ind], 
                          sdss_qso_data['err_g'][ind], sdss_qso_data['r'][ind],
                          sdss_qso_data['err_r'][ind], sdss_qso_data['i'][ind],
                          sdss_qso_data['err_i'][ind], sdss_qso_data['z'][ind], 
                          sdss_qso_data['err_z'][ind], sdss_qso_data['redshift'][ind], 
                          sdss_qso_data['ra'][ind], sdss_qso_data['dec'][ind], 
                          ra_deg_CRTS, dec_deg_CRTS, mjd_span, mjd_uniq_N, N_rows])
        colnames = ['CRTS_avg_mag','CRTS_avg_err', 'u', 'err_u', 'g', 'err_g',
                    'r', 'err_r', 'i', 'err_i', 'z','err_z', 
                    'redshift', 'ra_SDSS', 'dec_SDSS', 'ra_CRTS', 'dec_CRTS', 
                    'mjd_span', 'mjd_uniq_N', 'N_rows']

    data_qso_SDSS_CRTS= {}
    print '\nZipping  quasars...'
    
    for label, column in zip(colnames, datatable):
        data_qso_SDSS_CRTS[label] = column
    print 'I made a dictionary with data for ', len(data_qso_SDSS_CRTS['redshift']), \
          ' CRTS-SDSS cross-matched quasars'
    
    archive_SDSS_CRTS_qso = 'CRTS_SDSS_cross_matched_qso_'+catalog+'_catalog.txt'

    print '\nSaving the SDSS-CRTS cross-matched QSO catalog...' 
    print ' to ', archive_SDSS_CRTS_qso
     
    keys = colnames
    if catalog == 's82drw':
        DATA = np.column_stack((datatable))   
                                #  data_qso_SDSS_CRTS[keys[10]], data_qso_SDSS_CRTS[keys[9]], 
                               # data_qso_SDSS_CRTS[keys[8]],  data_qso_SDSS_CRTS[keys[7]], 
                               # data_qso_SDSS_CRTS[keys[6]],  data_qso_SDSS_CRTS[keys[5]],
                               # data_qso_SDSS_CRTS[keys[4]],  data_qso_SDSS_CRTS[keys[3]], 
                               # data_qso_SDSS_CRTS[keys[2]],  data_qso_SDSS_CRTS[keys[1]],
                               # data_qso_SDSS_CRTS[keys[0]]
        
        header=''
        # http://stackoverflow.com/questions/766141/reverse-a-string-in-python 
        for key in keys : 
            header= header+'{:<10}'.format(key[:10])+' '
        
        #fmt = ['%s', '%.4e', '%10.5f']
        np.savetxt(archive_SDSS_CRTS_qso, DATA, delimiter =' ', fmt='%11.5f'*8+'%6.i'+'%5.i'*2, header=header)
    
    if catalog == 'DB_QSO':
        DATA = np.column_stack((datatable))    
        
        header=''
        for key in keys : 
            header= header+'{:<10}'.format(key[:10])+' '

        #fmt = ['%s', '%.4e', '%10.5f']
        # old fmt '%11.5f'*12+'%6.i'+'%5.i'*2
        np.savetxt(archive_SDSS_CRTS_qso, DATA, delimiter =' ', fmt='%s  '*16, header=header)
  
    
    return data_qso_SDSS_CRTS

# Call all the necessary functions
#SDSS_idx, ra_deg_CRTS, dec_deg_CRTS, avg_mag, avg_err, sdss_star_data  = match_stars() 
crts = match_quasars(catalog='DB_QSO', purge_pickle=True)  # or 'master_qso',  's82drw'
#match_stars()