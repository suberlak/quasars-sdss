# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 13:07:52 2014

@author: suberlak


To match Chelsea's results of fitting SDSS quasars in stripe 82,  s82drw_r.dat
to the results of my fits of CRTS quasars in stripe 82  in   
../QSO_CRTS_analysis/javelin_CRTS_chain_results_err_rms.txt  (or ~err_w.txt)

1) Read from s82drw  file the  ra, dec,  log(sigma_hat) and log(tau) into the 
separate vectors. Recalculate sigma_hat into  sigma 

" short-term driving amplitude (called sigma or sigma-hat) can be related to 
the asymptotic structure function by SF_inf = sigma*sqrt(tau/365)" (MacLeod 2014, priv comm.)

Chelsea's program returns  log10(tau) and  log10(sigma_hat), so this has to be 
recalculated to get tau and sigma_hat 

2) Read in results of Javelin fits: filename, tau and sigma low, med, high values
(sampling of the posterior chain at different points of the 1D histogram) which 
come from javelin_chain_retrieve_CRTS.py  as 
np.column_stack((files_read, sigma_l, sigma_m, sigma_h, tau_l, tau_m, tau_h))

Values for both sigma and tau are the real ones, i.e. javelin  produces  from 
get_hpd  ln(tau) , and I took  np.exp()  of those values to get tau  

"""

import numpy as np

#################################
# load JAVELIN results for CRTS #
#################################


pre = '../'
dir_in = pre+'QSO_CRTS_analysis/'
CRTS_output= ['javelin_CRTS_chain_results_err_rms.txt', 'javelin_CRTS_chain_results_err_w.txt']

err_choice = 1 #  (or 1)

out_names = ['javelin_CRTS_err_rms_Chelsea_s82drw_r_compare.txt',  'javelin_CRTS_err_w_Chelsea_s82drw_r_compare.txt']
output = out_names[err_choice]

javelin_results = np.genfromtxt(dir_in + CRTS_output[err_choice], dtype="str")

print 'working file', CRTS_output[err_choice]
quasar_names_ra_dec = javelin_results[:,0]
sigma_l =  javelin_results[:,1].astype(np.float)
sigma_m =  javelin_results[:,2].astype(np.float)
sigma_h =  javelin_results[:,3].astype(np.float)
tau_l =  javelin_results[:,4].astype(np.float)
tau_m =  javelin_results[:,5].astype(np.float)
tau_h =  javelin_results[:,6].astype(np.float)

ra_jav_deg  =   np.empty(0, dtype=float)
dec_jav_deg =   np.empty(0, dtype=float)

#################################
# load CHELSEA results for SDSS #
#################################

lines = tuple(open('s82drw_r.dat', 'r'))

# the problem is that the last column sometimes is not separated by whitespace 
# from the preceding column. Thus I read the entire file as string
# and then, row by row, separate by whitespaces, and only take the elements that I need 
# (I am not using the 17/18 columns with mnu, #pts , anyways  )
# split 26-th row into separate entries, which are columns
# lines[26].split()
# 1-st element : ra , of the 26-th row 
# lines[26].split()[1]
# initialise data storing vectors 

ra_ch  = np.empty(0, dtype=float)
dec_ch = np.empty(0, dtype=float)
log_tau_med   = np.empty(0, dtype=float)
log_sigma_hat_med = np.empty(0, dtype=float)

# fill in vectors in values 

for i in range(3,len(lines)):   # from 3 because of the header
    ra_ch = np.append(ra_ch,float(lines[i].split()[1]))
    dec_ch= np.append(dec_ch,float(lines[i].split()[2]))
    log_tau_med   = np.append(log_tau_med,float(lines[i].split()[7])) 
    log_sigma_hat_med = np.append(log_sigma_hat_med,float(lines[i].split()[8])) 

# calculate sigma from sigma_hat for Chelsea results 
tau_ch = np.power(10,log_tau_med)
sigma_hat = np.power(10,log_sigma_hat_med)
sigma_ch = sigma_hat * np.sqrt(tau_ch / (2.0*365.0))

##############################
# CONVERT RA DEC TO DEGREES  #
##############################

def get_ra_dec_CRTS():
    """
    Extracting RA, DEC information from the QSO  name for the CRTS case 
    """
    ra_jav = np.empty(0, dtype=float)
    dec_jav= np.empty(0, dtype=float)
    
    for i in range(0,len(quasar_names_ra_dec)):
        ra_jav  = np.append(ra_jav,quasar_names_ra_dec[i][0:9])
        dec_jav = np.append(dec_jav, quasar_names_ra_dec[i][9:18])
    
    dec_jav_split = np.empty(0, dtype=str)
    for i in range(0,len(dec_jav)):
        dec_jav_split = np.append(dec_jav_split,  dec_jav[i][0:3]+' '+dec_jav[i][3:5]+' '+dec_jav[i][5:9])
    
    
    ra_jav_split = np.empty(0, dtype=str)
    for i in range(0,len(ra_jav)):
        ra_jav_split = np.append(ra_jav_split, ra_jav[i][0:2]+' '+ ra_jav[i][2:4] + ' ' + ra_jav[i][4:9])
    return ra_jav_split, dec_jav_split
     
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
    
ra_jav_split, dec_jav_split = get_ra_dec_CRTS()
ra_jav_deg, dec_jav_deg = convert_to_deg(ra_jav_split, dec_jav_split)

        
    
def diagnostic_print(qso_names, ra_split, dec_split, ra_deg, dec_deg, index):
    print qso_names[index]
    print ra_split[index]
    print dec_split[index]
    print ra_deg[index]
    print dec_deg[index]


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
    
    

chelsea_matching_rows = np.zeros_like(ra_jav_deg, dtype=int)
mask_mismatched = np.zeros_like(chelsea_matching_rows,dtype=bool)
matched_radius = np.zeros_like(ra_jav_deg,dtype=float)
# mask=np.zeros_like(chelsea_matching_rows,dtype=bool)

for i in range(0,len(ra_jav_deg)):
    matched_distances = ang_sep(ra_jav_deg[i],dec_jav_deg[i], ra_ch, dec_ch)
    indices = np.where(matched_distances <= 0.001)
    print 'ra_jav_deg row ', i 
    check = np.array(indices)
    dims = check.shape
    if (dims[1] == 1.0):
        matched_radius[i] = 0.001
        print 'Obj from javelin coords: ', ra_jav_deg[i], dec_jav_deg[i]
        print 'Obj from Chelsea coords: ', ra_ch[indices], dec_ch[indices]
        chelsea_matching_rows[i] = int(indices[0])
    else :  
        #chelsea_matching_rows[i] = 999999    # 1e6-1
        mask_mismatched[i] = True 
        print 'For object ra_jav_deg',ra_jav_deg[i] ,'dec_jav_deg',\
        dec_jav_deg[i],'There is are',dims[1],' matching entries from Chelsea results'
   
     
# secondary matching : allow bigger margin for those that were not matched : case by case
print '\n Secondary matching'

for i in range(0,len(ra_jav_deg[mask_mismatched])):
    ttl_index = np.where(mask_mismatched == True)[0][i]
    print '> Matching ra_jav_deg index number', ttl_index
    matched_distances = ang_sep(ra_jav_deg[mask_mismatched][i],dec_jav_deg[mask_mismatched][i], ra_ch, dec_ch)
    indices = np.where(matched_distances <= 0.001)
    check = np.array(indices)
    dims = check.shape
    radius = 0.001
    if(dims[1] == 0.0 ): 
        print 'No matched object for QSO at javelin coords ra:', \
        ra_jav_deg[mask_mismatched][i], ' dec: ',  dec_jav_deg[mask_mismatched][i]  
        while True:
            radius = radius + 0.001    
            # print 'Increase matching distance to ',radius
            indices = np.where(matched_distances <= radius)
            check = np.array(indices)
            dims = check.shape
            if(dims[1] ==1):
                print 'Matched to ' , ra_ch[indices], dec_ch[indices]
                print 'Matched radius is', radius 
                matched_radius[ttl_index] = radius
                chelsea_matching_rows[ttl_index] = int(indices[0])
                break
    if(dims[1] > 1.0):
        print 'Too many objects matched for QSO at javelin coords ra:',\
        ra_jav_deg[mask_mismatched][i], ' dec: ',  dec_jav_deg[mask_mismatched][i]
        while True:
            radius = radius - 0.0001    
            print 'Decrease matching distance to ',radius
            indices = np.where(matched_distances <= radius)
            check = np.array(indices)
            dims = check.shape
            if(dims[1] ==1.0):
                print 'Matched to ' , ra_ch[indices], dec_ch[indices]
                matched_radius[ttl_index] = radius
                chelsea_matching_rows[ttl_index] = int(indices[0])
                break

# checking whether they were all increased 
mk=np.array(np.where(mask_mismatched == True))
dims = mk.shape

if( len(np.where(mask_mismatched == True)[0]) == len(np.where(matched_radius > 0.001)[0])) :
    print '\n All mismatched objects had the matching radius increased to make a match'
    for i in range(0,dims[1]):
        if(np.where(matched_radius > 0.001)[0][i] == np.where(mask_mismatched == True)[0][i]):
            print 'The radius was increased from 0.001 to ', matched_radius[mask_mismatched][i]

            
#print ra_ch[chelsea_matching_rows], dec_ch[chelsea_matching_rows]
#print log_tau_med[chelsea_matching_rows], log_sigma_med[chelsea_matching_rows]

####################
# SAVE THE RESULTS #
####################


template = "{0:10}{1:10}{2:8}{3:10}{4:10}{5:12}{6:12}{7:12}{8:12}"
print template.format('ra_jav','dec_jav','ra_ch','dec_ch', 'tau_jav', 'tau_ch', 'sigma_jav', 'sig_ch','sig_ratio' )

col0 = [None]*len(ra_jav_deg)
col1 = np.zeros_like(ra_jav_deg)
col2 = np.zeros_like(ra_jav_deg)
col3= np.zeros_like(ra_jav_deg)
col4= np.zeros_like(ra_jav_deg)
col5= np.zeros_like(ra_jav_deg)
col6= np.zeros_like(ra_jav_deg)
col7= np.zeros_like(ra_jav_deg)
col8= np.zeros_like(ra_jav_deg)
col9= np.zeros_like(ra_jav_deg)

#for i in range(0,len(ra_jav_deg)):
#    if mask_mismatched[i] == False :  # avoids going through those which have more than one match
#        qso = quasar_names_ra_dec[i]
#        col0[i] = qso[0:18]        
#        col1[i] = ra_jav_deg[i]
#        col2[i] = ra_ch[chelsea_matching_rows][i]
#        col3[i] = dec_jav_deg[i]
#        col4[i] = dec_ch[chelsea_matching_rows][i]
#        col5[i] = tau_m[i]
#        col6[i] = np.power(10,log_tau_med[chelsea_matching_rows][i])
#        col7[i] = sigma_m[i]
#        col8[i] = np.power(10,log_sigma_med[chelsea_matching_rows][i])
#        col9[i] = col7[i] / col8[i]
#   
   
for i in range(0,len(ra_jav_deg)):
    qso = quasar_names_ra_dec[i]
    col0[i] = qso[0:18]        
    col1[i] = ra_jav_deg[i]
    col2[i] = ra_ch[chelsea_matching_rows][i]
    col3[i] = dec_jav_deg[i]
    col4[i] = dec_ch[chelsea_matching_rows][i]
    col5[i] = tau_m[i]
    col6[i] = tau_ch[chelsea_matching_rows][i]
    col7[i] = sigma_m[i]
    col8[i] = sigma_ch[chelsea_matching_rows][i]
    
    

DAT= np.column_stack((col0,col1,col2,col3,col4,col5,col6,col7,col8,col9))

print '\n We are done matching javelin and Chelseas results. We matched ', \
len(ra_jav_deg), ' Quasars, with ', len(np.where(mask_mismatched == True)[0]),\
' those that required wider matching radius than 0.001 radians. Results are saved to',\
output," with columns 'ra_jav','dec_jav','ra_ch','dec_ch', 'tau_jav', 'tau_ch', 'sigma_jav', 'sig_ch','sig_ratio' "

# sort the DAT column according to tau_javelin 
# newDAT=DAT[DAT[:,4].argsort()]

np.savetxt(output,DAT,fmt="%s")

# diagnostics ... 
diagnostic_print(quasar_names_ra_dec, ra_jav_split, dec_jav_split, ra_jav_deg, dec_jav_deg, 0  )
