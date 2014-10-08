# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 13:07:52 2014

@author: suberlak

To load Chelsea McLeod's results from s83drw_... dat file : one file for each
filter  ugriz  .  

s82drw_g.dat  is the file with Chelsea's  results, and we read the ra, dec, 
log(sigma) and log(tau) into the separate vectors. 

Then we read in javelin results, which are stored as fname , and boundaries for 
sigma as well as tau  


javelin outputs:
from the program javelin_chain_retrieve_sigma_tau_save.py  , which stores 

np.column_stack((files_read, sigma_l, sigma_m, sigma_h, tau_l, tau_m, tau_h))
(1) File name   (2)  Lower bound on hpd of sigma 0.16  (3) Middle value of hpd of sigma 0.50
(3) Upper bound on sigma , 0.84  (4),(5), (6) : same for tau . 

Values for both sigma and tau are the real ones, i.e. javelin  produces  from 
get_hpd  ln(tau) , and I took  np.exp()  of those values to get tau  


Chelsea's program stores:
log(tau)

Thus here I need to make  np.power(10,log_tau_med) , etc.  

"""
#from astropy import units as u
#from astropy.coordinates import SkyCoord

import numpy as np
#results = np.genfromtxt('s82drw_g.dat')



lines = tuple(open('s82drw_g.dat', 'r'))
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
log_sigma_med = np.empty(0, dtype=float)

# fill in vectors in values 

for i in range(3,len(lines)):   # (runs to len(lines)-1 )
    ra_ch = np.append(ra_ch,float(lines[i].split()[1]))
    dec_ch= np.append(dec_ch,float(lines[i].split()[2]))
    log_tau_med   = np.append(log_tau_med,float(lines[i].split()[7])) 
    log_sigma_med = np.append(log_sigma_med,float(lines[i].split()[8])) 


# comparison with the javelin results 
# javelin results come  from  
# 

javelin_results = np.genfromtxt('../QSO_try/fname_hpd-sigma_hpd-tau.txt', dtype="str")

quasar_names_ra_dec = javelin_results[:,0]
sigma_l =  javelin_results[:,1].astype(np.float)
sigma_m =  javelin_results[:,2].astype(np.float)
sigma_h =  javelin_results[:,3].astype(np.float)
tau_l =  javelin_results[:,4].astype(np.float)
tau_m =  javelin_results[:,5].astype(np.float)
tau_h =  javelin_results[:,6].astype(np.float)

# now need to extract  ra  and  dec from the name string
# first need to assess the format : anticipated 
# 123456.89+123456.8_chain.dat
# thus  name[0:9] is  the ra
# and  name[10:18] is the dec (including the +/-  sign)

ra_jav = np.empty(0, dtype=float)
dec_jav= np.empty(0, dtype=float)

for i in range(0,len(quasar_names_ra_dec)):
    ra_jav  = np.append(ra_jav,quasar_names_ra_dec[i][0:9])
    dec_jav = np.append(dec_jav, quasar_names_ra_dec[i][9:18])
    

# then need to convert ra and dec from h:m:s   and deg:m:s
# to  degrees.decimal 

# I assume they are using ICRS coordinates  

dec_jav_split = np.empty(0, dtype=str)
for i in range(0,len(dec_jav)):
    dec_jav_split = np.append(dec_jav_split,  dec_jav[i][0:3]+' '+dec_jav[i][3:5]+' '+dec_jav[i][5:9])


ra_jav_split = np.empty(0, dtype=str)
for i in range(0,len(ra_jav)):
    ra_jav_split = np.append(ra_jav_split, ra_jav[i][0:2]+' '+ ra_jav[i][2:4] + ' ' + ra_jav[i][4:9])
    
    

"""
Code below from http://www.bdnyc.org/2012/10/15/decimal-deg-to-hms/  
because astropy .coordinates does not seem to work! 
It FAILS  

astropy.test()

'4 failed, 5036 passed, 87 skipped, 12 xfailed in 171.09 seconds'

Ask someone to update astropy distribution!!

In the meantime use the code below  - does what is needed  
"""
def HMS2deg(ra='', dec=''):
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
    
    
ra_jav_deg  =   np.empty(0, dtype=float)
dec_jav_deg =   np.empty(0, dtype=float)

for i in range(0,len(dec_jav_split)):
    dec_jav_deg = np.append(dec_jav_deg,float(HMS2deg(dec=dec_jav_split[i])))


for i in range(0,len(ra_jav_split)):
    ra_jav_deg = np.append(ra_jav_deg,float(HMS2deg(ra=ra_jav_split[i])))
    
"""
Compare stages:

print quasar_names_ra_dec[0]
print ra_jav[0]
print dec_jav[0]
print ra_jav_split[0]
print dec_jav_split[0]
print ra_jav_deg[0]
print dec_jav_deg[0]
"""

# Perform matching : calculate the cartesian distance between the two positions, 
# and if the distance is smaller than a given value, assume we are speaking of the same object  

"""
John Ruan's  matching code 
"""

DEG_PER_HR = 360. / 24.             # degrees per hour
DEG_PER_MIN = DEG_PER_HR / 60.      # degrees per min
DEG_PER_S = DEG_PER_MIN / 60.       # degrees per sec
DEG_PER_AMIN = 1./60.               # degrees per arcmin
DEG_PER_ASEC = DEG_PER_AMIN / 60.   # degrees per arcsec
RAD_PER_DEG = np.pi / 180.             # radians per degree
def ang_sep(ra1, dec1, ra2, dec2):
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
    
"""
Code testing John's function - works well!
    
ra1 = np.array([0.1, 0.2, 0.3])
dec1 = np.array([0.1, 0.2, 0.3])
ra2 = np.array([0.1, 0.2, 0.3, 0.4])
dec2 = np.array([0.1, 0.2, 0.3, 0.4])

for i in range(0, len(ra1)):
     matched_distances = ang_sep(ra1[i], dec1[i], ra2, dec2)
     print matched_distances
"""

chelsea_matching_rows = np.zeros_like(ra_jav_deg, dtype=int)

for i in range(0,len(ra_jav_deg)):
    matched_distances = ang_sep(ra_jav_deg[i],dec_jav_deg[i], ra_ch, dec_ch)
    indices = np.where(matched_distances <= 0.01)
    print 'ra_jav_deg row ', i 
    
    if (len(indices) == 1.0):
        print 'Obj from javelin coords: ', ra_jav_deg[i], dec_jav_deg[i]
        print 'Obj from Chelsea coords: ', ra_ch[indices], dec_ch[indices]
        chelsea_matching_rows[i] = int(indices[0])
    else :  chelsea_matching_rows[i] = 0
        

#print ra_ch[chelsea_matching_rows], dec_ch[chelsea_matching_rows]
#print log_tau_med[chelsea_matching_rows], log_sigma_med[chelsea_matching_rows]

template = "{0:8}{1:8}{2:8}{3:8}{4:12}{5:12}{6:12}{7:12}"
print template.format('ra_jav','dec_jav','ra_ch','dec_ch', 'tau_jav', 'log_10_tau_ch', 'sigma_jav', 'log_sig_ch')

#print("%5.4s %5.4s %5.4s %5.4s %5.4s %5.4s %5.4s %5.4s" %('r_j','d_j','r_c','d_c', 't_j', 't_c', 's_j', 's_c'))
for i in range(0,len(ra_jav_deg)):
    print("%7.4f  %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f" % (ra_jav_deg[i], \
    ra_ch[chelsea_matching_rows][i], dec_jav_deg[i], dec_ch[chelsea_matching_rows][i],\
    np.exp(tau_m[i]),np.power(10,log_tau_med[chelsea_matching_rows][i]), sigma_m[i],log_sigma_med[chelsea_matching_rows][i]))