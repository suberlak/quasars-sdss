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

#  shorter chain results:  
# 
#  javelin_chain_results_new7_sigma_tau.txt
#
#

javelin_chains = 'javelin_chain_results_new99_sigma_tau.txt'
javelin_results = np.genfromtxt('../QSO_try/'+javelin_chains, dtype="str")
"""
Important : the quasar chain  names are assumed to be in the form 
231408.02-011355.4_chain.dat
it doesn't really matter what is after the main quasar name, 
but since ra and dec coords are pulled out of the name string,
the first characters in the name in quasar_names_ra_dec  list matter. Otherwise 
it'll crash. 

"""
print 'working file', javelin_chains
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

template = "{0:10}{1:10}{2:8}{3:10}{4:10}{5:12}{6:12}{7:12}{8:12}"
print template.format('ra_jav','dec_jav','ra_ch','dec_ch', 'tau_jav', 'tau_ch', 'sigma_jav', 'sig_ch','sig_ratio' )

#print("%5.4s %5.4s %5.4s %5.4s %5.4s %5.4s %5.4s %5.4s" %('r_j','d_j','r_c','d_c', 't_j', 't_c', 's_j', 's_c'))

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
    col6[i] = np.power(10,log_tau_med[chelsea_matching_rows][i])
    col7[i] = sigma_m[i]
    col8[i] = np.power(10,log_sigma_med[chelsea_matching_rows][i])
    col9[i] = col7[i] / col8[i]
   
   
    
    #print("%7.4f  %7.4f %7.4f %7.3f %10.3f %10.3f %10.3f %10.3f %10.3f" % (col1, \
    #col2, col3, col4,col5, col6 , col7, col8, col9))
    
    
output = 'javelin_chelsea_comparison3.txt'
DAT= np.column_stack((col0,col1,col2,col3,col4,col5,col6,col7,col8,col9))

print '\n We are done matching javelin and Chelseas results. We matched ', \
len(ra_jav_deg), ' Quasars, with ', len(np.where(mask_mismatched == True)[0]),\
' those that required wider matching radius than 0.001 radians. Results are saved to',\
output," with columns 'ra_jav','dec_jav','ra_ch','dec_ch', 'tau_jav', 'tau_ch', 'sigma_jav', 'sig_ch','sig_ratio' "

# sort the DAT column accoring to tau_javelin 
newDAT=DAT[DAT[:,4].argsort()]

np.savetxt(output,newDAT,fmt="%s")
#, fmt=['%s','%7.3f', '%7.3f', '%7.3f', '%7.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f']