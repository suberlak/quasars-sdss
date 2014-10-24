# -*- coding: utf-8 -*-
"""
Created on Sat Oct 18 12:53:41 2014

@author: suberlak

A program to preprocess SDSS qso files, to provide the input required by javelin  

Reads - in files that have 16 columns :
eg 
 
1199742

u (0-2),g (3-5),r (6-8),i (9-11),z (12-14)  ( three columns per band  for MJD, mag, mag_err), 
and ra_median (15) , dec_median (16)

and save as 1 file per band 

"""

import numpy as np

dir_choice=['QSO_try2/', 'QSO_S82/', 'QSO_SDSS_JAV/']

dir_input=dir_choice[1]
dir_output=dir_choice[2]
names=np.loadtxt(dir_input+'out_good.list',dtype=str)

ra_list = np.zeros_like(names, dtype=float)
dec_list = np.zeros_like(names, dtype=float)

for j in range(len(names)):   #
 
    address=dir_input+names[j]
    data=np.loadtxt(address)
    ra_list[j]  = data[0,15]
    dec_list[j] = data[0,16]
    # cond_good_data=np.empty_like(data[:,0],dtype=bool)    
    
    mjds = [0,0,0,0,0]
    mags = [0,0,0,0,0]
    mag_errs = [0,0,0,0,0]
    
    mjd_cols = [0,3,6,9,12]
    mag_cols = [1,4,7,10,13]
    err_cols = [2,5,8,11,14]
    prefix = ['u','g','r','i','z']
    
    for l in range(len(mag_cols)):
        # checking whether a given mag does not have a bad mesaurement 
         
        col = mag_cols[l]
        #print col
        cond_good_data=np.ones_like(data[:,0],dtype=bool)
        #print 'Original data length ', len(data[:,col])
        for k in range(len(data)):
            if( data[k,col] < 1.0  or data[k,col] > 40):
                # print data[k,col], 'is < 1.0 or > 40'
            
                cond_good_data[k] = False
        #print 'Filtered data length', len(data[cond_good_data, col])
        mjds[l] = data[cond_good_data,mjd_cols[l]]     
        mags[l] = data[cond_good_data,mag_cols[l]] 
        mag_errs[l] = data[cond_good_data,err_cols[l]]
 
        # saving to a file result for each band 
        col_data = np.column_stack((mjds[l],mags[l],mag_errs[l]))
        out_name = dir_output+prefix[l]+'_'+names[j]+'.txt'
        np.savetxt(out_name,col_data,fmt="%s")
        print 'Saved', out_name 
        
    print 'Done with quasar', j+1, 'out of', len(names), '\n'
     
names_ra_dec = np.column_stack((names, ra_list, dec_list))
filename = dir_output + 'name_ra_dec_list.txt'
np.savetxt(filename,names_ra_dec,fmt="%s")
print 'Saved list of all names, ra, and dec for identification to ', filename   

# save  filename ,  ra  , and dec   as mapped file 
    
