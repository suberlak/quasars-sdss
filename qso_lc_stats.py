# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 17:57:13 2014

@author: suberlak

Retrieve statistics for each lightcurve:

avg_mag ( to see whether the observed quasar is bright or dim  )
avg_err (whether the error isn't too big )
N_obs = number of rows = lightcurve length 
t_obs = timespan of observations
avg_Nobs_day

Parts are taken from   qso_var.py which calculates statistics PER DAY of 
observations, such as  mjd,  avg_day_mag , avg_day_mag_err, N_obs / day  , 
chi-sq of the departure of the actual mag on a given day from the avg_day_mag
But that code gave   out_123456.7+123456.7.txt   , with hundreds of lines per 
quasar (corresponding to the number of days for which a quasar was observed ).
Here I want only 1 number for the entire lightcurve 

Parts are taken from  qso_stars_analysis.py ,  where the out_... . txt  files 
were taken in, and the output was : directory_all-rows_mjd_mag_err_N_chiqs.npy  has 
all the content of input dumped into one million-line file , and 
directry_name_timespan_nobs.npy ,  where the qso_name , t_obs,  and N_obs  were 
stored. Use part of this code to also calculate avg_mag and avg_err, as well as 
avg_Nobs_day. 

Parts about reading in the .npy files are from  qso_stars_plot.py 

out_123456.7+124356.7.txt  files have structure 
 mjd            <mjd>        <mjd_err>   N_obs/ day  chi-sq 
 54010.1898     19.5804      0.0597      4.0000      1.2290
 54024.1557     19.4404      0.0991      4.0000      7.7362
 54028.2808     19.4942      0.0561      4.0000      1.5182

123456.7+123456.7.dat files have structure 

mjd         mag   mag_err
53552.42200 20.16 0.33 0
53552.42845 19.93 0.29 0
53552.43491 19.84 0.28 0


"""

import numpy as np
import matplotlib.pyplot as plt 


# list of several directories I used : 

dirs=['QSO','QSO_try', 'stars/0','stars_try']

# user input : set the two parameters below to be the working directories 
# where files are  located by choosing appropriate number, or add another dir 
# to the list above 

directory=dirs[0]+'/'
dir_name=dirs[0]
lc_stat_files =np.loadtxt(directory+'out_good.list',dtype=str)


# make an array for storing total timespan of obs per object, 
# as well as the total number of nights per object
nobs_object = np.zeros_like(lc_stat_files).astype(float) # store how many nights per object
timespan_obs = np.zeros_like(lc_stat_files).astype(float) # store what was the total 
					# timespan of observations per object
lc_length=np.zeros_like(lc_stat_files).astype(float)
avg_N_day = np.zeros_like(lc_stat_files).astype(float)
avg_mag_ttl= np.zeros_like(lc_stat_files).astype(float)
avg_err_ttl= np.zeros_like(lc_stat_files).astype(float)
avg_mjd_diff= np.zeros_like(lc_stat_files).astype(float)
mean_time_betw_obs = np.zeros_like(lc_stat_files).astype(float)
# check how many total rows we have to create lists of appropriate size:
cond_multline=np.empty_like(lc_stat_files,dtype=bool)
#n_rows = 0

# http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python 

# MUST INCLUDE  A MASK FOR EMPTY FILES  !!! 


for i in range(len(lc_stat_files)-1):
    address=directory+lc_stat_files[i]   # lc diags 
    address1=directory + lc_stat_files[i][4:22]+'.dat'  # the original lc  
    data=np.loadtxt(address)     
    lc_length[i] = sum(1 for line in open(address1))
    nobs_object[i] = sum(1 for line in open(address))
    avg_mag_ttl[i] = np.mean(data[:,1])
    avg_err_ttl[i] = np.mean(data[:,2])
    avg_N_day[i] = np.mean(data[:,3])
    if len(data) != data.size :  # handling  multline data 
        #n_rows += len(data)
        #print len(data)
        cond_multline[i] = True 
        mjdrange = data[-1,0]-data[0,0]  # gives number of days between first and 
	# np.min(data[:,0]) - np.max(data[:,0])	  # last observation
        # same as the function in  the line above - we can assume that rows are 
        # chronological 
        timespan_obs[i] = mjdrange
        mjd_diff = np.zeros(len(data[:,0]))
        
        # the first slowing down diagnostic 
        for j in range(len(data)-1):
            mjd_diff[j] = data[j+1,0] - data[j,0]
        avg_mjd_diff[i] = np.mean(mjd_diff)          

        # this last diagnostic slows it down considerably... 
        data1 = np.loadtxt(address1)
        mjd_diff1 = np.zeros(len(data1[:,0]))
        for j in range(len(data1)-1):
            mjd_diff1[j] = data1[j+1,0] - data1[j,0]
        mean_time_betw_obs[i] = np.mean(mjd_diff1) 
        
    else:  
        #n_rows += 1
        #print '1'
        cond_multline[i] = False 
        timespan_obs[i] = 0
        avg_mjd_diff[i] = 0

print 'We have ', len(lc_stat_files), 'lightcurve diagnostic "out" files' 


file1 = 'qso_name_timespan_nobs_1.npy'
stats = [lc_stat_files,timespan_obs, nobs_object, lc_length,avg_N_day,avg_mag_ttl,avg_err_ttl,avg_mjd_diff,mean_time_betw_obs]
np.save(file1,stats)
print 'I saved to ', file1, ' statistics for each object : (0) name, (1) total ',\
'timespan of observations (final mjd - start mjd),  (2) total number of averaged',\
' observations (number of days on which quasar was observed),',\
' (3) lightcurve length (number of rows  in the original QSO data)',\
'(4) average number of measurements per day, (5) total average magnitude of a quasar:,',\
' mean along the <day_mag> column, (6) total average error : mean along the error column ',\
'(7) average mean number of days between observing days (from out files) , ',\
'(8) mean time between individual observations (from original files)  '

#allrows = np.load(file1)