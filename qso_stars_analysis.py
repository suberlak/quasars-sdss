# Analysis module : pulls all the output files together and creates the input 
# for the plotting module 
# 
# A combined program: merger of qso_an, and stars_an  
# to analyse either CRTS quasars, preprocessed by qso_var , 
# or standard stars, preprocessed by str_var

# after performing all the averaging calculations on files located in the 
# directory with qso_var routine , we have output files as in out.list

# http://stackoverflow.com/questions/89228/calling-an-external-command-in-python

# with structure 

# ave_MJD_of_night  | ave mag |  weighted_error |  N_obs_per_night | chi-sq 

# read in all the files from out.list , plot global distributions of various 
# parameters


import numpy as np
import matplotlib.pyplot as plt 


# list of several directories I used : 

dirs=['QSO','QSO_try', 'stars/0','stars_try']

# user input : set the two parameters below to be the working directories 
# where files are  located by choosing appropriate number, or add another dir 
# to the list above 

directory=dirs[0]+'/'
dir_name=dirs[0]
names=np.loadtxt(directory+'out.list',dtype=str)

# make an array for storing total timespan of obs per object, 
# as well as the total number of nights per object
nobs_object = np.zeros_like(names).astype(float) # store how many nights per object
timespan_obs = np.zeros_like(names).astype(float) # store what was the total 
					# timespan of observations per object

# check how many total rows we have to create lists of appropriate size:
cond_multline=np.empty_like(names,dtype=bool)
n_rows = 0
for i in range(len(names)):
    address=directory+names[i]
    data=np.loadtxt(address)
       
    if len(data) != data.size :  # handling  multline data 
        n_rows += len(data)
        #print len(data)
        cond_multline[i] = True 
        mjdrange = data[-1,0]-data[0,0]  # gives number of days between first and 
	# np.min(data[:,0]) - np.max(data[:,0])	  # last observation
        # same as the function in  the line above - we can assume that rows are 
        # chronological 
        timespan_obs[i] = mjdrange
        nobs_object[i]= len(data)
    else:    # handling singleline data 
        n_rows += 1
        #print '1'
        cond_multline[i] = False 
        timespan_obs[i] = 0
        nobs_object[i]= 1

print 'We have ', n_rows, 'rows total, in ', len(names), 'out files' 


# make appropriate arrays 

mjd =  np.empty(0,dtype=float)
mag = np.empty(0,dtype=float)
err =  np.empty(0,dtype=float)
N =  np.empty(0,dtype=float)
chisq =  np.empty(0,dtype=float)


# fill in arrays with all the data from multiline 

for obj in names[cond_multline] : 
    address=directory+obj 
    data=np.loadtxt(address) 
    mjd = np.append(mjd,data[:,0])
    mag = np.append(mag,data[:,1])
    err = np.append(err,data[:,2])
    N = np.append(N,data[:,3]) 
    chisq = np.append(chisq,data[:,4])

# add the data from singleline 

for obj in names[np.logical_not(cond_multline)] : 
    address=directory+obj 
    data=np.loadtxt(address) 
    mjd = np.append(mjd,data[0])
    mag = np.append(mag,data[1])
    err = np.append(err,data[2])
    N = np.append(N,data[3])
    chisq = np.append(chisq,data[1])


file1 = dir_name + '_all-rows_mjd_mag_err_N_chisq.npy'
allrows = [mjd,mag,err,N,chisq]
np.save(file1,allrows)
print 'I saved to ', file1, 'all data for all objects : obs_date, mean mag,',\
' error, N, chi-squared, with ', len(mjd), 'rows' 

file2 = dir_name + '_name_timespan_nobs.npy'
stats = [names,timespan_obs, nobs_object]
np.save(file2,stats)
print 'I saved to ', file2, ' statistics for each object : name, total ',\
'timespan of observations, and total number of averaged observations, with ',\
len(names), 'rows'

#
#
# plotting is done with qso_stars_plot.py  
#
#


