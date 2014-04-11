# Calculate average magnitudes for each night, per object.
# We also calculate error by the  weighted mean of photometric error of 
# individual measurements from a given night.  
# Save the resulting average magnitude, weighted error, number of obs per night,
# and the reduced chi2 to a file named out_[inputfile].txt
#
# output file structure:
# ave_MJD_of_night  | ave mag |  weighted_error |  N_obs_per_night | chi-sq 

import numpy as np
np.seterr(invalid='raise')
import warnings 
import os 

directory='QSO_try/'
names_raw=np.loadtxt(directory+'file.list',dtype=str)


# check which files are empty
# and count how many nonempty files we have which meet that criterion

cond_notempty=np.empty_like(names_raw,dtype=bool)
for i in range(len(names_raw)):
    if os.stat(directory+names_raw[i])[6]!=0: 
        cond_notempty[i] = True
    else:
        cond_notempty[i] = False 

num_empty=np.sum(np.logical_not(cond_notempty))
num_notempty=np.sum(cond_notempty)

print 'Out of', len(names_raw), 'files, we have', num_notempty, 'not-empty files', \
'and therefore, ', num_empty, 'empty files'

print 'Files that were empty:', names_raw[np.logical_not(cond_notempty)]

names=names_raw[cond_notempty]

# of all notempty files 
# check if there are any files with only one line of measurement

cond_multline=np.empty_like(names,dtype=bool)

for i in range(len(names)):
    address=directory+names[i]
    data=np.loadtxt(address)
    
    if data.size == 4.0 :
         cond_multline[i] = False
    else:
         cond_multline[i] = True

num_multline = np.sum(cond_multline)
num_singleline = np.sum(np.logical_not(cond_multline))

# print 'After  checking for single-line measurements, we have '
print 'We have',  num_multline, ' objects with more than one measurement, and',\
 num_singleline, ' objects with only one measurement'


print '  ' 

print 'Performing calculations on files with more than one measurement...'


counter=0
# using only notempty multiline files
for obj in names[cond_multline]:
    address=directory+obj
    data=np.loadtxt(address)

    averages=np.zeros(shape=(len(data),3))

    mjd = data[:,0]
    mags = data[:,1]
    errs = data[:,2]
    days = data[:,0]
    days = [int(day) for day in days]
    days = np.unique(days)          # pulling out only the unique values 

    avg_mags = np.zeros_like(days).astype(float)
				    # preparing an array in advance to store 
    				    # is way more efficient
    avg_err = np.zeros_like(days).astype(float)
    chi2arr = np.zeros_like(days).astype(float)
    mjd_arr = np.zeros_like(days).astype(float)
    Nobs = np.zeros_like(days).astype(float)
    print 'obj= ',counter, 'For Quasar', obj

    # loop through days calculating mean, etc. 
    for i in range(len(days)):
        day = days[i]
        int_mjd = np.require(mjd,int)       # forcing mjd array -> integers
        condition = (int_mjd == day)        # finding where int(mjd) = day
        N = float(len(mags[condition]))            # number of obs in that night 
        Nobs[i] = N
        avgmag = np.average(mags[condition],weights=errs[condition])
        weights=1.0 / ( errs[condition] * errs[condition]) 
        avg_mags[i] = avgmag
        error = 1.0 / np.sqrt(np.sum(weights))
        avg_err[i] = error
        mean_mjd = np.mean(mjd[condition])
        mjd_arr[i] = mean_mjd 
        chi2 = np.sum(weights*(np.power((mags[condition]-avgmag),2.0))) 
        chi2arr[i] = chi2
        # print 'i = ', i, 'On day MJD', day, 'N obs=', N, 'avgmag=', avgmag, \
        # 'avg_err=',error, 'chi2=',chi2
    
    # save output of averaging of each file to a separate file 
    name_out=directory+'out_'+obj[:18]+'.txt'
    print 'Saved', name_out
    np.savetxt(name_out, np.column_stack((mjd_arr,avg_mags,avg_err,Nobs,\
    chi2arr)),fmt='%11.4f')
    counter += 1    
    print '  ' 

# using here only notempy single line files ( I  don't think that would have any
# practical use, but is necessary for consistent treatment of the dataset

print 'Now calculating single line measurements'
counter = 0
for obj in names[np.logical_not(cond_multline)]:
    address=directory+obj
    data=np.loadtxt(address)

    mjd = data[0]
    mag = data[1]
    err = data[2]
    Nobs=1
    print 'obj= ',counter, 'For Quasar', obj
    weight=1.0 / ( err * err) 
    error = 1.0 / np.sqrt(weight)
    chi2 = 1.0
      
    # save output 
    name_out=directory+'out_'+obj[:18]+'.txt'
    print 'Saved', name_out
    np.savetxt(name_out, np.column_stack((mjd,mag,error,Nobs,\
    chi2)),fmt='%11.4f')
    counter += 1    
    print '  ' 

# create a list of all output files, for other programs to use 
# command = 'ls '+directory+'out_* '+'> '+directory+'out.list'
# os.system(command)
# 
# unhelpful : because it includes the directory before the file name! 
# i think I'd have to create a list of length of all non-empty files , 
# and fill it with names as they are created , and then dump all of that 
# to a file called  file.list
#
# but I'm not sure if there would not be any troubles with formatting 

# for now, to get the list of output files, do in the shell 

# cd QSO/
# touch out.list
# ls out* > out.list
# cd ..


