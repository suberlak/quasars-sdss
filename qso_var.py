# Calculate average magnitudes for each night, per object.
# We also calculate error by the  weighted mean of photometric error of 
# individual measurements from a given night.  
# Save the resulting average magnitude, weighted error, number of obs per night,
# and the reduced chi2 to a file named out_[inputfile].txt

import numpy as np
np.seterr(invalid='raise')
import warnings 
import os 

directory='QSO_try/'
names=np.loadtxt(directory+'file.list',dtype=str)


# ignore warnings
#
#warnings.filterwarnings("error")

#def fxn():
#    warnings.warn("deprecated", DeprecationWarning)
#
#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")
#    fxn()


# instead of ignoring warnings of particular type, check which files are empty
# and count how many nonempty files we have which meet that criterion

cond_notempty=np.empty_like(names,dtype=bool)
for i in range(len(names)):
    if os.stat(directory+names[i])[6]!=0: 
        cond_notempty[i] = True
    else:
        cond_notempty[i] = False 

num_empty=np.sum(np.logical_not(cond_notempty))
num_notempty=np.sum(cond_notempty)

print 'Out of', len(names), 'files, we have', num_notempty, 'not-empty files', \
'and therefore, ', num_empty, 'empty files'

# check if there are any files with only one line of measurement
cond_multline=np.empty_like(names,dtype=bool)

for i in range(len(names)):
    address=directory+names[i]
    data=np.loadtxt(address)
    
    if len(data) == 4.0 :
         cond_multline[i] = False
    else:
         cond_multline[i] = True

num_multline = np.sum(cond_multline) - num_empty
num_singleline = np.sum(np.logical_not(cond_multline))

# print 'After  checking for single-line measurements, we have '
print 'We have',  num_multline, ' objects with more than one measurement, and',\
 num_singleline, ' objects with only one measurement'


# create a mask for notempty, multiline files 

cond_multi_notempty = np.empty_like(names,dtype=bool)
for i in range(len(names)) : 
    if  ( cond_notempty[i] == True ) and (cond_multline[i] == True) :
        cond_multi_notempty[i] = True
    else:
        cond_multi_notempty[i] = False 
num= np.sum(cond_multi_notempty)

print 'Out of ', len(names),'files total, we have ', num, \
'objects with more than one measurement'
print '  ' 

print 'Performing calculations on files with more than one measurement...'

counter=0
# using only notempty multiline files
for obj in names[cond_multi_notempty]:
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



# ADD HERE A CHECK THAT ADDRESSES FILES THAT HAVE ONLY ONE MEASUREMENT FOR A 
# QUASAR IN A DIFFERENT WAY: FOR   names[cond_notempty][180]  , I.E.
# QSO/000303.32+001019.6.dat   ,WE HAVE ONLY   data=[53643.27065, 21.90, 0.41, 0]
# SO THAT SAYING data[:,0]  is meaningless, with just one row / column there is 
# distinction into rows / columns 


# print which files were empty and not used 
print 'Files that were empty:', names[np.logical_not(cond_notempty)] 
     # np.logical_not negates the boolean values 

