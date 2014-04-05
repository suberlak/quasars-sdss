# Calculate average magnitudes for each night, per object.
# We also calculate error by the  weighted mean of photometric error of 
# individual measurements from a given night.  
# Save the resulting average magnitude, weighted error, number of obs per night,
# and the reduced chi2 to a file named out_[inputfile].txt

import numpy as np
np.seterr(invalid='raise')
import warnings 
import os 

directory='QSO/'
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
cond_notempty=np.empty_like(names,dtype=bool)
for i in range(len(names)):
    if os.stat(directory+names[i])[6]!=0: 
        cond_notempty[i] = True
    else:
        cond_notempty[i] = False

# how many nonempy files we have 

num_empty=len(names)- np.sum(cond_notempty)
num_notempty=np.sum(cond_notempty)

print 'Out of', len(names), 'files, we have', num_notempty, 'not-empty files', \
'and therefore, ', num_empty, 'empty files'

counter=0
# using only notempty files
for obj in names[cond_notempty]:
    address=directory+obj
    data=np.loadtxt(address)

    averages=np.zeros(shape=(len(data),3))

# ADD HERE A CHECK THAT ADDRESSES FILES THAT HAVE ONLY ONE MEASUREMENT FOR A 
# QUASAR IN A DIFFERENT WAY: FOR   names[cond_notempty][180]  , I.E.
# QSO/000303.32+001019.6.dat   ,WE HAVE ONLY   data=[53643.27065, 21.90, 0.41, 0]
# SO THAT SAYING data[:,0]  is meaningless, with just one row / column there is 
# distinction into rows / columns 
 

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
    
    print '  '
    # save output of averaging of each file to a separate file 
    name_out=directory+'out_'+obj[:18]+'.txt'
    np.savetxt(name_out, np.column_stack((mjd_arr,avg_mags,avg_err,Nobs,\
    chi2arr)),fmt='%11.4f')
    counter += 1    


# print which files were empty and not used 
print 'Files that were empty:', names[np.logical_not(cond_notempty)] 
     # np.logical_not negates the boolean values 

