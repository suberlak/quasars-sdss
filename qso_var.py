# Calculate average magnitudes for each night, per object.
# We also calculate error by the  weighted mean of photometric error of 
# individual measurements from a given night.  
# Save the resulting average magnitude, weighted error, number of obs per night,
# and the reduced chi2 to a file named out_[inputfile].txt

import numpy as np
np.seterr(invalid='raise')

names=np.loadtxt('QSO_try/file.list',dtype=str)

master = np.zeros((len(names),5))
# open the first one to try

for obj in names:
    address='QSO_try/'+obj
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
    Nobs = np.zeros_like(days).astype(float)
    print ' '
    print 'For Quasar', obj

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
        if N == 1.0 : 
            chi2dof = 1.0 
        else :
            chi2dof = np.sum(weights*(np.power((mags[condition]-avgmag),2.0))) / (N-1.0)
        chi2arr[i] = chi2dof
        print 'i = ', i, 'On day MJD', day, 'N obs=', N, 'avgmag=', avgmag, 'avg_err=',error, \
        'chi2dof=',chi2dof
    
    print '  '
    # save output of averaging of each file to a separate file 
    name_out='QSO_try/out_'+obj[:18]+'.txt'
    np.savetxt(name_out, np.column_stack((avg_mags,avg_err,Nobs,chi2arr)),fmt='%11.4f')


