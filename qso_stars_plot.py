# Plotting module - takes the output of the analysis module, and plots whatever
# is needed
#
# A program to read in the output of   qso_stars_analysis.py  , 
# and plot the relevant global properties of the sample

import numpy as np
import matplotlib.pyplot as plt

# Files that contain diagnostics are in diag.list  file, in the working directory
# but because of their naming scheme, if the user selects the working directory, 
# it will find the relevant .npy files automatically 

dirs=['QSO','QSO_try', 'stars/0','stars_try']
directory=dirs[2]+'/'
dir_name=dirs[2]

# always refer to qso_stars_analysis.py  for the names of the two diagnostic files

file1 = dir_name + '_all-rows_mjd_mag_err_N_chisq.npy'
file2 = dir_name + '_name_timespan_nobs.npy'

allrows=np.load(file1)
mjd= allrows[0]
mag= allrows[1]
err=allrows[2]
N = allrows[3]
chisq=allrows[4]

stats = np.load(file2) 
names=stats[0]
timespan_obs=stats[1]
timespan_obs=timespan_obs.astype(np.float)
nobs_object=stats[2]
nobs_object=nobs_object.astype(np.float)

# plot the distributions
# http://bespokeblog.wordpress.com/2011/07/11/basic-data-plotting-with-matplotlib-part-3-histograms/

# mean error 
plt.clf()
plt.hist(err,bins=40)
plt.title('Error distribution')
plt.xlabel('Mean weighted error')
plt.ylabel('Frequency')
fname=dir_name+'_mean_err.png'
plt.savefig(fname)


# number of observations

plt.clf()
plt.hist(N,bins=30)
plt.title('Distribution of the number of observations per night per object')
plt.xlabel('Number of observations')
plt.ylabel('Frequency')
fname=dir_name+'_N_obs.png'
plt.savefig(fname)

# average magnitude 

plt.clf()
plt.hist(mag,bins=30)
plt.title('Average magnitude distribution')
plt.xlabel('Magnitude')
plt.ylabel('Frequency')
fname=dir_name+'_mag.png'
plt.savefig(fname)

# error vs  average magnitude  as a colour- coded histogram
# http://oceanpython.org/2013/02/25/2d-histogram/
print 'plotting error vs average mag colour-coded histogram'
fig1 = plt.figure()
plt.plot(err,mag,'.r')
plt.xlabel('Error')
plt.ylabel('Magnitude')
nbins =100
H, xedges,yedges = np.histogram2d(err,mag,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlabel('Error')
plt.ylabel('Magnitude')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname=dir_name+'_mag_err.png'
plt.savefig(fname)


# chi-sq cumulative log-log
N_obs=np.unique(N)
print 'Plotting chi-sq cumulative'
for number in N_obs :
    plt.clf()
    zeromask = (chisq[N==number] != 0.0)  # avoid values of zero
    logchisq = np.log10(chisq[N==number][zeromask])
    plt.hist(logchisq,bins=200,  log=True, normed=True, cumulative=True)
    xmax = np.max(logchisq)
    plt.xlim((0,0.5*xmax))
    plt.title('Chi-sq cumulative distribution')
    plt.xlabel('log(chi-squared)')
    plt.ylabel('Probability')
    fname=dir_name+'_chisq_cum_N_loglog'+str(number)+'.png'
    plt.savefig(fname)

# chi-sq cumulative lin-log 
N_obs=np.unique(N)
print 'Plotting chi-sq cumulative'
for number in N_obs :
    plt.clf()
    zeromask = (chisq[N==number] != 0.0)  # avoid values of zero
    logchisq = np.log10(chisq[N==number][zeromask])
    plt.hist(logchisq,bins=200, normed=True, cumulative=True)
    xmax = np.max(logchisq)
    plt.xlim((0,0.5*xmax))
    plt.title('Chi-sq cumulative distribution')
    plt.xlabel('log(chi-squared)')
    plt.ylabel('Probability')
    fname=dir_name+'_chisq_cum_N_linlog'+str(number)+'.png'
    plt.savefig(fname)


# number of nights per quasar vs timespan of obs as a colour- coded histogram
# http://oceanpython.org/2013/02/25/2d-histogram/
print 'plotting number of nights vs timespan of obs per quasar as',\
' a  colour-coded histogram'
fig1 = plt.figure()
plt.plot(timespan_obs,nobs_object,'.r')
plt.xlabel('Timespan[days]')
plt.ylabel('Number of observations')
nbins =200.0
H, xedges,yedges = np.histogram2d(timespan_obs,nobs_object,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlabel('Timespan[days]')
plt.ylabel('Number of observations')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname=dir_name+'_N_t-obs200.png'
plt.savefig(fname)

