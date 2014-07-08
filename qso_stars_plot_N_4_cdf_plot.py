# Plotting module - takes the output of the analysis module, and plots whatever
# is needed
#
# A program to read in the output of   qso_stars_analysis.py  , 
# and plot the relevant global properties of the sample

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

# Files that contain diagnostics are in diag.list  file, in the working directory
# but because of their naming scheme, if the user selects the working directory, 
# it will find the relevant .npy files automatically 

dirs=['QSO','QSO_try', 'stars/0','stars_try']
directory=dirs[1]+'/'
dir_name=dirs[1]

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


# chi-sq cumulative lin-log 
N_obs=[4]
print 'Plotting chi-sq cumulative'
for number in N_obs :
    plt.clf()
    zeromask = (chisq[N==number] != 0.0)  # avoid values of zero
    logchisq = np.log10(chisq[N==number][zeromask])
    plt.hist(logchisq,bins=200, normed=True, cumulative=True)
    xmax = np.max(logchisq)
    #plt.xlim((0,0.5*xmax))
    plt.xlim((-2,2))
    nobjects=len(names)
    plt.title('Chi-sq cumulative distribution, N=4')
    plt.xlabel('log(chi-squared)')
    plt.ylabel('Probability')
    fname=dir_name+'_chisq_cum_N_linlog'+str(number)+'.png'
    

# overplot theoretical  cdf curve for k=3
df=3   # their name for k parameter ,  k = N-1 
x = np.linspace(chi2.ppf(0.0001,df),chi2.ppf(0.9999, df), 100)
# np.linspace(start,stop,number_of_samplings)

plt.plot(np.log10(x), chi2.cdf(x, df),'r-', lw=5, alpha=0.6, label='chi2 cdf')
plt.xlim((-2,2))
plt.savefig(fname)
