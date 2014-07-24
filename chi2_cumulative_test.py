# A test of my tools : 

# given the file which contains various measurements of quasars, each line stands
# for a separate day, with exactly four observations per night 


# calculate the average value of four measurements, 
# and thus calculate chi-squared for those days 

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

data = np.loadtxt('chi2testN4.dat')

err = data[:,0]
mag1 = data[:,1]
mag2 = data[:,2]
mag3 = data[:,3]
mag4 = data[:,4]

avg_mags = np.zeros_like(err).astype(float)
avg_err = np.zeros_like(err).astype(float)
chisq = np.zeros_like(err).astype(float)

for i in range(len(err)):
    avg = np.average(data[i,1:4])
    avg_mags[i] = avg 
    w = 1.0 / (err[i] * err[i])
    avg_err[i] = err[i] / np.sqrt(3.0)  
    chisq[i] = w * ( (mag1[i] - avg)**2.0 +  (mag2[i] - avg)**2.0 + (mag3[i] - avg)**2.0 + (mag4[i] - avg)**2.0 )

# chi-sq cumulative lin-log 
N_obs=[4]
print 'Plotting chi-sq cumulative'
for number in N_obs :
    plt.clf()
    #zeromask = (chisq != 0.0)  # avoid values of zero
    logchisq = np.log10(chisq)
    plt.hist(logchisq,bins=200, normed=True, cumulative=True)
    xmax = np.max(logchisq)
    #plt.xlim((0,0.5*xmax))
    plt.xlim((-2,2))
    plt.title('Chi-sq cumulative test distribution, N=4')
    plt.xlabel('log(chi-squared)')
    plt.ylabel('Probability')
    fname='test_chisq_cum_N_linlog_'+str(number)+'.png'
    

# overplot theoretical  cdf curve for k=3
df=3   # their name for k parameter ,  k = N-1 
x = np.linspace(chi2.ppf(0.0001,df),chi2.ppf(0.9999, df), 100)
# np.linspace(start,stop,number_of_samplings)

plt.plot(np.log10(x), chi2.cdf(x, df),'r-', lw=5, alpha=0.6, label='chi2 cdf')
plt.xlim((-2,2))
plt.savefig(fname)
