# after performing all the averaging calculations on files located in the 
# directory with qso_var routine , we have output files as in out.list

# http://stackoverflow.com/questions/89228/calling-an-external-command-in-python

# with structure 

# ave_MJD_of_night  | ave mag |  weighted_error |  N_obs_per_night | chi-sq 

# read in all the files from out.list , plot global distributions of various 
# parameters


import numpy as np
import matplotlib.pyplot as plt 

directory='QSO_try/'
dir_name='QSO_try'
names=np.loadtxt(directory+'out.list',dtype=str)

# check how many total rows we have to create lists of appropriate size:
cond_multline=np.empty_like(names,dtype=bool)
n_rows = 0
for i in range(len(names)):
    address=directory+names[i]
    data=np.loadtxt(address)
    if len(data) != data.size :
        n_rows += len(data)
        print len(data)
        cond_multline[i] = True 
    else:
        n_rows += 1
        print '1'
        cond_multline[i] = False 

print 'We have ', n_rows, 'rows total, in ', len(names), 'out files' 

# make appropriate arrays 

avg_mag = np.zeros(n_rows).astype(float)
err_m = np.zeros(n_rows).astype(float)
chi_sq = np.zeros(n_rows).astype(float)
N = np.zeros(n_rows).astype(float)
mjd = np.zeros(n_rows).astype(float)

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

# plot the distributions
# http://bespokeblog.wordpress.com/2011/07/11/basic-data-plotting-with-matplotlib-part-3-histograms/

# mean error 
plt.clf()
plt.hist(err,bins=30)
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

fig1 = plt.figure()
plt.plot(err,mag,'.r')
plt.xlabel('Error')
plt.ylabel('Magnitude')
nbins =200
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







