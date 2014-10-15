# A short program to load the quasar file from CRTS, and plot its lightcurve,
# to get a feel for what they look like 


import numpy as np
import matplotlib.pyplot as plt
from random import randint
import os 

dirs=['QSO','QSO_try', 'stars/0','stars_try']
choosenumber=1
directory=dirs[choosenumber]+'/'
dir_name=dirs[choosenumber]
names_raw=np.loadtxt(directory+'file.list',dtype=str)


#####  REMOVING  EMPTY  FILES AND THOSE WITH ONLY ONE MESAUREMENT

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

print 'Using  files with more than one measurement...'


########  CHOOSING THOSE WHICH HAVE A LOT OF OBS FOR LIGHTCURVE ##############

x = randint(0,len(names[cond_multline])-1)
obj=names[cond_multline][x]
print 'plotting day averaged lightcurve for  quasar ', obj
address=directory+obj
data=np.loadtxt(address)

mjd = data[:,0]
mags = data[:,1]
errs = data[:,2]
days = data[:,0]
days = [int(day) for day in days]
days = np.unique(days)          # pulling out only the unique values 
avg_mags = np.zeros_like(days).astype(float)
mjd_arr = np.zeros_like(days).astype(float)

# run through days 
for i in range(len(days)):
    day = days[i]
    int_mjd = np.require(mjd,int)       # forcing mjd array -> integers
    condition = (int_mjd == day)        # finding where int(mjd) = day
    avgmag = np.average(mags[condition],weights=errs[condition])
    avg_mags[i] = avgmag
    mean_mjd = np.mean(mjd[condition])
    mjd_arr[i] = mean_mjd 
  

plt.plot(mjd_arr-50000,avg_mags,'o')
objtitle='Quasar '+obj+' lightcurve'
plt.title(objtitle)
plt.xlabel('mjd-50000')
plt.ylabel('mag') 
plt.show()
