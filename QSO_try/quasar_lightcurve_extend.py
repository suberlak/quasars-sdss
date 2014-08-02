# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 12:54:46 2014

@author: suberlak
"""

# A program to extend the working lightcurve ten times, and see whether javelin 
# would have any problems with that 

import numpy as np
import matplotlib.pyplot as plt

fname = '231408.02-011355.4.dat'
data = np.loadtxt(fname)

mjd = data[:,0]
mags = data[:,1]
errs = data[:,2]

print 'The working file is ', fname, ' it has ', len(mags), ' flux measurements',\
  '(rows)' 
print ' '
# ANALYSING THIS LIGHTCURVE
# Check  what are the intervals between obs,  and how many obs there are per day
days = [int(day) for day in mjd]
days = np.unique(days)          # pulling out only the unique values 
Nobs = np.zeros_like(days).astype(float)

for i in range(len(days)-1):
    day = days[i]
    int_mjd = np.require(mjd,int)       # forcing mjd array -> integers
    condition = (int_mjd == day)        # finding where int(mjd) = day
    N = float(len(mags[condition]))            # number of obs in that night 
    Nobs[i] = N
    
#print 'N :', Nobs[0:10]

print 'For this quasar we calculate what are the intervals between observations',\
   ' and in identical intervals we duplicate the flux measured '
print ' ' 
   
mjd_interval = np.zeros_like(mjd).astype(float)
for i in range(len(mjd)-1):
    mjd_interval[i] = mjd[i+1]-mjd[i]

#print 'len(mjd_interval)=', len(mjd_interval) , 'few first elements: ',mjd_interval[0:10]

# EXTENDING THE LIGHTCURVE  : NEW MAGS 
#mjd_ext  =  mjd + mjd_interval
#mjd_new = mjd 
mags_new=mags
errs_new=errs

times = 9
for i in range(times-1):
    
    #mjd_new=np.append(mjd_new,mjd_ext)
    #mjd_ext += mjd_interval   
    mags_new = np.append(mags_new, mags)
    errs_new = np.append(errs_new, errs)

# EXTENDING THE LIGHTCURVE : MJD'S  
shape = times*len(mjd)
mjd_test = np.zeros(shape)

for i in range(len(mjd)):
    mjd_test[i] = mjd[i]

j=i
# print 'j=', j

for k in range(times-1):
    for i in range(len(mjd)):
        mjd_test[j+i+1] = mjd_test[j+i] + mjd_interval[i]
        #print 'mjd_test', j+i+1, 'is', mjd_test[j+i+1], 'i= ', i
    j=i+j+1
    # print 'j=', j

    
print 'Plotting the new mjds ordered' 
number2 = np.arange(len(mjd_test))
plt.plot(number2,mjd_test, '.-')
plt.xlabel('element number')
plt.ylabel('mjd')
plt.show()

print 'Plotting the new mags ordered'
number = np.arange(len(mags_new))
plt.plot(number, mags_new)
plt.xlabel('element number')
plt.ylabel('magnitude')
plt.show()


# TEST WHETHER THE SAME LENGTH  
if (len(mjd_test) == len(errs_new) == len(mags_new)) == True : 
    print 'Well done : all three columns have the same length. ' , \
      'The length of the extended lightcurve is', len(mjd_test), ' elements', \
      'As compared to the original one, of ', len(mjd), ' elements'
    print ' '


print 'Plotting below chunks of the new lightcurve, to show that indeed the shape of new',\
  'variations is the same for each copy-and-paste chunk,',\
  'and the time spacing between each new flux measurement is the same as in',\
  'the original lightcurve'
fig, ax = plt.subplots(1,figsize=(12,28))

ax1.set_title('Extended lightcurve')

mjd_new = mjd_test

ax1 =plt.subplot(811)
ax1.plot(mjd_new[0:363], mags_new[0:363] ,'.-')
ax1.set_ylabel('magnitude ')

ax2 = plt.subplot(812)
ax2.plot(mjd_new[363:2*363], mags_new[363:2*363] ,'.-')
ax2.set_ylabel('magnitude ')

ax3 = plt.subplot(813)
ax3.plot(mjd_new[2*363:3*363], mags_new[2*363:3*363] ,'.-')
ax3.set_ylabel('magnitude ')

ax4 = plt.subplot(814)
ax4.plot(mjd_new[3*363:4*363], mags_new[3*363:4*363] ,'.-')
ax4.set_ylabel('magnitude ')

ax5 = plt.subplot(815)
ax5.plot(mjd_new[4*363:5*363], mags_new[4*363:5*363] ,'.-')
ax5.set_ylabel('magnitude ')

ax6 = plt.subplot(816)
ax6.plot(mjd_new[5*363:6*363], mags_new[5*363:6*363] ,'.-')
ax6.set_ylabel('magnitude ')


ax6.set_xlabel('mjd')
plt.show()

fname1 = fname[:-4]+'_extended.dat'
 
print 'Saving the resultant lightcurve with ', len(mjd_test), ' elements to ',\
   'a file ', fname1
   
np.savetxt(fname1, np.column_stack((mjd_new,mags_new,errs_new)),fmt='%11.4f')
   
    



