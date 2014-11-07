# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 18:04:33 2014

@author: suberlak

A program to take the lightcurve stats from  qso_lc_stats.py ,  
and the output of comparison of javelin and Chelsea's result for each quasar. 
It pulls from the stats table all the relevant stats for the subset of 
quasars on which javelin was run, to plot the log(tau_difference) vs log(sigma_difference),
and colour - code the third dimension according to the value of the investigated
stat parameter for a given quasar,  in a quest of looking for a correlation 
that may show to what types of quasar lightcurves javelin is failing, and if
it is  failing systematically. 


"""
from math import e, pi, isinf
import numpy as np 
import matplotlib.pyplot as plt


"""
File with results of Chelsea's analysis
"""

dir_in = 's82drw/'
dir_out = 'QSO_SDSS_analysis/'
band=['u','g','r','i','z']
band_pos = 0
SDSS_javelin_files = 'javelin_SDSS_chelsea_comparison_'+band[band_pos]+'_band.txt'
javelin_chelsea_comparison =dir_in +SDSS_javelin_files
javch =  np.loadtxt(javelin_chelsea_comparison, dtype='str')

"""
File with the SDSS quasars statistics: one file for all bands 
"""

filename = 'QSO_SDSS_analysis/'+ 'SDSS_lc_stats_multiband.txt'
qso_stats =  np.loadtxt(filename) 

print 'We have ugriz stats for SDSS data from file ', filename, 'of length', len(qso_stats)
print 'We also have a comparison of Chelsea and javelin  SDSS fits for ',len(javch), ' quasars.' 

print 'Retrieving statistics for quasars from ', band[band_pos], 'band'


"""
Data from comparison of Javelin and Chelsea's fits 
"""
qso_name = javch[:,0]
ra_jav = javch[:,1].astype(np.float)
ra_ch= javch[:,2].astype(np.float)
dec_jav =javch[:,3].astype(np.float)
dec_ch = javch[:,4].astype(np.float)
tau_jav = javch[:,5].astype(np.float)
tau_ch = javch[:,6].astype(np.float)
sigma_jav =javch[:,7].astype(np.float)
sigma_ch = javch[:,8].astype(np.float)
sig_rat = javch[:,9].astype(np.float)


log_tau_ratio = np.log10(tau_jav / tau_ch)
log_sigma_ratio = np.log10(sigma_jav/sigma_ch)

print 'We are now retrieving stats data for each matched quasar...'

# initiate stats arrays

lc_length = np.zeros_like(ra_jav)
avg_mag_ttl = np.zeros_like(ra_jav)
avg_err_ttl = np.zeros_like(ra_jav)


print 'There are ', len(qso_name), ' quasars to match with stat data ' 

p = band_pos+1

for i in range(len(qso_name)):
    index = np.where(qso_stats[:,0] == qso_name[i].astype(float))[0][0]
    lc_length[i] = qso_stats[index,3*p]
    avg_mag_ttl[i] = qso_stats[index,3*p+1]
    avg_err_ttl[i] = qso_stats[index,3*p+2]
     
            
# good_indices= np.where(map(isinf,log_tau_ratio)==False)


# get rid of infinities : 
#for i in range(len(log_sigma_ratio)):
#    if map(isinf,log_sigma_ratio)[i] == True : 
#        log_sigma_ratio[i] = -1e20
#
#for i in range(len(log_tau_ratio)):
#    if map(isinf,log_tau_ratio)[i] == True : 
#        log_tau_ratio[i] = 1e20

'''
Ensuring that there are no infinite ratios
'''    
st = np.asarray(map(isinf,log_sigma_ratio),dtype=bool) 
lt = np.asarray(map(isinf,log_tau_ratio),dtype=bool) 
res = st+lt

ist = np.where(st == True)
ilt= np.where(lt ==True)

ires = np.where(res == True) # this np.boolean array should contain indices 
     # for which both st and lt  are True, i.e. either log_sigma_ratio  or log_tau_ratio  is infinite 

gi = -res   # good_indices 

non_inf = len(np.where(gi == True)[0])


print 'Out of ', len(ra_jav),' rows, we have ', non_inf, ' of those that do not',\
' have any infinities, and only those are saved in the txt output '

'''

Restricting the sample to choose only few...

Want: log_tau_ratio   , and  log_sigma_ratio   to be within  some bounds

available stats for all qso's:

- lc_length
- avg_mag_ttl
- avg_err_ttl

'''

cond = np.zeros(len(ra_jav),dtype=bool)

for i in range(len(ra_jav)):
    if (log_tau_ratio[i] < 0) and (log_tau_ratio[i] > -0.5) and (log_sigma_ratio[i] < 0) and (log_sigma_ratio[i] > -0.5): 
        cond[i] = True

print len(np.where(cond == True)[0]), 'Quasars meet the selection criteria'

# Finding the longest lc within those selected quasars:
index = np.where(lc_length == max(lc_length[cond]))
print '\nFrom those quasars, the longest is', qso_name[index][0], 'with', lc_length[index][0], 'rows.'
print 'Its <mag> =', avg_mag_ttl[index][0], 'and <avg_err_ttl>=', avg_err_ttl[index][0]

# Finding the brightest lc within the selected ones
index = np.where(avg_mag_ttl == min(avg_mag_ttl[cond]))
print '\nFrom those quasars, the brightest is', qso_name[index][0], 'with avg_mag',avg_mag_ttl[index][0] 
print 'Its length is', lc_length[index][0] , 'and <avg_err_ttl>=', avg_err_ttl[index][0]


'''

Some plotting ...

'''

# final indices for things that are both without  infinite ratios, and fulfill the
# selection criteria above 

new = np.zeros(len(ra_jav),dtype=bool)
for i in range(len(cond)): 
    new[i] = gi[i] and cond[i]
    
# Plot logs of ratios
    
x=log_tau_ratio[gi]
y=log_sigma_ratio[gi]
plt.plot(x,y,'.r')
nbins =1000
H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlim((-3,1))
plt.ylim((-2,1))
title = 'SDSS quasars '+band[band_pos]+' band Javelin vs Chelsea SDSS '+band[band_pos]+' band  '
plt.title(title)
plt.xlabel('log(tau_jav / tau_ch)')
plt.ylabel('log(sigma_jav / sigma_ch)')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname3='test_sdss_'+band[band_pos]+'_band-log_sigma_ratio_vs_log_tau_ratio2.png'
plt.savefig(fname3)


# Plot lc_length vs avg_mag... 
plt.clf()
title = 'SDSS quasars '+band[band_pos]+' selected '+ str(len(avg_mag_ttl[new]))
plt.title(title)
plt.xlabel('<magnitude>')
plt.ylabel('Lightcurve length [days]')
plt.scatter(avg_mag_ttl[new],lc_length[new],c=avg_err_ttl[new],cmap=plt.cm.coolwarm)
fname1='test_sdss_'+band[band_pos]+'_band-SEL_avg_mag_vs_lc_length1.png'
plt.savefig(fname1)


# plot lc_length  vs   avg_mag  as histogram 

x=avg_mag_ttl[new]
y=lc_length[new]
plt.plot(x,y,'.r')
nbins=100
H, xedges,yedges = np.histogram2d(x,y,bins=nbins)
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H)
fig3 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
#plt.xlim((-3,1))
#plt.ylim((-2,1))
title = 'SDSS quasars '+band[band_pos]+' selected '+ str(len(avg_mag_ttl[new]))
plt.title(title)
plt.xlabel('<magnitude>')
plt.ylabel('Lightcurve length [days]')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fname3='test_sdss_'+band[band_pos]+'_band-sel_avg_mag_vs_lc_length.png'
plt.savefig(fname3)

'''
Save all data 
'''


DAT = np.column_stack((qso_name,ra_jav,ra_ch,tau_jav,tau_ch,sigma_jav, sigma_ch,\
 lc_length, avg_mag_ttl, avg_err_ttl)) 

# sort the DAT column accoring to col(7) : log_tau_ratio : 
 
#newDAT=DAT[DAT[:,6].argsort()]

file_all = dir_out + 'qso_SDSS_jav_chelsea_and_stats_'+band[band_pos]+'_band_ALL.txt'
np.savetxt(file_all,DAT,fmt="%s")
print '\nSaving all the Chelsea vs Javelin SDSS comparison enriched with statistical',\
' data (even those where log_tau_ratio would lead to infinity), for ',len(ra_jav) ,'quasars',\
' for ', band[band_pos] ,' band to file ' , file_all 


'''
Save only selected data
'''

sel_DAT = np.column_stack((qso_name[new],ra_jav[new],ra_ch[new],tau_jav[new],tau_ch[new],\
sigma_jav[new], sigma_ch[new], lc_length[new], avg_mag_ttl[new], avg_err_ttl[new]))

file_sel = dir_out +  'qso_SDSS_jav_chelsea_compared_stats_'+band[band_pos]+'_band_SEL.txt'
np.savetxt(file_sel,sel_DAT,fmt="%s")
print '\nSaving all the Chelsea vs Javelin SDSS comparison enriched with statistical',\
' data (even those where log_tau_ratio would lead to infinity), for ',len(qso_name[new]) ,\
'quasars for ', band[band_pos] ,' band to file ' , file_sel 


mag = avg_mag_ttl[new]
avg_err_ttl[new]
lc = lc_length[new]
ind = np.where(mag < 18)
np.where(lc[ind] > 80)

lc[ind][41]

indices = [6,41]
for i in indices:
    print qso_name[new][ind][i],lc[ind][i],mag[ind][i],err[ind][i],tau_jav[new][ind][i],\
    tau_ch[new][ind][i],sigma_jav[new][ind][i],sigma_ch[new][ind][i]