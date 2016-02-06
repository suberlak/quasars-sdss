# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 12:22:55 2015

@author: suberlak
"""
import fsps
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 
sns.set_context("poster") 

##  PROBLEM 1

Gal_age = 12.5 # Age of Galaxy  in  Gyrs
Uni_age = 13.8 # Age of the Universe in Gyr 
sf_start = Uni_age - Gal_age
sp = fsps.StellarPopulation(compute_vega_mags=False, sfh=1, zmet=20,
                                dust_type=2, dust2=0.2, sf_start=sf_start)
                                
sdss_filters = fsps.find_filter('sdss')
tau_grid = np.logspace(-2,2,50)

A = Gal_age 
def mean_age(tau):
    return  A - tau * (1.0-np.exp(-A/tau)*(1.0+A/tau)) / (1.0-np.exp(-A/tau)) 

mean_age_grid = []
for tau in tau_grid:
    mean_age_grid.append(mean_age(tau))
    
   
# Find what tau is needed to provide appropriate colors for 
# red or blue galaxy  


# Assume solar metallicity, retrieve ur photometry at different times 
sp.params['zmet'] = 20 
ur_grid=[]

for tau in tau_grid:
    sp.params['tau'] = tau
    phot = sp.get_mags(bands=sdss_filters, tage=Uni_age)
    ur_grid.append(phot[0]-phot[3])

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)  
ax1.plot(mean_age_grid, ur_grid, 'bo', mean_age_grid, ur_grid, ':k')
ax1.axhline(y=1.6, color='blue', ls='--')
ax1.axhline(y=2.3, color='red', ls='--')
ax1.set_ylabel('u-r SDSS color')
ax1.set_xlabel('Mean age [Gyrs]')  
ax1.set_title('Color spread for solar metallicity')
plt.savefig('SSP_age_vs_ur_color_galaxies.png')
plt.show()


tau_red = np.interp(2.3, ur_grid[::-1], tau_grid[::-1])
tau_blue = np.interp(1.6, ur_grid[::-1], tau_grid[::-1])

print 'For u-r = 2.3, tau is ', tau_red
print 'For u-r = 1.6, tau is ', tau_blue
 
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(tau_grid, ur_grid, 'bo', tau_grid, ur_grid, ':k')
ax2.set_xscale('log')
ax2.axhline(y=2.3, xmin=1.0, xmax=9.0, color='red', ls='--')
ax2.axhline(y=1.6, xmin=5, xmax=11, color='blue', ls='--')
ax2.axvline(x=tau_red, ymin=2.1, ymax=2.5,  color='red', ls='--')
ax2.axvline(x=tau_blue,ymin=1.4, ymax=1.8 , color='blue', ls='--')
ax2.set_ylim(ymax=2.9)
ax2.set_ylabel('u-r SDSS color')
ax2.set_xlabel(r'$\tau$')
plt.savefig('SSP_vary_tau_solar_z.png')
plt.show()