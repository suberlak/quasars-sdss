# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 12:15:58 2015

@author: astronomy

Plot   mag   vs mag_error   for both stars and galaxies 


"""


import os
import numpy as np 
import matplotlib.pyplot as plt 


# Read in the LC files 

qso_or_star = 'star'

if qso_or_star == 'star':
    inDir =  './stars_CRTS_proc_err_w_good/' 
    objType= ' stars ' 

if qso_or_star == 'qso' :
    inDir = '../QSO_CRTS_proc_err_w_good/'
    objType = ' Quasars '


outDir = './sf_TRY/' 

if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 


inFiles = os.listdir(inDir)

# Initialise the plot here ...


mag_hold = np.empty(0,dtype=float)
err_hold = np.empty(0,dtype=float)

def plot_err(mag_hold, err_hold, ttl ):
    plt.clf()
    plt.xlim(10,25) 
    H, xedges,yedges = np.histogram2d(mag_hold,err_hold,bins=50)  
    H = np.rot90(H) ; H = np.flipud(H) ;  Hmasked = np.ma.masked_where(H==0,H)
    plt.pcolormesh(xedges, yedges, Hmasked)   # as  a color map 
    plt.title('Photometric error ' + ttl); plt.xlabel('Magnitude'); plt.ylabel(r'$\sigma _{M}$')
    plt.savefig('SF_error_plot_star.png') ; plt.show()


printcounter = 0

for i in range(len(inFiles)):
    file = str(inFiles[i])
    print '\nFile ',i, 'out of',  len(inFiles)
    # load the mjd, flux, and flux_error  from the read file 
    mjd,flx4,err = np.loadtxt(inDir+'%s' % (file),usecols=(0,1,2),unpack=True)
    
    mag_hold = np.append(mag_hold, flx4)
    err_hold = np.append(err_hold, err)
    
    if  (printcounter == 1000 ):
        nqso =  '- ' + str(i) + objType 
        plot_err(mag_hold, err_hold, nqso)
        printcounter = 0
    
    printcounter += 1     
    # Add this file to  the plot ... 
    
    
plot_err(mag_hold, err_hold, ttl = '- ' + str(i) + objType )


