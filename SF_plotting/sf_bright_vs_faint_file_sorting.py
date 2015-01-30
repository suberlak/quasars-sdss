# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 15:20:43 2015

Read in QSO and stars,  and separate them into two directories for each group,
depending on their average magnitude. 



"""


import os
import numpy as np 
import matplotlib.pyplot as plt 


# Read in the LC files 

qso_or_star = 'star'

if qso_or_star == 'star':
    inDir =  './stars_CRTS_proc_err_w_good/' 
    objType= ' stars ' 
    outDirB = './stars_bright/' 
    outDirF = './stars_faint/' 

if qso_or_star == 'qso' :
    inDir = '../QSO_CRTS_proc_err_w_good/'
    objType = ' Quasars '
    outDirB = './qso_bright/' 
    outDirF = './qso_faint/' 



# make the bright and faint dirs if need be ... 
if not os.path.exists(outDirB): os.system('mkdir %s' % outDirB) 
if not os.path.exists(outDirF): os.system('mkdir %s' % outDirF) 

inFiles = os.listdir(inDir)

printcounter = 0

faint_list = []
bright_list = []

for i in range(len(inFiles)):
    file = str(inFiles[i])
    
    # load the mjd, flux, and flux_error  from the read file 
    flx4 = np.loadtxt(inDir+'%s' % (file),usecols=(1,),unpack=True)
    
    # Split into groups depending on the brightness....
    
    if  18.0 < np.median(flx4) and np.median(flx4) < 19.0   :
        bright_list.append(file)
    if  19.0 < np.median(flx4) and np.median(flx4) < 20.0   :
        faint_list.append(file)
    
    if  (printcounter == 1000 ):
        print '\nFile ',i, 'out of',  len(inFiles)
        printcounter = 0
    
    printcounter += 1     
    # Add this file to  the plot ... 
    
np.savetxt(qso_or_star + '_faint_list.txt', faint_list, delimiter ='\n', fmt="%s")  
np.savetxt(qso_or_star + '_bright_list.txt', bright_list, delimiter ='\n', fmt="%s")  



