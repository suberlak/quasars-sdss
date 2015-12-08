# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 18:57:29 2015

@author: suberlak
"""

import os
import numpy as np


def get_stars_in_master_files():
    ''' A short function to read in the old star master files, 
    and check how many different stars are in 
    1) beginning 18 master files that had been used before alongside QSO
    2) all 128 master files . 
    This allows to check whether using only 18 random files out of 128 
    did not introduce any selection effects towards the available stars 
    
    Returns:
    -------
    a : list of names of stars present in the 128 master files 
    b : list of names of stars present in the first 18 master files, the
        same 18 files that were used in the analysis prior to making 
        one master per lightcurve 
        
    '''
    inDirStars   = 'sf_TRY/sf_stars/'
    inDirQSO = 'sf_TRY/sf_qso/'
    
    S = os.listdir(inDirStars)
    Q = os.listdir(inDirQSO)
    
    # Find out what stars are available for analysis in all 128 stellar master files 
    
    star_names = [] 
    
    for i in range(len(S)):
        print 'Processing ', S[i],'-', i, 'out of ', len(S)
        master = np.genfromtxt(inDirStars+S[i], dtype=str)
        names = master[:,3]
        uniq_names = np.unique(names)
        star_names.append(uniq_names)
    
    a = []
    extend = a.extend
    for s in star_names:
        extend(s)
    np.savetxt('Master_files_stars_128_names.txt', a, fmt='%s', delimiter=' ' )
    
    # Find out what stars were available for analysis in 18 master files...
    star_names_choice = []
    for i in range(len(Q)):
        print 'Processing ', S[i]
        master = np.genfromtxt(inDirStars+S[i], dtype=str)
        names = master[:,3]
        uniq_names = np.unique(names)
        star_names_choice.append(uniq_names)
    
    b = []
    extend = b.extend
    for s in star_names_choice:
        extend(s)
        
    np.savetxt('Master_files_stars_18_names.txt', b, fmt='%s', delimiter=' ')
    
    return a,b 

# read in the result of calculations 
a = np.genfromtxt('Master_files_stars_128_names.txt', dtype='str')
b = np.genfromtxt('Master_files_stars_18_names.txt', dtype='str')
# grab info about those stars from a catalog 
    
def get_stars_catalog():
    File = 'CRTS_SDSS_cross_matched_stars_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    stars_catalog = {}
    print 'zipping CRTS-SDSS stars catalog...'
    for label, column in zip(colnames, datatable.T):
        stars_catalog[label] = column
        
    return  colnames, stars_catalog

cols2 , star_cat= get_stars_catalog()
    
# plot their properties... 