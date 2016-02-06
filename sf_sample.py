# -*- coding: utf-8 -*-
"""
Created on Thu May 28 10:55:03 2015

@author: suberlak

A program to read-in a sample of QSO, Stars red, Stars blue, 
selected from my master files according to  SDSS_r < 20 ,  
CRTS_err < 0.3  ,  and log10(tau) < 1.7 . 

Further, I cut QSO     -1.5 < delflux < 1.5 , and 
              stars    -1.0 < delflux < 1.0 

Here I am reading in those files,  and plotting  them in various ways. 

"""

import numpy as np 
import matplotlib.pyplot as plt 

# READ IN THE FILES  

def read_file(obj):
    if obj == 'StarB' : 
        File = 'Sample_4_Stars_blue_314248_lines.txt'
    if obj == 'StarR' :
        File = 'Sample_4_Stars_red_373306_lines.txt'
    if obj == 'QSO' :
        File = 'Sample_4_QSO_721283_lines.txt'
        
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    data = {}
    print 'Zipping ', obj, ' data from ', File, ' ...'
    for label, column in zip(colnames, datatable.T):
        data[label] = column

    return data
    
starB = read_file('StarB') 
starR = read_file('StarR') 
QSO = read_file('QSO') 


