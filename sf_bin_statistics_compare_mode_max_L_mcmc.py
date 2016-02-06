# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:44:24 2015

@author: suberlak
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os 
import re 


pre_obj = ['StarB', 'StarR', 'QSO']
pre_mag = ['17-19', '18.5-19']
pre_corr = ['','_corr']


for obj in pre_obj : 
    for mag in pre_mag:
        for corr in pre_corr :
            pre = obj+'_bins_'+mag+corr
            fname = pre +'_N_max_sig_full_sig.txt'
            
            
            # check if the file was not already saved...
            if os.path.exists(fname) : 
                num_lines = sum(1 for line in open(fname))
                print num_lines
                if num_lines == 200 : 
                    print fname
                    print 'This one was already compared, moving on...\n'
            # if not, go ahead with the comparison        
            else:
                print 'Reading in ', fname
                #logL_bin, N, max_sig, full_sig = compare_sigma(pre)
                