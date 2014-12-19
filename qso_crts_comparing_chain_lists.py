# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 20:23:37 2014

@author: astronomy
"""

import numpy as  np

crts_in =  ['QSO_CRTS_processed_err_rms/', 'QSO_CRTS_processed_err_w/']

crts_chains = ['QSO_CRTS_err_rms_chains/', 'QSO_CRTS_err_w_chains/'  ]

ch = 0

files_in = np.loadtxt(crts_in[ch] + 'out.list', dtype=str)

chains_out = np.loadtxt(crts_chains[ch] + 'chain.list', dtype=str) 