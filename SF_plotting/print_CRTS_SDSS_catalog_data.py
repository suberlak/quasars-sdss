# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:53:10 2015

@author: suberlak
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
#mp.rcParams.update(mp.rcParamsDefault)
#mp.rcdefaults()

def get_qso_catalog(catalog):
    if catalog == 's82drw':
        File = 'CRTS_SDSS_cross_matched_qso_s82drw_catalog.txt'
    if catalog == 'DB_QSO':
        File = 'CRTS_SDSS_cross_matched_qso_DB_QSO_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    qso_catalog = {}
    print 'Zipping CRTS-SDSS quasars catalog from ', File, ' ...'
    for label, column in zip(colnames, datatable.T):
        qso_catalog[label] = column
    
    qso_names = np.genfromtxt('CRTS_SDSS_cross_matched_qso_names.txt', dtype=str)    
    for i in range(len(qso_names)):
        qso_names[i] = qso_names[i][4:-4]
    print 'Read in ', len(qso_catalog['redshift']), ', quasars from CRTS'
    return  colnames, qso_catalog, qso_names
    
    
colnames, qso_catalog, qso_names = get_qso_catalog('DB_QSO')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel('N')
nbins=50
    
mask = qso_catalog['r'] != 0

hist1, bin_edges = np.histogram(qso_catalog['CRTS_avg_m'][mask], bins = nbins, density=True)
bin_cen1 = (bin_edges[:-1] + bin_edges[1:])/2

# exclude those 17 objects that don't have SDSS r mag 
#hist2, bin_edges = np.histogram(qso_catalog['r'][qso_catalog['r'] != 0], bins=nbins, density=True)

hist2, bin_edges = np.histogram(qso_catalog['r'][mask], bins=nbins, density=True)
bin_cen2 = (bin_edges[:-1] + bin_edges[1:])/2

res = qso_catalog['CRTS_avg_m'][mask]-qso_catalog['r'][mask]

ax.plot(bin_cen1, hist1, color = 'red', ls='steps', label='CRTS avg LC  mag')
ax.plot(bin_cen2, hist2, color = 'blue', ls='steps', label='SDSS r mag')
ax.set_xlim(xmin=15, xmax=22)
ax.legend(loc='upper left')
plt.show()
