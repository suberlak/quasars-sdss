# -*- coding: utf-8 -*-
"""
Created on Wed Dec 24 01:28:47 2014

@author: astronomy

Plot the stats : output of qso_CRTS_stats_match_to_javelin.py

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from math import  isinf

ch = 0
dir_in_out = ['QSO_CRTS_analysis/', 'stars_CRTS_analysis/']
datafile = dir_in_out[ch]+'javelin_sigma_tau_plus_stats_matched_good_err.txt'

data = np.loadtxt(datafile, dtype='str')

fname = data[:,0]
lc_names = data[:,1]
sigma_m = data[:,2].astype(float)
tau_m = data[:,3].astype(float)
sigma_hat = data[:,4].astype(float)
mjd_span = data[:,5].astype(float)
mag_rms = data[:,6].astype(float)
mag_mean = data[:,7].astype(float)
err_mean = data[:,8].astype(float)
N_lines = data[:,9].astype(float)

ls = np.log10(sigma_m)

def histogram1D(x,title,xlabel, name):
    plt.clf()
    plt.hist(x, range())
    font = 20 
    
    xmin = np.percentile(x,5)
    xmax = np.percentile(x,95)
    
    fig1 = plt.figure(figsize=[10,8])
    gs = GridSpec(100,100,bottom=0.18,left=0.18,right=0.88)
    ax1 = fig1.add_subplot(gs[:,:90])   
    
    ax1.hist(x,100,range=[xmin,xmax],facecolor='gray',align='mid',alpha=0.5)
    
    plt.title(title,fontsize=font)
    plt.xlabel(xlabel,fontsize=font)
    plt.ylabel('Number of counts',fontsize=font )
    
    ax1.tick_params(axis='x', labelsize=font)
    ax1.tick_params(axis='y', labelsize=font)
    fname = dir_in_out+'CRTS_stars_hist_'+name+'.png'
    plt.savefig(fname)
    print 'We saved ', fname 
    
histogram1D(err_mean[ls<0],title='Mean error, $log(sigma)<0$',xlabel='Mean error [mag]',name='err_sm' )    
histogram1D(err_mean[ls>0],title='Mean error, $log(sigma)>0$',xlabel='Mean error [mag]',name='err_lg' ) 

histogram1D(N_lines[ls<0],title='Light curve length', xlabel = 'Length', name='N_lines_sm')
histogram1D(N_lines[ls>0],title='Light curve length', xlabel = 'Length', name='N_lines_lg')

histogram1D(mag_rms)

histogram1D(mag_mean)