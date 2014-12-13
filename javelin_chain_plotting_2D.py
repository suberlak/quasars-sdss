# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 17:01:19 2014

@author: suberlak

Plotting the Markov Chains from Javelin

INPUT : 

1) list of chains for a given band 
2) results of comparison of Chelsea's results to Javelin results 

OUTPUT : 

1) Plots of ln(tau)  vs ln(sigma) for each chain, together with overplotted 
crosshairs marking the position of 50% tau and sigma for javelin  (blue) and Chelsea (red).

NOTES:

The 50% hpd for Javelin was used to compare the results with Chelsea, so if it 
differs significantly from the plotted maximum of the chain distribution,
then 50% of hpd does not represent well the Javelin fit, and one needs to use 
a different value of sigma and tau to compare with Chelsea's fits.


"""
import numpy as np 
import matplotlib.pyplot as plt 
from math import  isinf

dir_choice = ['QSO_try/CRTS_chains_ALL/','QSO_SDSS_chains/','QSO_SDSS_analysis/','QSO_SDSS_chains/test/', 's82drw/','QSO_SDSS_chains/test/figs/', 'QSO_SDSS_chains/u_figs/']

dir_in = dir_choice[1]
dir_out = dir_choice[6]
band = 'u'

'''
The file with Chelsea's results matched with Javelin results, which gives values of 
sigma and tau that Chelsea got for a given quasar.

'''
chelsea_results_matched = dir_choice[4]+'javelin_SDSS_chelsea_comparison_u_band.txt'

chelsea= np.loadtxt(chelsea_results_matched,dtype=str)
qso_name = chelsea[:,0]
tau_jav=chelsea[:,5]
sig_jav=chelsea[:,7]
tau_che =chelsea[:,6]
sigma_ch = chelsea[:,8]

'''
NOTE : must make a chain_list_  ... .ls  file before running the program!
in QSO_SDSS_chains/  run :
ls ch_u_*.txt_chain.dat > chain_list_u.ls

'''
filename = dir_in + 'chain_list_'+band+'.ls'
files=np.loadtxt(filename,dtype=str)

# initialise storing vecfiles_rtors

sigma_l =  np.empty(0,dtype=float)
sigma_m =  np.empty(0,dtype=float)
sigma_h =  np.empty(0,dtype=float)
tau_l =  np.empty(0,dtype=float)
tau_m =  np.empty(0,dtype=float)
tau_h =  np.empty(0,dtype=float)
files_read = np.empty(0,dtype=str)

# load multiple chains 
for j in range(100):   #len(files)
    fchain = dir_in+files[j]
    qso_chain = files[j][5:-14]
    ind = np.where(qso_name == qso_chain)[0][0]
    sig_ch = sigma_ch[ind].astype(float)
    tau_ch = tau_che[ind].astype(float)
    t_j = tau_jav[ind].astype(float)
    s_j = sig_jav[ind].astype(float)
    flatchain= np.genfromtxt(fchain)
    fig1 = plt.figure()
    x=flatchain[:,0]
    y=flatchain[:,1]
    
    xinf = np.asarray(map(isinf,x),dtype=bool)
    yinf = np.asarray(map(isinf,y),dtype=bool)
    ttlinf = xinf + yinf
    # ttlwh = np.where(ttlinf == True)  list of good indices
    gi = -ttlinf  # good_indices 
    non_inf = len(np.where(gi == True)[0])
#    print 'Out of ', len(x),' rows, we have ', non_inf, ' of those that do not',\
#    ' have any infinities, and only those are used for plotting '
    
    plt.plot(x[gi],y[gi],'.r')
    nbins =100
    
    H, xedges,yedges = np.histogram2d(x[gi],y[gi],bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H==0,H)
    fig2 = plt.figure()
    plt.pcolormesh(xedges, yedges, Hmasked)
#    plt.xlim((-0.5,2))
#    plt.ylim((-1,1))
    title = 'Chain plot : '+band+' band for '+files[j]
    plt.title(title)
    plt.axhline(np.log(tau_ch),color='r',lw=2)
    plt.axvline(np.log(sig_ch),color='r',lw=2)
    plt.axhline(np.log(t_j),color='b',lw=2)
    plt.axvline(np.log(s_j),color='b',lw=2)
    
    plt.xlabel('ln (sigma)')
    plt.ylabel('ln (tau)')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    fname3=dir_out+'sdss_'+band+'_band-chain_plot_'+qso_chain+'.png'
    plt.savefig(fname3)

print 'We plotted ', j, 'chains, out of ', len(files), ' to sample the chain distribution'

#    flatchain_whole = np.copy(flatchain)
#    ndim = flatchain.shape[1]
#    hpd = np.zeros((3,ndim))
#    chain_len = flatchain.shape[0]
#    
#    """ 
#    The chain consists of two columns if we are fitting tau and sigma, which are natural 
#    logs of the true values - thus at the end we need to tae 
#    pct1sig are points at which we probe the chain, scaled to the length of the chain 
#    for the chain of length 5000, such values will be at 
#    positions 800,2500, 4200.  
#     """
#    pct1sig = chain_len * np.array([0.16,0.50,0.84])  
#    medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
#    vars = ["sigma", "tau"]
#    set_verbose=True  
#    
#    for i in xrange(ndim):
#        vsort = np.sort(flatchain[:,i])  # sorts the array along either sigma or tau dimension 
#        hpd[:,i] = vsort[medlowhig] # picks out values at the positions for the 
#                                    # points at 15%, 50%, and 84% of the maximum posterior distribution
#        if set_verbose :
#            print("HPD of %s"%vars[i])
#            if i < 2 :
#                # tau and sigma are stored as natural logs - other variables may not 
#                print("low: %8.3f med %8.3f hig %8.3f"%tuple(np.exp(hpd[:,i])))
#            else :
#                print("low: %8.3f med %8.3f hig %8.3f"%tuple(hpd[:,i]))
#                        
#    
#    sigma_lmh = hpd[:,0] 
#    tau_lmh = hpd[:,1]
#    
#    print 'HPD of sigma', np.exp(sigma_lmh)
#    print 'HPD of tau', np.exp(tau_lmh)
#
#    exp_sigma = np.exp(sigma_lmh)
#    exp_tau = np.exp(tau_lmh)
#    
#    sigma_l = np.append(sigma_l,exp_sigma[0])
#    sigma_m = np.append(sigma_m,exp_sigma[1])
#    sigma_h = np.append(sigma_h,exp_sigma[2])
#    tau_l = np.append(tau_l, exp_tau[0])
#    tau_m = np.append(tau_m, exp_tau[1])
#    tau_h = np.append(tau_h, exp_tau[2])
#    quasar_name = files[j][5:]           
#    print 'band', band, quasar_name
#    files_read=np.append(files_read,quasar_name)
#    
### save all the information to output file
#
#fout = dir_out + 'javelin_SDSS_chain_results_'+band+'_band_TEST.txt'
#DAT= np.column_stack((files_read, sigma_l, sigma_m, sigma_h, tau_l, tau_m, tau_h))

# sort the DAT column accoring to QSO names 
newDAT=DAT[DAT[:,0].argsort()]

np.savetxt(fout,newDAT, delimiter=" ", fmt="%s")

print 'We saved the result to file ', fout