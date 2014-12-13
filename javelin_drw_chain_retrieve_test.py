# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 17:01:19 2014

@author: suberlak

A new program that does not require Javelin to run 

Loading chains bit .

REQUIREMENTS:

Before running you need to update the chain list file (in terminal): 

> ls new7_*_chain.dat > chains_list_new8.ls  

NOTE:

The chain name has to have the form  

xxx7_123456.78+123456.7_chain.dat , 

because the ra and dec of quasar is taken from the name string
> quasar_name = files[j][5:23]

OUTPUT: 

A text file like javelin_chain_results_new8_sigma_tau.txt    ,  with columns

QSO_name           sigma_l,       sigma_m,       sigma_h,       tau_l,        tau_m,       tau_h            
222857.83-010734.8 0.147325389744 0.178587132035 0.227150026128 70.2886684856 116.70354396 209.732145982

NOTE : 

THE SCRIPT TAKES QUASAR NAME FROM THE INPUT FILE TO IDENTIFY A CHAIN . 
FOR SDSS QSO'S IT IS quasar_name = files[j][5:]    
BUT FOR CRTS IT IS DIFFERENT  !!!  


"""
import numpy as np 
import matplotlib.pyplot as plt 
from math import  isinf

dir_choice_short = ['qso_drw_chains/','qso_drw_chains/no_prior/', 'qso_drw_analysis/', 'qso_drw_analysis/figs_w_prior/', 'qso_drw_analysis/figs_no_prior/']
dir_choice_medium = ['qso_drw_medium_LC_chains/with_prior/','qso_drw_medium_LC_chains/no_prior/', 'qso_drw_medium_analysis/']

dir_in = dir_choice_medium[1]  # with prior : 0,3,2  |   no prior :  1,4,2
#dir_figs = dir_choice[4]
dir_out = dir_choice_medium[2] 

fout = dir_out + 'javelin_drw_test_chain_results_no_prior_all.txt'


'''
NOTE : must make a chain_list_  ... .ls  file before running the program!

'''
filename = dir_in + 'chains_err_all.ls'
files=np.loadtxt(filename,dtype=str)

# storing arrays 
files_read = np.empty(0,dtype=str)
sigma_max =  np.empty(0,dtype=float)
tau_max   =  np.empty(0,dtype=float)

sigma_l =  np.empty(0,dtype=float)
sigma_m =  np.empty(0,dtype=float)
sigma_h =  np.empty(0,dtype=float)
tau_l =  np.empty(0,dtype=float)
tau_m =  np.empty(0,dtype=float)
tau_h =  np.empty(0,dtype=float)

error_lines =[]

# load multiple chains 
for j in range(len(files)):   # len(files)
    fchain = dir_in+files[j]
    flatchain= np.genfromtxt(fchain)
    
    if files[j][-14:] =='.dat_chain.dat' : drw_name = files[j][3:-14]
    else  :
        if  files[j][-9:] =='chain.dat'  : drw_name = files[j][3:-10]  
        else : 
            print 'PROBLEM !!' 
            error_lines.append[j]
          
    print  '\n',drw_name
    
    # finding the maximum of the 2D  distribution  
#    x=flatchain[:,0]   # ln sigma 
#    y=flatchain[:,1]   # ln tau 
    sigma = np.exp(flatchain[:,0]) 
    tau = np.exp(flatchain[:,1])

     # I make x and y  log10()  because I want to plot in log10()
    
    x=np.log10(sigma)  # log10(sigma)
    y=np.log10(tau)    # log10(tau)
    
    xinf = np.asarray(map(isinf,x),dtype=bool)
    yinf = np.asarray(map(isinf,y),dtype=bool)
    ttlinf = xinf + yinf
    # ttlwh = np.where(ttlinf == True)  list of good indices
    gi = -ttlinf  # good_indices 
    non_inf = len(np.where(gi == True)[0])
    
    # INITIALIZE THE HISTOGRAM  HERE, to find max values ...   
    nbins =100
    H, xedges,yedges = np.histogram2d(x[gi],y[gi],bins=nbins)
    a,b = np.where(H == H.max())    
    x_max = xedges[a[0]]
    y_max = yedges[b[0]]
    
    sigma_max = np.append(sigma_max,np.power(10.0,x_max))
    tau_max = np.append(tau_max,np.power(10.0,y_max))
    
    
    # the Javelin way : finding positions .... 
    ndim = flatchain.shape[1]
    hpd = np.zeros((3,ndim))
    chain_len = flatchain.shape[0]
    
    
    """ 
    The chain consists of two columns if we are fitting tau and sigma, which are natural 
    logs of the true values - thus at the end we need to tae 
    pct1sig are points at which we probe the chain, scaled to the length of the chain 
    for the chain of length 5000, such values will be at 
    positions 800,2500, 4200.  
     """
    pct1sig = chain_len * np.array([0.16,0.50,0.84])  
    medlowhig  =pct1sig.astype(np.int32) # expresses the pointers above as integers
    vars = ["sigma", "tau"]
    set_verbose=True  
    
    for i in xrange(ndim):
        vsort = np.sort(flatchain[:,i])  # sorts the array along either sigma or tau dimension 
        hpd[:,i] = vsort[medlowhig] # picks out values at the positions for the 
                                    # points at 15%, 50%, and 84% of the maximum posterior distribution 
    ln_sigma_lmh = hpd[:,0] 
    ln_tau_lmh = hpd[:,1]
    
    print 'HPD of sigma', np.exp(ln_sigma_lmh)
    print 'HPD of tau', np.exp(ln_tau_lmh)

       
    
    sigma_lmh = np.exp(ln_sigma_lmh)
    tau_lmh = np.exp(ln_tau_lmh)
   
    log_sigma_lmh = np.log10(sigma_lmh)
    log_tau_lmh = np.log10(tau_lmh)
   
    # values below are just true values 
   
    sigma_l = np.append(sigma_l,sigma_lmh[0])
    sigma_m = np.append(sigma_m,sigma_lmh[1])
    sigma_h = np.append(sigma_h,sigma_lmh[2])
    tau_l = np.append(tau_l, tau_lmh[0])
    tau_m = np.append(tau_m, tau_lmh[1])
    tau_h = np.append(tau_h, tau_lmh[2])
 
    
    files_read=np.append(files_read,drw_name)
    
    
#    plot on the chain representation all the values found....
    
#    H = np.rot90(H)
#    H = np.flipud(H)
#    Hmasked = np.ma.masked_where(H==0,H)
#    fig2 = plt.figure()
#    plt.pcolormesh(xedges, yedges, Hmasked)
#    
#    title = 'Chain plot for '+files[j][3:-14]
#    plt.title(title)
#    plt.axhline(y_max,color='b',lw=3)
#    plt.axvline(x_max,color='b',lw=3)
#    
#    colors = ['orange','red','purple']
#    for k in range(3):
#        plt.axvline(log_sigma_lmh[k],color=colors[k],lw=1)   # low 
#        plt.axhline(log_tau_lmh[k],color=colors[k],lw=1)
#   
#    plt.xlabel(r'$\log_{10}{\,\sigma}$',fontsize=16)
#    plt.ylabel(r'$\log_{10}{\,\tau}$',fontsize=16)
#    cbar = plt.colorbar()
#    cbar.ax.set_ylabel('Counts')
#    fname3=dir_figs+'chain_2D_'+files[j][3:-14]+'.png'
#    plt.savefig(fname3)
#    print 'Saving  ', fname3
#    
    
## save all the information to output file

DAT= np.column_stack((files_read, sigma_max,  tau_max, sigma_l, sigma_m,sigma_h, tau_l,tau_m,tau_h ))

# sort the DAT column accoring to QSO names 
#newDAT=DAT[DAT[:,0].argsort()]

np.savetxt(fout,DAT, delimiter=" ", fmt="%s")
print 'We retrieved the  javelin fitting results for ', len(files_read) ,'chains'
print 'We saved the result to file ', fout

if error_lines == '' : print 'All lines fulfilled name-extracting ways' 