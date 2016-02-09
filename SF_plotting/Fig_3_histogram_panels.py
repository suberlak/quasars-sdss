import matplotlib.pyplot as plt
import numpy as np
import matplotlib 
import os
import datetime

matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=15) 


#  this is the old poster_hist_r_cut_qso_starsB_mag_tau_grid.py  program
#  all updated and developed from Fig_3_histogram_panels.ipynb 

b = 'r_cut'
outDir = os.path.join(os.getcwd(), 'Fig_3_data', datetime.datetime.now().strftime('%Y-%m-%d')+ '_histogram_table/')

def plot2Chistograms(chiQSO, chiSTAR, Xmin, Xmax, Ymin, Ymax, Xlabel, Ylabel, ax, bins=20, title=''):
    limits = [(Xmin, Xmax, Ymin, Ymax)]
    labels = [Xlabel, Ylabel]
    ax.set_xlim(Xmin, Xmax)
    ax.set_ylim(Ymin, Ymax)
    ax.set_xlabel(Xlabel, fontsize=12)
    ax.set_ylabel(Ylabel, fontsize=12)
    #plt.tick_params(axis='both', which='major', labelsize=15)
    xTitle = Xmin + 0.05*(Xmax-Xmin)
    yTitle = Ymax - 0.2*(Ymax-Ymin)
    ax.text(xTitle, yTitle, title, fontsize=12)

    # plot a histogram
    ax.hist(chiSTAR, bins=bins, normed=True, facecolor='blue', histtype='stepfilled', alpha=0.4)
    ax.hist(chiQSO, bins=bins, normed=True, facecolor='red', histtype='stepfilled', alpha=0.2)

    # plot the robust width of both distributions
    stdev_rob_QSO = 0.7414 *(np.percentile(chiQSO,75) - np.percentile(chiQSO,25) )
    stdev_rob_S = 0.7414 *(np.percentile(chiSTAR,75) - np.percentile(chiSTAR,25) )
    
   # add the titles of plots...  
    xTitle = Xmin + 0.65*(Xmax-Xmin)
    yTitle = Ymax - 0.2*(Ymax-Ymin)
    StarSigmaG = r'$'+str(stdev_rob_S)[:4]+'$'
    ax.text(xTitle, yTitle, StarSigmaG, fontsize=12)
    
    
    xTitle = Xmin + 0.65*(Xmax-Xmin)
    yTitle = Ymax - 0.35*(Ymax-Ymin)
    QSOSigmaG = r'$'+str(stdev_rob_QSO)[:4]+'$'
    ax.text(xTitle, yTitle, QSOSigmaG, fontsize=12)
    
    
Min_arr = [17, 18, 18.5]
Max_arr = [18, 18.5, 19]
tau_min_arr = [0,   2.3, 2.8, 3.2]
tau_max_arr = [1.7, 2.5, 3.0, 3.4]
xlims_arr = [5,10,10,10]

# Initialise figure 
fig, axs = plt.subplots(4,3, figsize=(8, 8))
fig.subplots_adjust(wspace=0.46, hspace=0.36, left=0.12, right=0.94, bottom=0.05, top=0.95)



for i in range(len(tau_max_arr)):  # 
    for j in range(len(Min_arr) ):  # 
       

        datafileS = outDir+b+'_'+str(Min_arr[j])+'-'+str(Max_arr[j])+'_'+'starsB'+'_mi_tau_ei-log_tau_'+\
                        str(tau_min_arr[i])+'-'+str(tau_max_arr[i])+'.txt'
        vS = np.loadtxt(datafileS, unpack=True)
        chiS = vS[0]/vS[2]
        chiSok = chiS[np.abs(chiS)<5]
        
        datafileQ = outDir+b+'_'+str(Min_arr[j])+'-'+str(Max_arr[j])+'_'+'qso'+'_mi_tau_ei-log_tau_'+\
                str(tau_min_arr[i])+'-'+str(tau_max_arr[i])+'.txt'
        vQ = np.loadtxt(datafileQ, unpack=True)
        chiQ = vQ[0]/vQ[2]
        chiQok = chiQ[np.abs(chiQ)<xlims_arr[i]]

        # plot histograms
        Xlabel = '$\chi = \Delta mag / error$'
        Ylabel = '$n / (N\Delta_{bin})$'
        Xmin = -xlims_arr[i]
        Xmax = xlims_arr[i]
        bins = 40 
        title= r'$ '+' '+ str(tau_min_arr[i])+'-'+str(tau_max_arr[i])+'$'
        plot2Chistograms(chiQok, chiSok, Xmin=Xmin, Xmax=Xmax, Ymin=0.0, 
                             Ymax=0.55, Xlabel=Xlabel, Ylabel=Ylabel, ax=axs[i,j],bins=bins,  
                         title=title)

#name = 'poster_r_cut_qso_starsB_histogram_grid.png'
name = 'Fig_3_histogram_panels.png'
if (name is None):
    plt.show() 
else:
    print 'saving plot to:', name
    plt.savefig(name, bbox_inches='tight')
