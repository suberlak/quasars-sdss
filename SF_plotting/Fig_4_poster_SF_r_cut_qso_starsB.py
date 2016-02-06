import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from scipy.optimize import curve_fit
'''
A program to make the Fig.4 , which is the 
comparison of three magnitude bins for SF

It needs results of Fig_4_SF_three_panels.ipynb,
where I calculate SF and error for three magnitude bins, applying 
appropriate correction factors. 

'''
def model_sf(t, sf_inf=0.25, tau = 1.0):
    br = 1.0-np.exp(-t/tau)
    sf = sf_inf * np.power(br,0.5)
    return sf
        

rcParams['ytick.labelsize'] = 25
rcParams['xtick.labelsize'] = 25
rcParams['axes.labelsize'] = 35
rcParams['axes.linewidth'] = 3
rcParams['font.size'] = 25
rcParams.update({'figure.autolayout': False})
rcParams['legend.numpoints'] = 1

Min_arr = [17, 18, 18.5]
Max_arr = [18, 18.5, 19]
obj_arr = ['qso', 'starsB']
fc_arr =[0.864, 1.091, 1.298]    # Correction factors based on blue stars 
labels_arr = ['Quasars', 'Blue Stars']
colors = ['black','blue']


b = 'r_cut'
dir = 'crucial_results/Feb_2/'
fig,ax = plt.subplots(3,1, figsize=(12,12), sharex=True)
fig.subplots_adjust(hspace=0)
for i in range(len(Min_arr)):
    for j in range(len(obj_arr)):
        obj = obj_arr[j]

        fname =  b+'_'+str(Min_arr[i])+'-'+str(Max_arr[i])+'_'+obj+'_fc-'+str(fc_arr[i])+'_mean_tau_sig_approx_err.txt'
        print fname
        d =  np.loadtxt(dir+fname, dtype=float)
        mean_tau = d[:,0]
        sig_approx = d[:,1]
        sig_err = d[:,2]
        ax[i].scatter(np.log10(mean_tau), sig_approx, alpha=0.3, label=labels_arr[j], color=colors[j])
        ax[i].errorbar(np.log10(mean_tau), sig_approx,sig_err, linestyle='None',color=colors[j], alpha=0.3)#
    
 
        if obj == 'qso' : 
            # Calculate the model DRW fit for QSO
        
            xdata = mean_tau
            sf = sig_approx
            popt, pcov = curve_fit(model_sf, xdata, sf)
            y = model_sf(xdata, sf_inf=popt[0], tau = popt[1]) # tau 1 year in days 

            # Fold-in the error to the model SF , plot 
            # both folded and not-folded version \
            err_sig = sig_err
            y_fold = np.sqrt((y ** 2.0)+ (err_sig ** 2.0) )
            ax[i].plot(np.log10(xdata), y_fold , lw=3, c = 'green', ls='--')
            ax[i].plot(np.log10(xdata), y , lw=3, c = 'orange', ls='--')

            # text = r'$ \mathrm{Model:}\ \tau=%.3f \,\mathrm{days} , \ SF_{\infty}=%.3f \,\mathrm{mag}$'%(popt[1],popt[0])
            # ax[i].text(x=0.75, y=0.3,s = text )

    ax[i].text(1.0, 0.25, 'mag: '+str(Min_arr[i])+'-'+str(Max_arr[i]) )
    ax[i].text(1.0, 0.15, 'fc: '+str(fc_arr[i]))
    ax[i].set_ylabel(r'$SF$')
    ax[i].set_ylim(-0.05,0.45)
    ax[i].set_xlim(0.5, 3.7)
    ax[i].grid() 
    ax[i].hlines(y=0.05, xmin =0.5, xmax=1.7, color='green', lw = 4, linestyle = '--' , alpha = 0.8  )
    # styles http://matplotlib.org/examples/lines_bars_and_markers/line_styles_reference.html 
    #axs[i].grid(axis='x')

    ax[i].set_yticks([0,0.1,0.2,0.3,0.4])
    ax[i].set_yticklabels(['0.0','0.1', '0.2', '0.3', '0.4'])
    
axbox = ax[0].get_position()
x_value=0.3
y_value=0.15
legend = ax[0].legend(fancybox=True,loc=(axbox.x0 + x_value, axbox.y0 - y_value))

#legend.get_frame().set_edgecolor('1.0')
legend.get_frame().set_alpha(0.5)
plt.rc('legend',**{'fontsize':6})
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')

ax[-1].set_xlabel(r'$log_{10} (\Delta {t})$ [days]') 
    
#ax[-1].set_xlabel(r'$\log_{10}{\tau}$', fontsize=20)

plt.savefig('Fig_4_SF_QSO_starsB_r_cut_fc-'+str(fc_arr[0])+'-'+str(fc_arr[1])+'-'+str(fc_arr[2])+'_NEW.png')
