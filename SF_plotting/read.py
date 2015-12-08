import numpy as np
from numpy import percentile as pc
import matplotlib.pyplot as plt
import os
#import corner
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('text', usetex=True)
import re


#### varying input #####
## edit this part

#filename of the chain you want to analyse
pre= 'StarR_bins_18.5-19' # 'QSO_bins_18.5-19' # 'StarR_bins_18.5-19' # 'QSO_bins_18.5-19_corr'  # 'StarB_bins_18.5-19' # 'QSO_bins_18.5-19'
      # 'StarR_bins_17-19' #'QSO_bins_17-19' #'StarB_bins_17-19_corr' # 'QSO_bins_17-19_corr'  
      # 'StarB_bins_17-19' #
print 'Reading ...' , pre 
# filename of bin files (i.e. beginning of chain files)
# so that the bin number could be read from the fname

start = '_bin_'  # what precedes the bin number
end = '_xi_ei' # what is after bin number 



chain_dir = pre+'_chains/'
chains = os.listdir(chain_dir)


def get_chain_info(chainfile):
    
    
    #plot = True # if you want the corner plot
    #plot_save = True # if you want to save the corner plot
    #plot_filename = chainfile[:-4]+'.pdf'
    
    # prints the chain variation - way to check if we converged
    # will cause errors if steps<100 (if your chain did not get to 100 yet)
    print_chain_check = False 
    
    # prints values with 1sigma confidence levels
    print_errors = False
    
    # the default burn_in, if you see that posteriors have
    # strange beheavioir (are not 2D gaussian, see the plot), try and increase it
    # but if you go past~250/300 you should worry about chain length
    burn_in = 150 
    
    #########################
    
    
    #chain_dic = {'sig': 1, 'mu': 0}
    arr = np.loadtxt(chainfile)
    nwalkers = 2000
    #ndim = 2
    #steps = len(arr)/nwalkers
    chain = arr[burn_in*nwalkers:]
    
    if print_errors:
    	print "results with 1sigma confidence level:"
    	print "   mu:" , '%.4g' % pc(chain[:,0], 50),', '+str('%.4g' % pc(chain[:,0], 16))+'-'+str('%.4g' % pc(chain[:,0], 84))
    	print "sigma:" , '%.4g' % pc(chain[:,1], 50), '%.4g' % pc(chain[:,1], 16),'-', '%.4g' % pc(chain[:,1], 84)
    
    #else:
    #	print "mu, sigma:", '%.4g' %np.mean(chain[:,0]), '%.4g' %np.mean(chain[:,1])
    
    
    if print_chain_check:
    	print "if values below vary past 3 sig. digits, worry about the length"
    	burn = [100, 150, 200, 250, 300]
    	for burnin in burn:
    
    		chain_temp = arr[nwalkers*burnin:] # 
    		print burnin, np.mean(chain_temp[:,0]), np.mean(chain_temp[:,1])
    
    mu_chain = pc(chain[:,0], 50)
    sig_chain = pc(chain[:,1], 50)
 
    return mu_chain, sig_chain
    #if plot:#
    #	corner.corner(chain, labels=['mu','sigma'], quantiles=[0.16,0.5,0.84],
    #	                title_args={"fontsize": 12},
    #	               plot_datapoints=True, fill_contours=True, color='green')
    
    #if plot_save and plot:
    #	plt.savefig(plot_filename)
    #plt.show()
    

bin_N_chains = np.array([float(re.search('%s(.*)%s' % (start,end),chains[i]).group(1)) for i in range(len(chains))])



chainfile = pre+'_chains_res.txt'
if os.path.exists(chainfile):
    flag = 1
    print 'Using an existing file...'
    saved = np.genfromtxt(chainfile)
    N_saved = saved[:,0]
    
    N = list(N_saved)
    mu = list(saved[:,1])
    sig = list(saved[:,2])
    
    # those Ns that are in the file, but not yet read 
    N_to_read = bin_N_chains[~np.in1d(bin_N_chains, N_saved)]
 
else: # if such file does not exist need to initiate storage arrays 
    flag = 0
    N = []
    mu = []
    sig = []
    
 
for i in range(len(chains)):
    result = re.search('%s(.*)%s' % (start,end),chains[i]).group(1)
    N_bin = float(result)
   
    # only read values if has  not yet
    # been read 
    if flag == 1 and np.any(N_saved == N_bin):
        print 'chain #', i+1, ' was already read'
    if (flag == 1 and  np.any(N_to_read == N_bin)) or flag == 0 : 
        print 'reading chain #', i
        N.append(N_bin)
        
        mu_bin, sig_bin = get_chain_info(chain_dir+chains[i])
        mu.append(mu_bin)
        sig.append(sig_bin)
    
d = np.column_stack((N,mu,sig))    
np.savetxt(pre+'_chains_res.txt', d, delimiter = ' ', fmt='%s')
print 'Saved all chain N, sigma to ', pre+'_chains_res.txt'