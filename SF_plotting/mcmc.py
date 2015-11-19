import numpy as np
import emcee
import sys 

#### varying input #####
## edit this part
args = sys.argv
bin_name = args[1]
dir_input=args[2]  # dir with bins 
dir_output=args[3]

File = dir_input+bin_name

nproc = 2 # number of threads (on how many processors do you want this to run)

step_check = 25 # print out chain progress every .. steps
#########################


datatable = np.genfromtxt(File)
xi_orig = datatable[:,0]  # data points
ei_orig = datatable[:,1]  # error points 

def gaussgauss_logL(params):

	mu, sigma = params 
	if sigma < 0:  # assuming sigma cannot be negative!! 
		return -np.inf
	else:
		#ndim = len(np.broadcast(sigma, mu).shape)
 
		xi = xi_orig#.reshape(xi_orig.shape + tuple(ndim * [1]))
		ei = ei_orig#.reshape(ei_orig.shape + tuple(ndim * [1]))
		s2_e2 = sigma ** 2 + ei ** 2
		return -0.5 * np.sum(np.log(s2_e2) + (xi - mu) ** 2 / s2_e2)#,  -1 - ndim)



n_dim = 2 # dimension of parameter space, meaning we are sampling for 2 parameters
nwalkers = 2000

par_in = np.array([0., 0.3]) # where about to start the walkers, not that important
sigg_in = np.array([5., 1.]) # spread in initial position

pp0 = np.transpose(np.array([np.random.normal(par, sigg, size=nwalkers) for par, sigg in zip(par_in, sigg_in)]).reshape((n_dim, nwalkers)))
pp0[:,1] = np.fabs(pp0[:,1]) # so that the sigma is not negative

# create MC sampler

sampler = emcee.EnsembleSampler(nwalkers, n_dim, gaussgauss_logL, threads=nproc)


# steps of the chain, that gives n_it*nwalkers of actual MC steps
# checked that at least for the bins provided this is a reasonable length
n_it = 500 

chainfile = dir_output+bin_name+'_chain_25.dat' #output chain file name
print 'Running for file ', File
f = open(chainfile, "w")
f.close()
count = 0

# this runs the sampler and writes down the chain in realtime
# you can access the chain at any time without stopping the sampling 
# by accessing the file directly
for result in sampler.sample(pp0, iterations=n_it, storechain=False):
    position = result[0]
    f = open(chainfile, "a")
    for k in range(position.shape[0]):
            for m in range(position[k].shape[0]):
                f.write(str(position[k][m])+"\t")
            f.write("\n")
    f.close()
    if count%step_check == 0: 
    	print count

    count +=1

# reset the sampler so that the next one would not use the same ... 
sampler.reset() 



