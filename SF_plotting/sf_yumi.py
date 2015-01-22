#execfile('r.py')

import scipy.optimize.minpack as sc
from FileIO import * 
import os, sys, numpy, scipy
from pylab import *

#objType = sys.argv[1]

# List of LC files
inDir =  './stars_CRTS_proc_err_w_good_TRY/' 
outDir = '/sf/' 
if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 

inFiles = os.listdir(inDir)

outf = open(outDir+'sf_fitted_param_exp.txt','w')
outf.write('#ID RA DEC SF_inf SF_inf_err Tau(days) Tau_err Chi2 N\n')
for i in range(len(inFiles)):
#for i in range(1):
    file = str(inFiles[i])

 #   mjd,flx4,err = numpy.loadtxt(inDir+'out_211838.67+004035.5.txt',usecols=(0,1,2),unpack=True)
    mjd,flx4,err = numpy.loadtxt(inDir+'%s' % (file),usecols=(0,1,2),unpack=True)
    add = numpy.where(numpy.isnan(flx4) == True)[0]

    # Sort by mjd
    mjd = numpy.delete(mjd,add); flx4 = numpy.delete(flx4,add); err = numpy.delete(err,add)
    ind = numpy.argsort(mjd)
    mjd = mjd[ind]; flx4 = flx4[ind]; err = err[ind]
    if len(mjd) < 3: continue  # There should be at least 3 meaningful detections to get more than just 1 SF point
    
    # Calculate tau, delmag2 and delmagerr for each pair (k,j), where tau=t_j-t_k (j>k) 
    delflx2  = []; delflxerr = []; tau = []
    for j in range(len(mjd)-1):
        for k in range(j+1):     
            tau.append(mjd[j+1]-mjd[k])  # j from 1 and k<j
            noise2 = err[k]**2+err[j+1]**2  # delflx_(k)^2 + delflx_(j+1)^2
            delflx2.append((flx4[k]-flx4[j+1])**2)  # (flx4_(k) - flx4_(j+1))^2
            delflxerr.append((noise2)**0.5)  # sqrt(delflx_(k)^2 + delflx_(j+1)^2)
    tau = numpy.array(tau); delflx2 = numpy.array(delflx2); delflxerr = numpy.array(delflxerr) 

    int0 = numpy.argsort(tau)
    tau = tau[int0]; delflx2 = delflx2[int0]; delflxerr = delflxerr[int0]


    # Fitting function; exponential model; sf = sf_inf*(1-e^(-t/Tau))^0.5 
    fp = lambda v, x: v[0]*numpy.power((1.-numpy.exp(-numpy.abs(x)/v[1])),0.5)

    # Error function
    e  = lambda v, x, y, dy: (fp(v,x)-y)/dy  # error function for the leastsq 

    # Initial guess
    v0 = [0.2, 100.]    # initial guess for exp

    # leastsq
    v_whole, cov, info, mesg, ier = sc.leastsq(e, v0, args=(tau,numpy.sqrt(delflx2),delflxerr),full_output=True) 
    if ier != 1: print 'Fail in fitting SF for ' + str(inFiles[i])
    if cov == None:
        t0 = 0.0; t1 = 0.0
    else:
        t0 = numpy.sqrt(cov[0][0]); t1 = numpy.sqrt(cov[1][1])

    expected = v_whole[0]*numpy.power((1.-numpy.exp(-numpy.abs(tau)/v_whole[1])),0.5) # Expected Value
    chi2 = numpy.sum((numpy.sqrt(delflx2)-expected)**2/expected)

    WriteColumns(outDir+'tau%s.dat' % (file),((reshape(tau,(len(tau),1))),), header='# SF')
    WriteColumns(outDir+'sf%s.dat' % (file),((reshape(numpy.sqrt(delflx2),(len(tau),1))),), header='# SF')  
    outf.write('%.5f %.5f %.5f %.5f %.5f %i\n' % (v_whole[0],t0,v_whole[1],t1,chi2,len(tau))) 

outf.close()

#
##"""
#figure(1,facecolor='w')
#plot(tau,numpy.sqrt(delflx2),'*k')
##errorbar(tau,numpy.sqrt(delflx2),yerr=delflxerr,fmt='-')
#plot(tau,fp(v_whole,tau),'g-',label='Whole')
##legend()
#title('ID = %s' % (file))
#xlabel('tau')
#ylabel('SF')
#xlim(0.1,4000)
#show()
##"""
#
#import pdb
#pdb.set_trace()
