PAPERS:

Challenged one: Graham et al. 2011

Basis : Ivezic et al. 2014 (Optical variability of quasars : a damped random walk)
Ivezic et al. 2014 (Optical selection of quasars : SDSS and LSST)

Method: MacLeod et al. 2012 (A description of quasar variability measured using 
repeated SDSS and POSS imaging)
MacLeod et al. 2011 (Quasar selection based on photometric variability)
MacLeod et al. 2010 ( Modeling the time variability of SDSS stripe 82 quasars as 
a damped random walk)




PLOTTING STANDARD STARS 

1) I loaded the stripe82calibStars_v2.6.dat  data from 
http://www.astro.washington.edu/users/ivezic/sdss/catalogs/stripe82.html   
using the stars1.py code, which ignores the comment lines (#), ignores the 
first column (which only has a keyword), and saves as my_array1.npy , which is 
much easier to work with 

REMEMBER:
array[row,column] counting from 0 ! 

COLUMNS:
0	RA  print arr[0:5,0]  [ 308.5002136  308.5000916  308.5000916  308.5018921  308.5046082]
1	Dec print arr[0:5,1]  [-1.227721  -1.2402642 -1.2171559 -1.1764922 -1.174189 ]
2	RArms   [ 0.032  0.015  0.009  0.005  0.028]
3	Decrms [ 0.032  0.015  0.024  0.024  0.038]
4 	Ntot  [ 4.  4.  4.  4.  4.]
5	Ar    [ 0.587  0.596  0.579  0.559  0.556]
6	u  Nobs 
7	u  mmed 
8 	u  mmu 
9 	u  msig
10	u  mrms
11	u  mchi2
12	g  Nobs 
13	g  mmed
14	g  mmu
15 	g  msig
16	g  mrms
17 	g  mchi2
18 	r  Nobs 
19	r  mmed 
20	r  mmu
21 	r  msig
22	r  mrms
23	r  mchi2
24	i  Nobs 
25	i  mmed
26 	i  mmu
27	i  msig
28 	i  mrms
29	i  mchi2
30 	z  Nobs 
31	z  mmed
32	z  mmu
33 	z  msig
34	z  mrms
35	z  mchi2


# 1) RA Dec RArms Decrms: the mean position and its rms per coordinate,
#     this is J2000, decimal degrees for RA and Dec, and arcsec for rms
#     NB: standard errors can be computed as rms/sqrt(Ntot)
#
# 2) Ntot: the total number of epochs
# 3) Ar: the Schlegel, Finkbeiner & Davis (1998) ISM extinction value in 
#    the r band; extinction in other bands can be computed as [Rv=3.1]: 
#     Am = Cm*Ar, with Cm=(1.873, 1.377, 0.758, 0.537) for m=(ugiz) 
# 4) and then in each band (ugriz, in this order):
#      (Nobs mmed mmu msig mrms mchi2), which are: 
#      the total number of observations in this band
#      the median magnitude 
#      the mean magnitude
#      the standard error for the mean (1.25 larger for the median)
#      the root-mean-square scatter
#      chi2 per degree of freedom (computed using the mean magnitude)
#


2) In the new set : plot1.py , I load the my_array1.npy  array, which consists 
of 36 columns, and 1006849 rows of data  (  check by arr.shape ). I plot  rms vs. 
magnitude for u, g,r,i ,z bands   separately, and altogether 

3) When comparing any band to Fig. 10 of Graham et al. 2011 I don't get the 
increasing slope when considering the full sample - the slope is only visible 
if I take a small subset of the SDSS data. 

--------------

I'm not sure whether I was to plot the standard stars, or something else... 
Have a look at lightcurves:

http://www.astro.washington.edu/users/ivezic/sdss/catalogs/S82variables.html

"here is a little blurb for Chris about your work, and request for help from you 
at the end. 

Chris, Chelsea and I convinced ourselves that Stripe 82 sampling was not
introducing significant bias in best-fit tau estimates after she did analysis
summarized in fig. 5 in
http://www.astro.washington.edu/users/ivezic/Publications/macleod2010.pdf
The smallest input tau she used was 100 days in observer’s frame, 
and the Graham et al. 54 days time scale in rest frame corresponds to 
152 days in observers frame. Therefore, Stripe 82 sampling ought be able to 
uncover evidence for such short time scales. It would be good to repeat this 
analysis with even shorter input tau scales (we have to run into a problem with 
stripe 82 sampling at some short scale). I would propose that first you try to
reproduce this plot as it was done by Chelsea, and then extend the tau
sampling range all the way to tau = 10 days.

Chelsea, could you please help Chris with your tools that you used to make
that figure? 


PS
OK, a few more details about the big picture: - the above will demonstrate that 
SDSS data are good enough to recover
  time scales claimed by Graham et al. (but you did *not*)

- Chris will download CRTS light curves for standard stars (which we know 
  are not variable) from  (those that are bright enough to be in CRTS)
http://www.astro.washington.edu/users/ivezic/sdss/catalogs/stripe82.html

and produce rms vs. magnitude plot using CRTS data (to show that their 
fig. 10 is bogus - that variation is not due to quasar variability but due to
their errors)

- Chris will download CRTS light curves for quasars from your Stripe 82 
sample, available from 
http://www.astro.washington.edu/users/ivezic/cmacleod/qso_dr7/Southern.html

and will determine tau and SFinfinity for them using *CRTS* data; we will 
show that the CRTS based values are very different from the SDSS-based
values (because their data are of lower quality). " 

"Chelsea:
Hello,

This sounds like a good plan (sorry for the delay in responding).

I've attached the fortran code (qso.f) which fits a group of light curves with 
format as a DRW, and compiles as

g77 qso.f -ffixed-line-length-none -o qsoexe

The reference for this code is Kozlowski et al. (2010). The light curves needs 
the format: time, magnitude, error (no header lines)

But first, the code needs a couple input files (see attached examples for a list
 of two light curves):

1. process.dat: first character (in header line) needs to be the number of light 
curves you are fitting. Then one row for each light curve, where the first column
 is the number of light curve data points. Then, columns for ra, dec, redshift, 
and abs mag, which don't matter for the fits (just use the values in the example
 file).  Then the light curve name.

2. file called 'input': just the number 1, since you are only fitting one band.

Currently, the fortran code will only accept magnitudes with values between -30 
and 30 and spaced at least
0.5 days (or whatever time units you use) apart.  You can change these settings 
on lines 101 and 104 in qso.f.

The code evaluates the goodness of fit of the DRW via the parameters chi^2DRW
and the likelihood ('Plike'). The following table (see footnotes) explains how 
to make the cuts in fit quality:
http://www.astro.washington.edu/users/cmacleod/qso/S82var_format_m
or you can also just see the paper.

The code outputs a bunch of fort.* files:

   general statistics            --> fort.39
   PRH method                    --> fort.40
   Forecasting  method           --> fort.41
   monte carlo PRH               --> fort.50
   monte carlo forecasting       --> fort.51
   bad epochs list               --> fort.80

 The ones of interest are fort.39 (containing some of the fit quality parameters)
 and fort.40 (containing DRW parameters, tau etc.):

fort.40 columns (rows are same as in process.dat): 
----
edge      ; edge flag -- 1 on edge of grid, 2 too close to edge to properly 
centroid peak
mu        ; average magnitude (best estimate from PR method)
maxlike   : max likelihood                               
minchi    ; min chi^2                                    
minchired ; min chi^2/dof                                
sigma     ; best sigma (log)                             
sig_err_lo; 1-sigma lower limit on sigma (log) (Bayesian)
sig_err_hi; 1-sigma upper limit on sigma (log) (Bayesian)
tau        ; best tau (log)                             
tau_err_lo ; 1-sigma lower limit on tau (log) (Bayesian)
tau_err_hi ; 1-sigma upper limit on tau (log) (Bayesian)
tau_med    ; median tau (log) (Bayesian)
sig_med     ; median sigma (log) (Bayesian)

and columns in fort.39:
----
index    ,$;index
z_pr     ,$;redshift
distmod  ,$;dist modulus
mi_pr    ,$;M_i
muinit   ,$;initial estimate of mean mag
chiconst ,$;chi^2 for fitting LC as constant
chilinear,$;chi^2 for fitting it as a linear trend
Pnoise   ,$;Pnoise = likelihood of noise solution 
Pinf     ,$;Pinf   = likekihood of tau->infinity  
Plike    ,$;Plike  = PR likelihood                
chinoise ,$;chi^2 of noise                        
chiinf   ,$;chi^2 of tau->infintiy                
chi_pr   ,$;chi^2 of PR result                    
signoise ,$;log(sigma) of noise                   
siginf   ,$;log(sigma) of inf                     
npts       ;# points in lightcurve


So, to produce Fig 5, I
1) generated light curves using the method in Kelly et al. 2009 (see appendix), 
for a certain range in tau and SFinf, and adopting a certain cadence and errors.  
2) fit these simulated LCs as a DRW, using this attached code, and made cuts in 
light curve quality (excluding the blue and red dots in fig 5) before comparing 
the input and output parameters.
"





