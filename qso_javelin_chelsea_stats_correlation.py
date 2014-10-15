# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 18:04:33 2014

@author: suberlak

A program to take the lightcurve stats from  qso_lc_stats.py ,  
and the output of comparison of javelin and Chelsea's result for each quasar. 
It pulls from the stats table all the relevant stats for the subset of 
quasars on which javelin was run, to plot the log(tau_difference) vs log(sigma_difference),
and colour - code the third dimension according to the value of the investigated
stat parameter for a given quasar,  in a quest of looking for a correlation 
that may show to what types of quasar lightcurves javelin is failing, and if
it is  failing systematically. 

"""

