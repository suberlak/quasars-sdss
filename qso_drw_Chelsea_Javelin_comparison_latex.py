# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 14:08:37 2014

@author: astronomy

Read in results of analysis of Chelsea and Javelin files,  

and produce a nice latex table : I make a separate program for this part because

obtaining the stats takes about 20 sec of computational time 
(including making the 2D histograms, etc)

"""

import numpy as np 
from matrix2latex import matrix2latex
        
        
results_ch_arr = np.loadtxt('drw_chelsea_results_stats.txt', dtype = 'str')
print 'Chelsea results : S (1,2,3), M (1,2,3) , L(1,2,3)', results_ch_arr

col0 = (r'$\log(\tau_{fit})_{med}- \log(\tau_{0})$',r'$\log(\tau_{fit})_{rms}$' ,\
 r'$\log(\sigma_{fit})_{med} -\log(\sigma_{0})$',  r'$\log(\sigma_{fit})_{rms}$',\
 r'$\log(\hat{\sigma}_{fit})_{med} -\log(\hat{\sigma}_{0})$',  r'$\log(\hat{\sigma}_{fit})_{rms}$',\
 r'$\log(K_{fit})_{med} -\log(K_{0})$',  r'$\log(K_{fit})_{rms}$')
 
#
A = []
A.append(col0)
for line in results_ch_arr:
    A.append(line)
A_arr = np.asarray(A)
A_arrT = A_arr.transpose()
  
hr = [['','Short', 'Short', 'Short', 'Medium','Medium','Medium', 'Long','Long','Long'],\
 ['','err1','err2','err3','err1','err2','err3','err1','err2','err3']]
fc = ['$%g$', '$%.3f$', '$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$']
t = matrix2latex(A_arrT, headerRow=hr, caption='Chelsea results, $x_{0}$ is the true value, rms = $(\log(x_{fit,75\%}) - \log(x_{fit,25\%} ) 0.7413$ ',formatColumn=fc)
print t

results_jav_arr = np.loadtxt('drw_javelin_results_stats.txt', dtype = 'str')
results_jav_arrT = results_jav_arr.transpose()
print 'Javelin results : S(wp 1,2,3) (np 1,2,3) , M (wp 1,2) (np 1,2) '

A = []
A.append(col0)

res_jav = results_jav_arrT[2:]
res_javT = res_jav.transpose()

for line in res_javT:
    A.append(line)
A_arr = np.asarray(A)
A_arrT = A_arr.transpose()

prior_row = results_jav_arrT[0]
header_row2 = []
header_row2.append(' ')
for word in prior_row:
    header_row2.append(word)

error_row = results_jav_arrT[1]
header_row3 = []
header_row3.append(' ')
for word in error_row:
    header_row3.append(word)
    
hr = [[' ', 'Short','Short','Short','Short','Short','Short','Medium','Medium','Medium','Medium'], header_row2, header_row3]
fc = ['$%g$', '$%.3f$', '$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$','$%.3f$']
u = matrix2latex(A_arrT, headerRow=hr, caption='Javelin results',formatColumn=fc)
print u

p = open('latex_preamble.txt','r+')
preamble_text = p.read()
p.close()

f = open('drw_results_ch_jav.tex', 'w')
f.write(preamble_text)
lines = [r'\title{DRW analysis }', r'\begin{document}', r'\maketitle', r'\section{Chelsea results}']
f.writelines(lines)
f.write(t)
f.write('\n')

f.write('\section{Javelin results}')
f.write(u)
lines=['\n', r'\end{document}', '\n']
f.writelines(lines)
f.close()