#  load and plot  stripe82calibStars_v2.6.dat

# " produce rms vs. magnitude plot using CRTS data (to show that their [Graham et al. 2011] fig. 10 is bogus - that variation # is not due to quasar variability but due to their errors)"


import numpy as np

import matplotlib.pyplot as plt

arr=np.load('my_array1.npy')

# u  mmu (mean magnitude)  vs  u  rms (root mean square) 

# plt.title('u-band rms vs. magnitude')
# plt.xlabel('Mean u magnitude')
# plt.ylabel('Root-mean-square scatter')
# plt.scatter(arr[1:1000,8], arr[1:1000,10])
# plt.xlim(15,)  # print only those dimmer than 15th mag 
# plt.savefig('u_band_mean_mag.png')
# plt.show()

# plot ugriz as separate plots
nrows=1006848
band_names=['u','g','r','i','z']
band_cols_mmu=[8,14,20,26,32] # columns with mmu : mean magnitude
band_cols_med=[7,13,19,25,31] # columns with median magnitude 

# col_med + 3 : rms
# col_mmu + 2 : rms  

i=0
for col in band_cols_mmu:
    # plt.xlim(15,)  # print only those dimmer than 15th mag 
    plt.clf()
    plt.title(band_names[i]+'-band rms vs. mean magnitude')
    plt.xlabel('Mean '+band_names[i]+' magnitude')
    plt.ylabel('Root-mean-square scatter')
    plt.scatter(arr[0:nrows,col],arr[0:nrows,col+2])
    fname='mean_'+str(nrows)+'_rows_'+band_names[i]+'.png'
    plt.savefig(fname)
    #plt.show()
    i=i+1
