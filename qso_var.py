# read in the list of files, open the file from the list, and read the content
# average the values of the second column and do some operations on the third 
# column for those rows which have the same few numbers of the first column (i.e.
# the same 5 digits of the MJD, which signifies the same day). 

import numpy as np
np.seterr(invalid='raise')

names=np.loadtxt('QSO_try/file.list',dtype=str)


# open the first one to try

address='QSO_try/'+names[2]
data=np.loadtxt(address)

# want to include something that saves as a separate list the names of empty 
# files, but I have no idea how... Those links : 
# https://docs.python.org/2/library/warnings.html
# 
# 
# do not help 

j=0
averages=np.zeros(shape=(len(data),3))

# ask - how do I create an empty array to which I can then append something, 
# without knowing it's size in advance ?

mjd=[]
mag_day=[]
mag_m=[]
mag_m_err=[]

while j != 0:
    for i in range(len(data[:,0]):
        if str(data[i+1,0])[:5] == str(data[i,0])[:5] :
           mag_day=np.append(mag_day,data[i,1]
    
    
mjd=np.append(mjd,data[i,0])


# ?? how do I save those few values of magnitude for the same night, to 
# calculate their averages ? 

