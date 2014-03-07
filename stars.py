# read in the SDSS Stripe 82 for Standard Stars, 
# ignore the commented lines,
# and ignore the first column, which only has 
# 'CALIBRATED' keyword

import numpy as np

filename='stripe82calibStars_v2.6.dat'

rawtext = open(filename,'r').read().splitlines()

#splitlist = []
my_array = np.zeros((1006808,36))
counter = 0
j=0
for i in range(1000):
    if rawtext[i].startswith("#") == False:
        rawtext_linesplit = rawtext[i].split()
        splitline = map(float,rawtext_linesplit[1:])  # avoid the first column 
        #splitlist.append(splitline)
	my_array[j,:] = splitline
        j += 1
	counter += 1

comment = 0
for  i in range(1000):
	if rawtext[i].startswith("#") == True:
           comment += 1
print comment

print counter

#np.save('my_array1.npy',my_array)


# now my_array has a structure   (rows,columns)  
