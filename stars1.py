# read in the SDSS Stripe 82 for Standard Stars, 
# ignore the commented lines,
# and ignore the first column, which only has 
# 'CALIBRATED' keyword

import numpy as np

filename='stripe82calibStars_v2.6.dat'

rawtext = open(filename,'r').read().splitlines()


# probe the first hundred lines and check how many comment lines there are 
comment=0
for  i in range(100):
	if rawtext[i].startswith("#") == True:
           comment += 1
print comment

# cut out all the comment lines
rawtext=rawtext[comment:]

# initiate an array only for the number of rows that include numbers
my_array = np.zeros((len(rawtext),36))
counter = 0
i=0
for i in range(len(rawtext)):
    rawtext_linesplit = rawtext[i].split()
    splitline = map(float,rawtext_linesplit[1:])  # avoid the first column 
    my_array[i,:] = splitline
    counter += 1



print counter

np.save('my_array1.npy',my_array)


# now my_array has a structure   (rows,columns)  
