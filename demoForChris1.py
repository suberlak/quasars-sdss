# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 13:29:15 2014

@author: bmmorris
"""
import numpy as np

filename = 'data.txt'

rawtext = open(filename,'r').read().splitlines()

splitlist = []

for i in range(len(rawtext)):
    if rawtext[i].startswith("#") == False:
        rawtext_linesplit = rawtext[i].split()
        splitline = map(float,rawtext_linesplit[1:])  # avoid the first column 
        splitlist.append(splitline)

a =  np.array(splitlist)[:,:]  # as rows instead of columns

col1 = np.copy(a[:,0])  # extract the first column , taking first element of each row




'''
In [12]: a
Out[12]: 
array([[ 1.,  2.,  3.],
       [ 2.,  3.,  4.],
       [ 5.,  6.,  7.]])

In [13]: run untitled1.py

In [14]: print a
[[ 2.  3.]
 [ 3.  4.]
 [ 6.  7.]]

In [15]: print a
[[ 2.  3.]
 [ 3.  4.]
 [ 6.  7.]]

In [16]: print a[:,0
   ....: ]
[ 2.  3.  6.]

In [17]: print a[:,0]
[ 2.  3.  6.]

In [18]: print a[:,1]
[ 3.  4.  7.]

In [19]: c = a[:,0]

In [20]: a
Out[20]: 
array([[ 2.,  3.],
       [ 3.,  4.],
       [ 6.,  7.]])

In [21]: c
Out[21]: array([ 2.,  3.,  6.])

In [22]: c[0] = 0

In [23]: a 
Out[23]: 
array([[ 0.,  3.],
       [ 3.,  4.],
       [ 6.,  7.]])

In [24]: c = np.copy(a[:,0])

In [25]: c
Out[25]: array([ 0.,  3.,  6.])

In [26]: a
Out[26]: 
array([[ 0.,  3.],
       [ 3.,  4.],
       [ 6.,  7.]])

In [27]: a[0,0] = 10

In [28]: a
Out[28]: 
array([[ 10.,   3.],
       [  3.,   4.],
       [  6.,   7.]])

In [29]: c
Out[29]: array([ 0.,  3.,  6.])

In [30]: a
Out[30]: 
array([[ 10.,   3.],
       [  3.,   4.],
       [  6.,   7.]])

In [31]: b = np.copy(a)

In [32]: b
Out[32]: 
array([[ 10.,   3.],
       [  3.,   4.],
       [  6.,   7.]])

In [33]: a
Out[33]: 
array([[ 10.,   3.],
       [  3.,   4.],
       [  6.,   7.]])

In [34]: b[0] 
Out[34]: array([ 10.,   3.])
'''
