#!/usr/bin/env python
#coding: utf8

from pylab import *
from sys import argv

if len(argv) < 2:
    filename = "output.txt"
else:
    filename = argv[1]
    
data = loadtxt(filename, delimiter=",")
N = data.shape[0]

mass = array([1.988435e30, 1.8988e27, 5.685e26, 8.6625e25, 1.0278e26, 1.314e22])

print data.shape

CoM = zeros((N,3))
for i in xrange(N):
    CoM[i,0] = sum(mass * data[i, 0::6]) / sum(mass)
    CoM[i,1] = sum(mass * data[i, 1::6]) / sum(mass)
    CoM[i,2] = sum(mass * data[i, 2::6]) / sum(mass)
    print CoM[i,:]

colors = ["Yellow", "Red", "Orange", "Blue", "Green", "Black"]

for i in xrange(6):
    print 6*i, 6*i+1
    scatter(data[:,6*i] - CoM[:,0], data[:,6*i+1] - CoM[:,1], 1, color=colors[i])
show()
