#!/usr/bin/env python
#coding: utf8

from pylab import *

data = loadtxt("output.txt", delimiter=",")
N = data.shape[1]

mass = array([1.988435e30, 1.8988e27, 5.685e26, 8.6625e25, 1.0278e26, 1.314e22])

CoM = zeros((2,N))
for i in xrange(N):
    CoM[0,i] = sum(mass * data[0::2, i]) / sum(mass)
    CoM[1,i] = sum(mass * data[1::2, i]) / sum(mass)

colors = ["Yellow", "Red", "Orange", "Blue", "Green", "Black"]

for i in xrange(6):
    scatter(data[2*i,:] - CoM[0,:], data[2*i+1,:] - CoM[1,:], 1, color=colors[i])
show()
