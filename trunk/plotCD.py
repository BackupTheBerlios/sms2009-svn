#!/usr/bin/env python

from numpy import loadtxt,shape,linspace
from pylab import plot,show,clf,show,ion

# Import
d=loadtxt('data')
dims = shape(d)

clf()
ion()
x = linspace(0,1,dims[1])
ph,=plot(x,d[0,:],'k')
for i in range(1,dims[0]):
    ph.set_ydata(d[i,:])
    show()
