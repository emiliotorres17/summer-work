#!/usr/bin/env python3
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
call(['gfortran', '-o', 'test', 'gnu_plot.f90', 'plotter.f90'])
call(['./test'])
data    = np.loadtxt('data.dat')
x       = data[:,0]
y       = data[:,1]
plt.plot(x,y,'r')
plt.show()

