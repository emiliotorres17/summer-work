#!/usr/bin/env python3
import os
import sys
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
from ales_post.plot_settings    import plot_setting

M       = 128
dx      = 1./float(M)
y       = loadtxt('phi-data.dat')[:,0]
phi     = loadtxt('phi-data.dat')[:,1]
v       = loadtxt('v-data.dat')
u       = loadtxt('u-data.dat')
plot_setting()
plt.plot(y, phi, 'r', lw=1.5)
plt.title('$\phi$')
plt.ylim([0.0, 2.0])
plt.grid(True)
plt.show()

y   = linspace(0.0, 1.0, M+1)
plt.plot(y,v, 'r', lw=1.5)
plt.title('$v$ velocity')
plt.grid(True)
plt.show()

y   = linspace(-0.5*dx, 1.0+0.5*dx, M+2)
plt.plot(y, u, 'r', lw=1.5)
plt.title('$u$ velocity')
plt.grid(True)
plt.show()

y       = loadtxt('pressure-data.dat')[:,0]
pres    = loadtxt('pressure-data.dat')[:,1]
plt.plot(y, pres, 'r', lw=1.5)
plt.title('Pressure')
plt.grid(True)
plt.show()
