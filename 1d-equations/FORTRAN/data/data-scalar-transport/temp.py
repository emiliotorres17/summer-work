#!/usr/bin/env python3
import os
import sys
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
from ales_post.plot_settings    import plot_setting

call(['clear'])
M       = 128
dx      = 1./float(M)
y       = loadtxt('phi-data.dat')[:,0]
phi     = loadtxt('phi-data.dat')[:,1]
v       = loadtxt('v-data.dat')
u       = loadtxt('u-data.dat')
const   = loadtxt('constants.dat')
time    = const[0]
plot_setting()

phi_field   = True
v_field     = True
u_field     = False
p_field     = False

phi0        = zeros(len(y))
for i in range(0, len(y)):
    phi0[i]        = -8.0*(y[i]-0.5)**2.0 + 2.0 
phi0[0]     = -phi0[1]
phi0[M+1]   = -phi0[M]

if phi_field is True:
    plt.plot(y, phi0, 'k--', lw=3.0, label='$\\varphi_{0}$')
    plt.plot(y, phi, 'r', lw=1.5,\
            label='$\\varphi(t) = %4.2f$'    %(time))
    plt.title('$\phi$')
    plt.legend(loc=0)
    plt.grid(True)
    plt.savefig('phi-%4.2f.png'     %(time))
    plt.show()

if v_field is True:
    y   = linspace(0.0, 1.0, M+1)
    plt.plot(y,v, 'r', lw=1.5)
    plt.title('$v$ velocity')
    plt.grid(True)
    plt.show()

if u_field is True:
    y   = linspace(-0.5*dx, 1.0+0.5*dx, M+2)
    plt.plot(y, u, 'r', lw=1.5)
    plt.title('$u$ velocity')
    plt.grid(True)
    plt.show()

if p_field is True:
    y       = loadtxt('pressure-data.dat')[:,0]
    pres    = loadtxt('pressure-data.dat')[:,1]
    plt.plot(y, pres, 'r', lw=1.5)
    plt.title('Pressure')
    plt.grid(True)
    plt.show()
