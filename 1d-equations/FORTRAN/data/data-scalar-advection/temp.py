#!/usr/bin/env python3
"""========================================================================
Purpose:
    Temporary file to evaluate the different field outputs.

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Preamble                                                                #
#=========================================================================#
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import os
import sys
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------#
# User defined functions                                                  #
#-------------------------------------------------------------------------#
from ales_post.plot_settings    import plot_setting
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    M       = 1024
    dx      = 1./float(M)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    y       = loadtxt('phi-data.dat')[int(0*(M+2)):int(1*(M+2)),0]
    y0      = loadtxt('phi-data.dat')[int(0*(M+2)):int(1*(M+2)),0]
    phi0    = loadtxt('phi-data.dat')[int(0*(M+2)):int(1*(M+2)),1]
    y1      = loadtxt('phi-data.dat')[int(1*(M+2)):int(2*(M+2)),0]
    phi1    = loadtxt('phi-data.dat')[int(1*(M+2)):int(2*(M+2)),1]
    y2      = loadtxt('phi-data.dat')[int(2*(M+2)):int(3*(M+2)),0]
    phi2    = loadtxt('phi-data.dat')[int(2*(M+2)):int(3*(M+2)),1]
    y3      = loadtxt('phi-data.dat')[int(3*(M+2)):int(4*(M+2)),0]
    phi3    = loadtxt('phi-data.dat')[int(3*(M+2)):int(4*(M+2)),1]
    v       = loadtxt('v-data.dat')
    u       = loadtxt('u-data.dat')
    const   = loadtxt('constants.dat')
    time    = const[0]
    #---------------------------------------------------------------------#
    # Plotting flags                                                      #
    #---------------------------------------------------------------------#
    phi_field   = True
    v_field     = True
    u_field     = False
    p_field     = False
    #---------------------------------------------------------------------#
    # Initial conditions                                                  #
    #---------------------------------------------------------------------#
    phiE        = zeros(len(y))
    for i in range(0, len(y)):
        if y[i] < 0.5:
            phiE[i] = 0.0
        else:
            phiE[i] = 1.0
    phiE[0]     = -phiE[1]
    phiE[M+1]   = 2.0-phiE[M]
    #---------------------------------------------------------------------#
    # Phi field                                                           #
    #---------------------------------------------------------------------#
    plot_setting()
    if phi_field is True:
        plt.plot(y0, phi0, 'k--', lw=3.0, label='$\zeta(y,t=0)$')
        plt.plot(y1, phi1, 'r', lw=1.5, label='$\zeta(y,t=1)$')
        plt.plot(y2, phi2, 'b', lw=1.5, label='$\zeta(y,t=2)$')
        plt.plot(y3, phi3, 'g', lw=1.5, label='$\zeta(y,t=3)$')
        plt.legend(loc=0)
        plt.grid(True)
        plt.savefig('zeta-high-dif%4.2f.png'     %(time))
        plt.show()
    #---------------------------------------------------------------------#
    # v field                                                             #
    #---------------------------------------------------------------------#
    if v_field is True:
        y   = linspace(0.0, 1.0, M+1)
        plt.plot(y,v, 'r', lw=1.5)
        plt.title('$v$ velocity')
        plt.grid(True)
        plt.show()
    #---------------------------------------------------------------------#
    # u field                                                             #
    #---------------------------------------------------------------------#
    if u_field is True:
        y   = linspace(-0.5*dx, 1.0+0.5*dx, M+2)
        plt.plot(y, u, 'r', lw=1.5)
        plt.title('$u$ velocity')
        plt.grid(True)
        plt.show()
    #---------------------------------------------------------------------#
    # p field                                                             #
    #---------------------------------------------------------------------#
    if p_field is True:
        y       = loadtxt('pressure-data.dat')[:,0]
        pres    = loadtxt('pressure-data.dat')[:,1]
        plt.plot(y, pres, 'r', lw=1.5)
        plt.title('Pressure')
        plt.grid(True)
        plt.show()
