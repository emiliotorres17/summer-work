#!/usr/bin/env python3
"""========================================================================
Plotting velocity for Poiseulle flow 
========================================================================"""
import sys
from numpy import *
import matplotlib.pyplot as plt
if __name__ == '__main__':
    M       = 64
    u       = loadtxt('u-temp.dat')[1:M+1,32]
    dx      = 1.0/float(M)
    G       = 1.0
    mu      = 0.1
    h       = 1.0
    norm    = (G*h**2.0)/(2.0*mu)
    u       = u/norm
    y       = linspace(0.5*dx, 1.0-0.5*dx, M) 
    x       = linspace(0.0, 1.0, M+1) 
    u_exact = G/(2.0*mu)*y*(h-y)
    u_exact = u_exact/norm
    #---------------------------------------------------------------------#
    # Font settings                                                       #
    #---------------------------------------------------------------------#
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    SMALL_SIZE = 14
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 12
    plt.rc('font',      size=SMALL_SIZE)            # controls default text sizes
    plt.rc('axes',      titlesize=SMALL_SIZE)       # fontsize of the axes title
    plt.rc('axes',      labelsize=MEDIUM_SIZE)      # fontsize of the x and y labels
    plt.rc('xtick',     labelsize=SMALL_SIZE)       # fontsize of the tick labels
    plt.rc('ytick',     labelsize=SMALL_SIZE)       # fontsize of the tick labels
    plt.rc('legend',    fontsize=SMALL_SIZE)        # legend fontsize
    plt.rc('figure',    titlesize=BIGGER_SIZE)      # fontsize of the figure title
    #---------------------------------------------------------------------#
    # Font settings                                                       #
    #---------------------------------------------------------------------#
    plt.plot(y, u, 'r', lw=1.5,label='Simulation')
    plt.plot(y, u_exact, 'go', markevery=8, lw=1.5, label='Exact')
    plt.grid(True)
    plt.ylabel('$\\frac{u}{Gh^{2}/2\\mu}$')
    plt.xlabel('$y/h$')
    plt.tight_layout()
    plt.legend(loc=0)
    plt.savefig('u-poiseuille.png')
    plt.close()
