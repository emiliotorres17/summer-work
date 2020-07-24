#!/usr/bin/env python3
"""========================================================================
Purpose:
    Post processing the results for the Taylor-Green vortex problem 

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Python packages                                                         #
#=========================================================================#
import os 
import sys
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    sep     = os.sep
    pwd     = os.getcwd()
    IC_path = pwd + '%c..%cIC-data%c'               %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    xv      = loadtxt(IC_path + 'v-x.dat')
    yv      = loadtxt(IC_path + 'v-y.dat')
    v       = loadtxt('v-temp.dat')
    xu      = loadtxt(IC_path + 'u-x.dat')
    yu      = loadtxt(IC_path + 'u-y.dat')
    u       = loadtxt('u-temp.dat')
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    M       = 128
    t       = 1.00483E+00 
    nu      = 0.001
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
    plt.tight_layout()
    #---------------------------------------------------------------------#
    # v-velocity plot                                                     #
    #---------------------------------------------------------------------#
    yint    = 64
    v       = v[yint,:]
    vE      = -sin(xv)*cos(yv[yint])*exp(-2.*nu*t) 
    plt.plot(xv, v, 'r', lw=1.5, label='Simulation')
    plt.plot(xv, vE, 'bo', markevery=8, lw=1.5, label='Exact')
    plt.grid(True)
    plt.xlabel('$0\leq x \leq 2\pi$')
    plt.ylabel('v')
    plt.xlim([0.0, 2.0*pi])
    plt.legend(loc=0)
    plt.savefig('v-velocity.png')
    plt.close()
    #---------------------------------------------------------------------#
    # u-velocity plot                                                     #
    #---------------------------------------------------------------------#
    xint    = 64
    u       = u[:,xint]
    uE      = cos(xu[xint])*sin(yu)*exp(-2.*nu*t) 
    plt.plot(yu, u, 'r', lw=1.5, label='Simulation')
    plt.plot(yu, uE, 'bo', markevery=8, lw=1.5, label='Exact')
    plt.grid(True)
    plt.xlabel('$0\leq y \leq 2\pi$')
    plt.ylabel('u')
    plt.xlim([0.0, 2.0*pi])
    plt.legend(loc=0)
    plt.savefig('u-velocity.png')
    plt.close()
