#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to perform the post processing for the 
    finite difference solution to the wave equation.

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
import numpy as np
import matplotlib.pyplot as plt
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    sol25   = np.loadtxt('c-25.dat')
    sol50   = np.loadtxt('c-50.dat')
    sol1    = np.loadtxt('c-1.dat')
    sol     = np.loadtxt('c-temp.dat')
    IC      = np.loadtxt('IC.dat')
    u       = IC[:,0]
    x       = IC[:,1]
    plt.plot(x, sol25, 'b--', lw=1.5, label = 'c=0.25')
    plt.plot(x, sol50, 'k--', lw=1.5, label = 'c=0.50')
    plt.plot(x, sol1, 'r--', lw=1.5, label = 'c=1.0')
    plt.legend(loc=0)
    plt.grid(True)
    plt.savefig('u-velocity.png')

    sys.exit(0)

