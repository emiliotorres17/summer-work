#!/usr/bin/env python3
"""========================================================================
Purpose:
    Perform the post processing for the Crank-Nickolsen and BTCS time
    stepping methods.

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
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    #---------------------------------------------------------------------#
    # Path variables                                                      #
    #---------------------------------------------------------------------#
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + '%cdata%c'      %(sep, sep)  
    media_path  = pwd + '%cmedia%c'     %(sep, sep)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    data    = loadtxt(data_path + 'u-exact.dat')
    data2   = loadtxt(data_path + 'u-IC.dat')
    data3   = loadtxt(data_path + 'u-approx.dat')
    data4   = loadtxt(data_path + 'u-cn.dat')
    u       = data[:,1]
    x       = data[:,0]
    uFD     = data3[:,1]
    xFD     = data3[:,0]
    uCN     = data4[:,1]
    xCN     = data4[:,0]
    #---------------------------------------------------------------------#
    # Plot settings                                                       #
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
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    plt.plot(xFD, uFD, 'ko--', lw=1.5, label =  'BTCS')
    plt.plot(xCN, uCN, 'bo--', lw=1.5, label =  'CN')
    plt.plot(x, u, 'r', lw=1.5, label = 'Exact')
    plt.legend(loc=0)
    plt.grid(True)
    plt.savefig(media_path + 'CN-BTCS-1d-heat.png')
    plt.close()

    print('**** Successful run ****') 
    sys.exit(0)
