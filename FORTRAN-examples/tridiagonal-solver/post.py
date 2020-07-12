#!/usr/bin/env python3
"""========================================================================
Purpose:
    Post processing for BTCS example.

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
    # Path variables                                                      #
    #---------------------------------------------------------------------#
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + '%cdata%c'              %(sep, sep)
    media_path  = pwd + '%cmedia%c'             %(sep, sep)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    u_0         = loadtxt(data_path + 'u-IC.dat')
    u_f         = loadtxt(data_path + 'u-final.dat') 
    u_18        = loadtxt(data_path + 'u-18.dat') 
    u_36        = loadtxt(data_path + 'u-36.dat') 
    u_54        = loadtxt(data_path + 'u-54.dat') 
    u_72        = loadtxt(data_path + 'u-72.dat') 
    u_90        = loadtxt(data_path + 'u-90.dat') 
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    N           = 64
    y           = linspace(0.0, 0.4, N+1) 
    #---------------------------------------------------------------------#
    # Plotting settings                                                   #
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
    plt.plot(u_0, y, 'k',  lw=1.5, label='$t=0$')
    plt.plot(u_18, y, 'b', lw=1.5, label='$t=0.18$')
    plt.plot(u_36, y, 'g', lw=1.5, label='$t=0.36$')
    plt.plot(u_54, y, 'm', lw=1.5, label='$t=0.54$')
    plt.plot(u_72, y, 'c', lw=1.5, label='$t=0.72$')
    plt.plot(u_90, y, 'y', lw=1.5, label='$t=0.90$')
    plt.plot(u_f, y, 'r',  lw=1.5, label='$t=1.0$')
    plt.legend(loc=0)
    plt.grid(True)
    plt.xlabel('$u$ (m/sec)')
    plt.ylabel('$y$ (m)')
    plt.savefig(media_path + 'laaason-example-plot.png')
    plt.close()


    print('**** Successful run *****')
    sys.exit()
