#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to perform the post processing for the
    SOR method that we are going to implement into NS code.

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
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#       
    # Main preamble                                                       #
    #---------------------------------------------------------------------#       
    call(['clear'])
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + '%cdata%c'          %(sep, sep)
    #---------------------------------------------------------------------#       
    # Loading data                                                        #
    #---------------------------------------------------------------------#       
    phi         = loadtxt(data_path + 'SOR-data.dat')
    res         = loadtxt(data_path + 'residual-data.dat')
    #---------------------------------------------------------------------#       
    # Domain variables                                                    #
    #---------------------------------------------------------------------#       
    M           = int(phi.shape[0]-1)
    x           = linspace(-1.0, 1.0, M+1)
    y           = linspace(-1.0, 1.0, M+1)
    [X, Y]      = meshgrid(x,y)
    iters       = linspace(1, len(res), len(res))
    #---------------------------------------------------------------------#       
    # 3D surf plot                                                        #
    #---------------------------------------------------------------------#       
    fig         = plt.figure()
    ax          = fig.gca(projection='3d')
    surf        = ax.plot_surface(X, Y, phi, cmap='jet')
    ax.set_zlim(-0.8, 0.6)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    #---------------------------------------------------------------------#       
    # Iteration versus line plot                                          #
    #---------------------------------------------------------------------#       
    plt.semilogy(iters, res, 'r', lw=1.5, label='$L_{\infty} of residual')
    plt.grid(True)
    plt.ylim([1e-12, 1e07])
    plt.show()

    print('**** Successful run ****')
    sys.exit(0)

