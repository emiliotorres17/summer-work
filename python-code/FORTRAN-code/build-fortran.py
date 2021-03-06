#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is build and execute the FORTRAN code, and 
    various plots.

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
import time
import numpy as np
import matplotlib.pyplot as plt
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #=====================================================================#
    # Main preamble                                                       #
    #=====================================================================#
    call(['clear'])
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + '%c..%cdata%c'                  %(sep, sep, sep) 
    media_path  = pwd + '%c..%cmedia%c'                 %(sep, sep, sep) 
    #=====================================================================#
    # Performance flags                                                   #
    #=====================================================================#
    fortran_flag    = True
    vorticity_flag  = True
    #=====================================================================#
    # Executing FORTRAN                                                   #
    #=====================================================================#
    if fortran_flag is True:
        call(['gfortran', '-o', 'ns_solve', 'navier_stokes.f90'])
        print('Compiled ---> navier_stokes.f90')
        time.sleep(3)
        call(['./ns_solve'])
    #=====================================================================#
    # Post processing                                                     #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    omega   = np.loadtxt(data_path + 'vorticity.dat')
    omega   = np.transpose(omega)
    dim     = omega.shape
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    x       = np.linspace(0.0, 1.0, dim[0])
    y       = np.linspace(0.0, 1.0, dim[1])
    [X,Y]   = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Plot settings                                                       #
    #---------------------------------------------------------------------#
    lower   = 0.0
    upper   = 1.0
    dp      = (upper-lower)/10.0
    con     = [0.0, -0.5, 0.5, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0, -4.0, -5.0]
    con     = np.array([-5.0, -4.0, -3.0, -2.0, -1.0, -0.5, 0.0, 0.5,\
                        1.0, 2.0, 3.0])
    #---------------------------------------------------------------------#
    # Generating contour plots                                            #
    #---------------------------------------------------------------------#
    if vorticity_flag is True:
        cnt     = plt.contour(X,Y, omega, con,\
                                cmap='jet', extend='both')
        for c in cnt.collections:
            c.set_edgecolor('face')
        plt.colorbar()
        plt.savefig(media_path + 'vorticity-' + str(dim[0]-1) + '.png')
        plt.close()
