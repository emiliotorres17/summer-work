#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to build and execute the FORTRAN 
    Navier-Stokes solver.
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
    data_path   = pwd + '%cFORTRAN-data%c'      %(sep, sep)
    if os.path.exists(data_path) is False:
        os.mkdir(data_path)
    #---------------------------------------------------------------------#
    # Generating executables                                              #
    #---------------------------------------------------------------------#
    call1   = call(['gfortran', '-c', 'precision_m.f90'])
    if call1 == 1:
        print('**** Compiling error --> precision_m.f90')
        sys.exit(1)
    call2   = call(['gfortran', '-c', 'navier_stokes_module.f90'])
    if call2 == 1:
        print('**** Compiling error --> navier_stokes_module.f90')
        sys.exit(1)
    call3   = call(['gfortran','-o', 'ns_solve', 'navier_stokes_solver2.f90',\
            'navier_stokes_module.f90', 'precision_m.f90'])
    if call3 == 1:
        print('**** Compiling error --> navier_stokes_solver2.f90')
        sys.exit(1)
    #---------------------------------------------------------------------#
    # Executing FORTRAN                                                   #
    #---------------------------------------------------------------------#
    call(['./ns_solve'])

    print('**** Sucessful run ****')
    sys.exit(0)

