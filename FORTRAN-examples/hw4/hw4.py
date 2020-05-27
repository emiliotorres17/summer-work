#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to solve the 2D advection equation.

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
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Delta                                                                   #
#-------------------------------------------------------------------------#
def delta(
        X,
        Y):

    """ Calculating the delta function """
    #---------------------------------------------------------------------#
    # Delta                                                               #
    #---------------------------------------------------------------------#
    x0  = 1.5
    y0  = 3.0
    D   = np.sqrt((X-x0)**2.0 + (Y-y0)**2.0)
    
    return D
#-------------------------------------------------------------------------#
# P(x,y)                                                                  #
#-------------------------------------------------------------------------#
def pfunc(
        X,
        Y):

    """ Calculating the P(x,y) function """
    #---------------------------------------------------------------------#
    # P(x,y)                                                              #
    #---------------------------------------------------------------------#
    D   = delta(X,Y)
    if D <= 1.0:
        P = np.cos(0.5*np.pi*D)
        print('here')
    else:
        P = 0.0

    return P 
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    C           = 1.0
    D           = 0.5
    dx          = 0.02
    dy          = 0.02
    dt          = 0.01
    xf          = 5.0
    yf          = 5.0
    tf          = 6.0
    time        = np.arange(0.0, tf+dt, dt)
    x           = np.arange(0.0, xf+dy, dx)
    y           = np.arange(0.0, yf+dx, dy)
    [X1, X2]    = np.meshgrid(x,y)
    u           = np.zeros((len(x),len(y)))
    unew        = np.zeros((len(x),len(y)))
    #---------------------------------------------------------------------#
    # Initial condition                                                   #
    #---------------------------------------------------------------------#
    for j in range(0,len(y)):
        for i in range(0,len(x)):
            u[i,j]  = pfunc(x[i], y[j])
            print(u[i,j])
    #---------------------------------------------------------------------#
    # Plotting the initial conditions                                     #
    #---------------------------------------------------------------------#
    print(np.amax(u))
    c1      = 1.0
    upper   = c1
    lower   = 0.0
    cnt     = plt.contourf(X1, X2, np.transpose(u), np.arange(lower, upper,\
                            (upper-lower)/1000.0), cmap='jet', extend='both')
    for c in cnt.collections:
        c.set_edgecolor('face')
    plt.colorbar(cnt)
    plt.savefig('IC.png')
    plt.clf()
    #---------------------------------------------------------------------#
    # Marching thru time                                                  #
    #---------------------------------------------------------------------#
    for t in range(1, len(time)):
        for j in range(1,len(y)):
            for i in range(1, len(x)):
                unew[i,j]  = u[i,j] - dt*(C*(u[i,j]-u[i-1,j])/dx + D*(u[i,j]-u[i,j-1])/dy)
        unew[0,:]   = unew[-1,:]
        unew[:,0]   = unew[:,-1]
        u           = np.copy(unew)
    #------------------u---------------------------------------------------#
    # Plotting the t = 0.2                                                #
    #---------------------------------------------------------------------#
    print(np.amax(u))
    #cnt     = plt.contour(X1, X2, np.transpose(u), np.arange(lower, upper,\
    #                        (upper-lower)/500.0), cmap='jet', extend='both')
    cnt     = plt.contour(X1, X2, np.transpose(u), 10, cmap='jet', extend='both')
    #for c in cnt.collections:
    #    c.set_edgecolor('face')
    plt.colorbar(cnt)
    plt.savefig('test.png')
    plt.clf()