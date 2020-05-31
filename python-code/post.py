#!/usr/bin/env python3
import os
import sys
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Importing data                                                      #
    #---------------------------------------------------------------------#
    u       = np.loadtxt('data/u-vel_field.txt')[int(40*20-20):(40*20), 0:20]
    v       = np.loadtxt('data/u-vel_field.txt')[int(40*20-20):(40*20), 0:20]
    print(u.shape)         
    print(u[:,0])
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    M               = 20                                            # nx
    N               = 20                                            # ny
    lx              = 1.0                                           # length x-dir
    ly              = 1.0                                           # length y-dir
    dx              = lx/M                                          # x step size
    dy              = ly/N                                          # y step size
    [Xgrid, Ygrid]  = np.meshgrid(np.linspace(dx/2, lx-dx/2, M),\
                                    np.linspace(dy/2,ly-dy/2,N))
    #---------------------------------------------------------------------#
    # Calculating vorticity                                               #
    #---------------------------------------------------------------------#
    vort    = np.gradient(v, dx, edge_order=2)[0] -\
                    np.gradient(u, dx, edge_order=2)[1] 
    #---------------------------------------------------------------------#
    # plot u-velocity field                                               #
    #---------------------------------------------------------------------#
    upper   = 10.0
    lower   = -10.0
    dp      = (upper-lower)/100.0
    cnt     = plt.contourf(Xgrid, Ygrid, vort, np.arange(lower, upper+dp, dp),\
                    cmap='jet', extend='both')
    for c in cnt.collections:
        c.set_edgecolor('face')
    plt.colorbar()
    plt.show()
