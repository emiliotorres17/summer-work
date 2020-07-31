#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to calculate resolve 1D Navier-Stokes
    equation (i.e., u(y)).

Author:
    Emilio Torres 
========================================================================"""
#=========================================================================#
# Preamble                                                                #
#=========================================================================#
#-------------------------------------------------------------------------#
# Python preamble                                                         #
#-------------------------------------------------------------------------#
import os
import sys
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Pressure calculator                                                     #
#-------------------------------------------------------------------------#
def pressure_calculator(
        num,
        Dx,
        Dt,
        Vstar,
        P_old,
        max_error):

    """ Calculating the pressure """
    #---------------------------------------------------------------------#
    # Convergence variables                                               #
    #---------------------------------------------------------------------#
    err     = 1.
    #---------------------------------------------------------------------#
    # Preallocating space                                                 #
    #---------------------------------------------------------------------#
    P       = copy(P_old)
    res     = zeros(num+2)
    rhs     = zeros(num+2)
    rhs2    = zeros(num+2)
    #---------------------------------------------------------------------#
    # Calculating RHS                                                     #
    #---------------------------------------------------------------------#
    for j in range(0,len(Vstar)):
        rhs[j]  = 0.5*(Dx/Dt)*(v[j] - v[j-1])
        rhs2[j] = (1./Dx)*(1./Dt)*(v[j] - v[j-1])
    #---------------------------------------------------------------------#
    # Iteration loop                                                      #
    #---------------------------------------------------------------------#
    while err > max_error:
        for j in range(1,M+1):
            P[j] = 0.5*(P[j+1] + P[j-1]) - rhs[j]  
        for j in range(1,M+1):
            res[j] = (1/Dx**2.0)*(P[j+1] - 2.*P[j] + P[j-1]) - rhs2[j] 
        err = amax(abs(res))

    return P
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
    data_path   = pwd + '%cdata%c'          %(sep,sep)
    media_path  = pwd + '%cmedia%c'         %(sep,sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    M       = 128
    nu      = 0.05
    dx      = 1.0/float(M)
    rdx     = 1./dx
    rdx2    = 1./(dx**2.)
    dt      = 0.25*dx**2.0/nu
    t       = 0.0
    tfinal  = 15.0
    #---------------------------------------------------------------------#
    # Iteration variables                                                 #
    #---------------------------------------------------------------------#
    gs_error = 1.0e-08
    #---------------------------------------------------------------------#
    # Wall velocities                                                     #
    #---------------------------------------------------------------------#
    ub      = 0.0
    ut      = 0.0
    vt      = 0.0
    vb      = 0.0
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    v       = zeros(M+1)
    vstar   = zeros(M+1)
    u       = zeros(M+2)
    unew    = zeros(M+2)
    p_old   = zeros(M+2)
    #---------------------------------------------------------------------#
    # Defining boundary conditions                                        #
    #---------------------------------------------------------------------#
    u[0]        = 2.0*ub-u[1] 
    u[M+1]      = 2.0*ut-u[M] 
    unew[0]     = 2.0*ub-unew[1] 
    unew[M+1]   = 2.0*ut-unew[M] 
    v[0]        = vb
    v[M]        = vt
    #---------------------------------------------------------------------#
    # Time loop                                                           #
    #---------------------------------------------------------------------#
    count = 0
    while t < tfinal:
        #-----------------------------------------------------------------#
        # Updating time values and counters                               #
        #-----------------------------------------------------------------#
        t       += dt
        count   += 1
        #-----------------------------------------------------------------#
        # Calculating the u at the next time step                         #
        #-----------------------------------------------------------------#
        for j in range(1,M+1):  
            u1      = 0.5*(u[j] +  u[j+1])
            u2      = 0.5*(u[j] +  u[j-1])
            unew[j] = u[j] - dt*rdx*(u1*v[j] - u2*v[j-1]) \
                        + nu*dt*rdx2*(u[j+1] - 2.*u[j] + u[j-1])\
                        + dt*10.0
        unew[0]     = 2.*ub - unew[1]
        unew[M+1]   = 2.*ut - unew[M]
        u           = copy(unew)
        #-----------------------------------------------------------------#
        #-----------------------------------------------------------------#
        # Calculating v star                                             #
        #-----------------------------------------------------------------#
        for j in range(1, M):
            v1          = 0.5*(v[j+1] + v[j])
            v2          = 0.5*(v[j] + v[j-1])
            vstar[j]    = v[j] - dt*rdx*(v1**2.0 - v2**2.0) + \
                                dt*rdx2*nu*(v[j+1] - 2.*v[j] + v[j-1]) 
        #-----------------------------------------------------------------#
        # vstar boundary conditions                                       #
        #-----------------------------------------------------------------#
        vstar[0]    = vb
        vstar[M]    = vt
        #-----------------------------------------------------------------#
        # Calculating pressure                                            #
        #-----------------------------------------------------------------#
        p       = pressure_calculator(M, dx, dt, vstar, p_old, gs_error)
        p_old   = copy(p)
        #-----------------------------------------------------------------#
        # Calculating the v at the next time step                         #
        #-----------------------------------------------------------------#
        for j in range(1,M):
            v[j] = vstar[j] - (dt/dx)*(p[j+1]-p[j])
        #-----------------------------------------------------------------#
        # v boundary conditions                                           #
        #-----------------------------------------------------------------#
        v[0]    = vb
        v[M]    = vt
        #-----------------------------------------------------------------#
        # u boundary conditions                                           #
        #-----------------------------------------------------------------#
        u[0]    = 2.*ub - u[1]
        u[M+1]  = 2.*ut - u[M]
        #-----------------------------------------------------------------#
        # Print count                                                     #
        #-----------------------------------------------------------------#
        if count > 5:
            print('time --> %10.5e'                     %(t))
            print('maximum v vel. --> %10.5e'           %(amax(abs(v))))
            print('maximum u vel. --> %10.5e'           %(amax(abs(u))))
            count   = 0
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
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    y   = linspace(0.5*dx, 1.0-0.5*dx, M)
    uE  = 10./(2.*nu)*y*(1.-y)
    plt.plot(y, uE/(10./(2.*nu)), 'k--', lw=3.0, label='Exact')
    plt.plot(y, u[1:M+1]/(10.0/(2.*nu)), 'r', lw=1.5, label='Simulation')
    plt.legend(loc=0)
    plt.grid(True)
    plt.xlabel('$0 \leq y \leq$')
    plt.ylabel('$u(y)$')
    plt.savefig(media_path +  'u-poiseuille.png')
    plt.show()
    plt.close()
            
