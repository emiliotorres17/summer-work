#!/usr/bin/env python3
"""========================================================================
Purpose:
    py_471proj3.py: Spring 2005 MAE471 2D Lid-Driven Cavity Navier-Stokes
    Project Converted to Python

Author:
    Jon Baltzer
========================================================================"""
#=========================================================================#
# Preamble                                                                #
#=========================================================================#
__author__    = "Jon Baltzer"
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import os
from subprocess import call
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#-------------------------------------------------------------------------#
# Provisional time step                                                   #
#-------------------------------------------------------------------------#
def ddt_prov(M, N, dx, dy, nu, u, v, ul, ur, ub, ut, vl, vr, vb, vt):

    """ calculates time rate of change of provisional velocity field u*,v* 
        (without pressure term) """

    dx2 = dx*dx
    dy2 = dy*dy

    Lu = np.zeros((N+1,M))
    Lv = np.zeros((N,M+1))
    Nx = np.zeros((N+1,M))
    Ny = np.zeros((N,M+1))

    # linear terms for u

    # Lower left wall
    Lu[1][1] = nu*(u[1][2]-2.*u[1][1]+ul[1])/dx2
    Lu[1][1] += nu*(u[2][1]-3.*u[1][1]+2.*ub[1])/dy2

    for i in range(2, M-1):
        Lu[1][i] = nu*(u[1][i+1]-2.*u[1][i]+u[1][i-1])/dx2
        Lu[1][i] += nu*(u[2][i]-3.*u[1][i]+2.*ub[i])/dy2

    # Lower right wall
    Lu[1][M-1] = nu*(ur[1]-2.*u[1][M-1]+u[1][M-2])/dx2
    Lu[1][M-1] += nu*(u[2][M-1]-3.*u[1][M-1]+2.*ub[M-1])/dy2

    for j in range(2, N):
        # left wall
        Lu[j][1] = nu*(u[j][2]-2.*u[j][1]+ul[j])/dx2
        Lu[j][1] += nu*(u[j+1][1]-2.*u[j][1]+u[j-1][1])/dy2

        for i in range(2, M-1):
            Lu[j][i] = nu*(u[j][i+1]-2.*u[j][i]+u[j][i-1])/dx2
            Lu[j][i] += nu*(u[j+1][i]-2.*u[j][i]+u[j-1][i])/dy2

        # right wall
        Lu[j][M-1] = nu*(ur[j]-2.*u[j][M-1]+u[j][M-2])/dx2
        Lu[j][M-1] += nu*(u[j+1][M-1]-2.*u[j][M-1]+u[j-1][M-1])/dy2

    # Upper left wall
    Lu[N][1] = nu*(u[N][2]-2.*u[N][1]+ul[N])/dx2
    Lu[N][1] += nu*(2.*ut[1]-3.*u[N][1]+u[N-1][1])/dy2

    for i in range(2, M-1):
        Lu[N][i] = nu*(u[N][i+1]-2.*u[N][i]+u[N][i-1])/dx2
        Lu[N][i] += nu*(2.*ut[i]-3.*u[N][i]+u[N-1][i])/dy2

    # Upper right wall
    Lu[N][M-1] = nu*(ur[N]-2.*u[N][M-1]+u[N][M-2])/dx2
    Lu[N][M-1] += nu*(2.*ut[M-1]-3.*u[N][M-1]+u[N-1][M-1])/dy2


    # linear terms for v

    # lower left
    Lv[1][1] = nu*(v[1][2]-3.*v[1][1]+2.*vl[1])/dx2
    Lv[1][1] += nu*(v[2][1]-2.*v[1][1]+vb[1])/dy2

    for i in range(2, M):
        Lv[1][i] = nu*(v[1][i+1]-2.*v[1][i]+v[1][i-1])/dx2
        Lv[1][i] += nu*(v[2][i]-2.*v[1][i]+vb[i])/dy2

    # Lower right
    Lv[1][M] = nu*(2.0*vr[1]-3.*v[1][M]+v[1][M-1])/dx2
    Lv[1][M] += nu*(v[2][M]-2.*v[1][M]+vb[M])/dy2

    for j in range(2, N-1):
        # left wall
        Lv[j][1] = nu*(v[j][2]-3.*v[j][1]+2.*vl[j])/dx2
        Lv[j][1] += nu*(v[j+1][1]-2.*v[j][1]+v[j-1][1])/dy2

        for i in range(2, M):
            Lv[j][i] = nu*(v[j][i+1]-2.*v[j][i]+v[j][i-1])/dx2
            Lv[j][i] += nu*(v[j+1][i]-2.*v[j][i]+v[j-1][i])/dy2

        # right wall
        Lv[j][M] = nu*(2.*vr[j]-3.*v[j][M]+v[j][M-1])/dx2
        Lv[j][M] += nu*(v[j+1][M]-2.*v[j][M]+v[j-1][M])/dy2

    # Upper left wall
    Lv[N-1][1] = nu*(v[N-1][2]-3.*v[N-1][1]+2.*vl[N-1])/dx2
    Lv[N-1][1] += nu*(vt[1]-2.*v[N-1][1]+v[N-2][1])/dy2

    for i in range(2, M):
        Lv[N-1][i] = nu*(v[N-1][i+1]-2.*v[N-1][i]+v[N-1][i-1])/dx2
        Lv[N-1][i] += nu*(vt[i]-2.*v[N-1][i]+v[N-2][i])/dy2

    # Upper right wall
    Lv[N-1][M] = nu*(2.*vr[N-1]-3.*v[N-1][M]+v[N-1][M-1])/dx2
    Lv[N-1][M] += nu*(vt[M]-2.*v[N-1][M]+v[N-2][M])/dy2



    # nonlinear terms for u

    for j in range(1, N+1):
        # left wall
        Nx[j][1] = (0.25*(u[j][2]+u[j][1])**2 - 0.25*(u[j][1]+ul[j])**2)/dx

        for i in range(2, M-1):
            Nx[j][i] = (0.25*(u[j][i+1] + u[j][i])**2 - 0.25*(u[j][i] + u[j][i-1])**2)/dx

        # right wall
        Nx[j][M-1] = (0.25*(ur[j]+u[j][M-1])**2 - 0.25*(u[j][M-1]+u[j][M-2])**2)/dx

    for i in range(1, M):
        Nx[1][i] += (0.25*(u[2][i]+u[1][i])*(v[1][i+1]+v[1][i]) - 0.5*(ub[i])*(vb[i]+vb[i+1]))/dy

    for j in range(2, N):
        for i in range(1, M):
            Nx[j][i] += (0.25*(u[j+1][i]+u[j][i])*(v[j][i+1]+v[j][i]) - 0.25*(u[j][i]+u[j-1][i])*(v[j-1][i]+v[j-1][i+1]))/dy

    for i in range(1, M):
        Nx[N][i] += (0.5*(ut[i])*(vt[i+1]+vt[i]) - 0.25*(u[N][i]+u[N-1][i])*(v[N-1][i]+v[N-1][i+1]))/dy


    # nonlinear terms for v

    for i in range(1, M+1):
        # left wall
        Ny[1][i] = (0.25*(v[2][i]+v[1][i])**2 - 0.25*(v[1][i]+vb[i])**2)/dy

    for j in range(2, N-1):
        for i in range(1, M+1):
            Ny[j][i] = (0.25*(v[j+1][i] + v[j][i])**2 - 0.25*(v[j][i] + v[j-1][i])**2)/dy

    for i in range(1, M+1):
        # right wall
        Ny[N-1][i] = (0.25*(vt[i]+v[N-1][i])**2 - 0.25*(v[N-1][i]+v[N-2][i])**2)/dy

    for j in range(1, N):
        Ny[j][1] += (0.25*(u[j][1]+u[j+1][1])*(v[j][2]+v[j][1]) - 0.5*(ul[j]+ul[j+1])*(vl[j]))/dx

        for i in range(2, M):
            Ny[j][i] += (0.25*(u[j][i]+u[j+1][i])*(v[j][i+1]+v[j][i]) - 0.25*(u[j][i-1]+u[j+1][i-1])*(v[j][i]+v[j][i-1]))/dx

        Ny[j][M] += (0.5*(ur[j]+ur[j+1])*(vr[j]) - 0.25*(u[j][M-1]+u[j+1][M-1])*(v[j][M]+v[j][M-1]))/dx

    return [Lu, Lv, Nx, Ny]
#-------------------------------------------------------------------------#
# Calculating pressure                                                    #
#-------------------------------------------------------------------------#
def calcpress(M, N, dx, dy, maxerror, ustar, vstar, ul, ur, ub, ut, vl, vr, vb, vt):

    """ Calculating the pressure """

    gsiter = 0
    normnew = 0.
    diffnorm = 0.
    normerror = 2.*maxerror      # normerror is relative error in norm of diagonal
    rhs = 0.
    dx2 = dx*dx
    dy2 = dy*dy

    P = np.zeros((N+1,M+1))
    Pold = np.zeros((N+1,M+1))

    #while gsiter<5000:
    while((normerror > maxerror and gsiter < 5000) or (gsiter < 30)):
        
        #normold = normnew
        normnew = 0.
        diffnorm = 0.        

        gsiter += 1

        # lower left
        rhs = (ustar[1][1]-ul[1])/dx + (vstar[1][1]-vb[1])/dy
        P[1][1] = 1./(1./dx2+1./dy2)*(P[1][2]/dx2+P[2][1]/dy2-rhs)

        # lower
        for i in range(2, M):
            rhs = (ustar[1][i]-ustar[1][i-1])/dx + (vstar[1][i]-vb[i])/dy
            P[1][i] = 1./(2./dx2+1./dy2)*((P[1][i+1]+P[1][i-1])/dx2+P[2][i]/dy2-rhs)

        # lower right
        rhs = (ur[1]-ustar[1][M-1])/dx + (vstar[1][M]-vb[M])/dy
        P[1][M] = 1./(1./dx2+1./dy2)*(P[1][M-1]/dx2+P[2][M]/dy2-rhs)

        for j in range(2, N):
            # left
            rhs = (ustar[j][1]-ul[j])/dx + (vstar[j][1]-vstar[j-1][1])/dy
            P[j][1] = 1./(1./dx2+2./dy2)*(P[j][2]/dx2+(P[j+1][1]+P[j-1][1])/dy2-rhs)

            # away from the boundaries
            for i in range(2, M):
                rhs = (ustar[j][i]-ustar[j][i-1])/dx + (vstar[j][i]-vstar[j-1][i])/dy
                P[j][i] = 1./(2./dx2+2./dy2)*((P[j][i+1]+P[j][i-1])/dx2+(P[j+1][i]+P[j-1][i])/dy2-rhs)

            # right
            rhs = (ur[j]-ustar[j][M-1])/dx + (vstar[j][M]-vstar[j-1][M])/dy
            P[j][M] = 1./(1./dx2+2./dy2)*(P[j][M-1]/dx2+(P[j+1][M]+P[j-1][M])/dy2-rhs)

        # upper left
        rhs = (ustar[N][1]-ul[N])/dx + (vt[1]-vstar[N-1][1])/dy
        P[N][1] = 1./(1./dx2+1./dy2)*(P[N][2]/dx2+P[N-1][1]/dy2-rhs)

        # upper
        for i in range(2, M):
            rhs = (ustar[N][i]-ustar[N][i-1])/dx + (vt[i]-vstar[N-1][i])/dy
            P[N][i] = 1./(2./dx2+1./dy2)*((P[N][i+1]+P[N][i-1])/dx2+P[N-1][i]/dy2-rhs)

        # upper right
        rhs = (ur[N]-ustar[N][M-1])/dx + (vt[M]-vstar[N-1][M])/dy
        P[N][M] = 1./(1./dx2+1./dy2)*(P[N][M-1]/dx2+P[N-1][M]/dy2-rhs)

        for j in range(1,N+1):
            for i in range(1,M+1):
                normnew += P[j][i]**2
                diffnorm += (P[j][i]-Pold[j][i])**2

        # normnew = 
        # normerror = abs((normnew-normold)/normnew)
        normerror = diffnorm/normnew

        Pold = np.copy(P)    

    print("%16.5e: %5i iterations" % (normerror, gsiter))

    return P
#=========================================================================#
# Main code                                                               #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Path variables                                                      #
    #---------------------------------------------------------------------#
    call(['clear'])
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + '%cdata%c'          %(sep, sep) 
    media_path  = pwd + '%cmedia%c'         %(sep, sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    M   = 20                                # nx
    N   = 20                                # ny
    lx  = 1.0                               # length x-dir
    ly  = 1.0                               # length y-dir
    dx  = lx/M                              # x step size
    dy  = ly/N                              # y step size
    #---------------------------------------------------------------------#
    # Writing and plotting flags                                          #
    #---------------------------------------------------------------------#
    write_vel_field         = True
    write_terms             = True
    write_P                 = True
    write_divergence_field  = True
    plot_vel_field          = True
    #---------------------------------------------------------------------#
    # Preallocation variables                                             #
    #---------------------------------------------------------------------#
    u       = np.zeros((N+1,M))             # u-velocity field
    v       = np.zeros((N,M+1))             # v-velocity filed
    ustar   = np.zeros((N+1,M))             # u-star non pressure term
    vstar   = np.zeros((N,M+1))             # v-star non pressure term
    ul      = np.zeros((N+1))               # u-left wall boundary velocity
    ur      = np.zeros((N+1))               # u-right wall boundary velocity
    vl      = np.zeros((N+1))               # v-left wall boundary velocity  
    vr      = np.zeros((N+1))               # v-right wall boundary velocity 
    ub      = np.zeros((M+1))               # u-bottom wall velocity
    ut      = np.zeros((M+1))               # u-top wall velocity
    vb      = np.zeros((M+1))               # v-bottom wall velocity
    vt      = np.zeros((M+1))               # v-top wall velocity
    #---------------------------------------------------------------------#
    # Defining the error variables                                        #
    #---------------------------------------------------------------------#
    maxerror    = 1.e-12            # not sure this specific variables
    diverg      = 0.
    maxdiverg   = 0.
    #---------------------------------------------------------------------#
    # set up BC's, velocity components at boundaries are all set to 0     #
    # until now -- lid-driven cavity                                      #
    # I think these are the left and right boundaries
    #---------------------------------------------------------------------#
    #---------------------------------------------------------------------#
    # set up BC's, velocity components at boundaries are all set to 0     #
    # until now -- unidirectional flow                                    #
    #---------------------------------------------------------------------#
    for i in range(1,M+1):
        ut[i] = 1.                  # setting the u-vel of top wall to 1 
        ub[i] = 0.0                 # setting the u-vel of bottom wall to 1
    for j in range(1,N+1):
        vr[j] = 0.0                 # not sure of these are correct
        vl[j] = 0.0                 # I think it should vl and ur
    #---------------------------------------------------------------------#
    # Setting the time step variables                                     #
    #---------------------------------------------------------------------#
    nu      = 0.01                  # viscosity
    dt      = 0.25*dx*dx/nu/2.      # time step
    tmax    = 40.0                  # final time
    #---------------------------------------------------------------------#
    # velocity output files                                               #
    #---------------------------------------------------------------------#
    if write_vel_field:
        #-----------------------------------------------------------------#
        # u-velocity field                                                #
        #-----------------------------------------------------------------#
        if os.path.exists(data_path + 'u-vel_field.txt'):
            os.remove(data_path + 'u-vel_field.txt')
        u_fvelfld   = open(data_path + 'u-vel_field.txt', 'w')
        #-----------------------------------------------------------------#
        # v-velocity field                                                #
        #-----------------------------------------------------------------#
        if os.path.exists(data_path + 'v-vel_field.txt'):
            os.remove(data_path + 'v-vel_field.txt')
        v_fvelfld   = open(data_path + 'v-vel_field.txt', 'w')
    #---------------------------------------------------------------------#
    # NEED TO COMMENT                                                     #
    #---------------------------------------------------------------------#
    if write_terms:
        if os.path.exists(data_path + 'terms.txt'):
            os.remove(data_path + 'terms.txt')
        fterms      = open(data_path + 'terms.txt', 'w')
    #---------------------------------------------------------------------#
    # Pressure output files                                               #
    #---------------------------------------------------------------------#
    if write_P:
        if os.path.exists(data_path + 'P.txt'):
            os.remove(data_path + 'P.txt')
        fP      = open(data_path + 'P.txt', 'w')
    #---------------------------------------------------------------------#
    # Divergence field                                                    #
    #---------------------------------------------------------------------#
    if write_divergence_field:
        if os.path.exists(data_path + 'divg_field.txt'):
            os.remove(data_path + 'divg_field.txt')
        fdiverg = open(data_path + 'divg_field.txt', 'w')
    #---------------------------------------------------------------------#
    # Setting the cell centered mesh grid                                 #
    #---------------------------------------------------------------------#
    [Xgrid, Ygrid]      = np.meshgrid(np.linspace(dx/2, lx-dx/2, M),\
                                    np.linspace(dy/2,ly-dy/2,N))
    #---------------------------------------------------------------------#
    # main time stepping loop                                             #
    #---------------------------------------------------------------------#
    for nstep in range(1, int(1.0000001*tmax/dt)+1):
        [Lu, Lv, Nx, Ny]    = ddt_prov(M, N, dx, dy, nu, u, v, ul, ur, ub,\
                                            ut, vl, vr, vb, vt)
        #-----------------------------------------------------------------#
        # u-star and v-star solution                                      #
        #-----------------------------------------------------------------#
        for j in range(1, N+1):
            for i in range(1, M):
                ustar[j][i] = u[j][i] + dt*(-Nx[j][i] + Lu[j][i])
        for j in range(1, N):
            for i in range(1, M+1):
                vstar[j][i] = v[j][i] + dt*(-Ny[j][i] + Lv[j][i])
        #-----------------------------------------------------------------#
        # Calculating the pressure term                                   #
        #-----------------------------------------------------------------#
        P                   = calcpress(M, N, dx, dy, maxerror, ustar,\
                                            vstar, ul, ur, ub, ut, vl,\
                                            vr, vb, vt)
        #-----------------------------------------------------------------#
        # u and solution                                                  #
        #-----------------------------------------------------------------#
        for j in range(1, N+1):
            for i in range(1, M):
                u[j][i] = ustar[j][i] - (P[j][i+1]-P[j][i])/dx
        for j in range(1, N):
            for i in range(1, M+1):
                v[j][i] = vstar[j][i] - (P[j+1][i]-P[j][i])/dy
        #-----------------------------------------------------------------#
        # NEED TO COMMENT                                                 #
        #-----------------------------------------------------------------#
        if write_terms:
            for j in range(1, N+1):
                for i in range(1, M):
                    fterms.write("%16.5e" % Lu[j][i])
                fterms.write("\n")
            fterms.write("\n")
            for j in range(1, N):
                for i in range(1, M+1):
                    fterms.write("%16.5e" % Lv[j][i])
                fterms.write("\n")
            fterms.write("\n")
            for j in range(1, N+1):
                for i in range(1, M):
                    fterms.write("%16.5e" % Nx[j][i])
                fterms.write("\n")
            fterms.write("\n")
            for j in range(1, N):
                for i in range(1, M+1):
                    fterms.write("%16.5e" % Ny[j][i])
                fterms.write("\n")
            fterms.write("\n")
        #-----------------------------------------------------------------#
        # Writing pressure solution                                       #
        #-----------------------------------------------------------------#
        if write_P:
            for j in range(1, N+1):
                for i in range(1, M+1):
                    fP.write("%16.5e" % P[j][i])
                fP.write("\n")
            fP.write("\n")
        #-----------------------------------------------------------------#
        # output cell-centered (interpolated) u-velocity field values     #
        #-----------------------------------------------------------------#
        if write_vel_field or plot_vel_field:
            uinterp = np.zeros((N,M))
            vinterp = np.zeros((N,M))
            for j in range(1, N+1):
                # left
                velcent = (ul[j]+u[j][1])/2.
                uinterp[j-1][0] = velcent
                if write_vel_field: 
                    u_fvelfld.write("%16.5e" % velcent)
                # interior
                for i in range(2, M):
                    velcent = (u[j][i]+u[j][i-1])/2.
                    uinterp[j-1][i-1] = velcent
                    if write_vel_field: 
                        u_fvelfld.write("%16.5e" %velcent)
                # right
                velcent = (ur[j]+u[j][M-1])/2.
                uinterp[j-1][M-1] = velcent
                if write_vel_field:
                    u_fvelfld.write("%16.5e" %velcent)
                if write_vel_field:
                    u_fvelfld.write("\n")
            if write_vel_field: 
                u_fvelfld.write("\n")
        #-----------------------------------------------------------------#
        # output cell-centered (interpolated) v-velocity field values     #
        #-----------------------------------------------------------------#
            for i in range(1, M+1):
                velcent = (v[1][i]+vb[i])/2.
                vinterp[0][i-1] = velcent
                if write_vel_field:
                    v_fvelfld.write("%16.5e" % velcent)
            if write_vel_field: 
                v_fvelfld.write("\n")
            for j in range(2, N):
                for i in range(1, M+1):
                    velcent = (v[j][i]+v[j-1][i])/2.
                    vinterp[j-1][i-1] = velcent
                    if write_vel_field:
                        v_fvelfld.write("%16.5e" % velcent)
                if write_vel_field:
                    v_fvelfld.write("\n")
            for i in range(1, M+1):
                velcent = (v[N-1][i]+vt[i])/2.
                vinterp[N-1][i-1] = velcent
                if write_vel_field:
                    v_fvelfld.write("%16.5e" % velcent)
            if write_vel_field:
                v_fvelfld.write("\n")
            if write_vel_field:
                v_fvelfld.write("\n")
        #-----------------------------------------------------------------#
        # Plotting the velocity field                                     #
        #-----------------------------------------------------------------#
        if plot_vel_field:
            fig, ax = plt.subplots()
            q = ax.quiver(Xgrid,Ygrid,uinterp,vinterp)
            ax.set_aspect('equal')
            plt.savefig(media_path + 'plot-' +  str(nstep) + '.png')
            plt.close()
        #-----------------------------------------------------------------#
        # compute cell divergence                                         #
        #-----------------------------------------------------------------#
        #-----------------------------------------------------------------#
        # bottom section                                                  #
        #-----------------------------------------------------------------#
        # bottom left
        if write_divergence_field:
            fdiverg.write("%16.5e" % diverg)
        # across the bottom
        for i in range(2, M):
            diverg = (u[1][i]-u[1][i-1])/dx + (v[1][i]-vb[i])/dy
            if abs(diverg)>abs(maxdiverg):
                maxdiverg = diverg
            if write_divergence_field:
                fdiverg.write("%16.5e" % diverg)
        # bottom right
        diverg = (ur[1]-u[1][M-1])/dx + (v[1][M]-vb[M])/dy
        if abs(diverg)>abs(maxdiverg):
            maxdiverg = diverg
        if write_divergence_field:
            fdiverg.write("%16.5e" %(diverg))
        if write_divergence_field:
            fdiverg.write("\n")
        #-----------------------------------------------------------------#
        # Middle section                                                  #
        #-----------------------------------------------------------------#
        for j in range(2, N):
            # middle left
            diverg = (u[j][1]-ul[j])/dx + (v[j][1]-v[j-1][1])/dy
            if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
            if write_divergence_field: fdiverg.write("%16.5e" % diverg)
            # across the middle
            for i in range(2, M):
                diverg = (u[j][i]-u[j][i-1])/dx + (v[j][i]-v[j-1][i])/dy
                if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
                if write_divergence_field: fdiverg.write("%16.5e" % diverg)
            # middle right
            diverg = (ur[j]-u[j][M-1])/dx + (v[j][M]-v[j-1][M])/dy
            if abs(diverg)>abs(maxdiverg): 
                maxdiverg = diverg
            if write_divergence_field: 
                fdiverg.write("%16.5e" %(diverg))
            if write_divergence_field: 
                fdiverg.write("\n")
        #-----------------------------------------------------------------#
        # Top section                                                     #
        #-----------------------------------------------------------------#
        # top left
        diverg = (u[N][1]-ul[N])/dx + (vt[1]-v[N-1][1])/dy
        if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
        if write_divergence_field: fdiverg.write("%16.5e" % diverg)
        # across the top
        for i in range(2, M):
            diverg = (u[N][i]-u[N][i-1])/dx + (vt[i]-v[N-1][i])/dy
            if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
            if write_divergence_field: fdiverg.write("%16.5e" % diverg)
        # top right
        diverg = (ur[M]-u[N][M-1])/dx + (vt[M]-v[N-1][M])/dy
        if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
        if write_divergence_field: fdiverg.write("%16.5e" % diverg)
        if write_divergence_field: fdiverg.write("\n")
        if write_divergence_field: fdiverg.write("\n")

        print("%5i/%5i | Maximum velocity divergence: %16.5e" \
                    % (nstep, int(1.0000001*tmax/dt), maxdiverg))
    #---------------------------------------------------------------------#
    # Closing files                                                       #
    #---------------------------------------------------------------------#
    if write_vel_field:
        u_fvelfld.close()
        v_fvelfld.close()
    if write_terms:
        fterms.close()
    if write_P:
        fP.close()
    if write_divergence_field:
        fdiverg.close()

