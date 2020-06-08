#!/usr/bin/env python3
"""========================================================================
Purpose:
    py_471proj3.py: Spring 2005 MAE471 2D Lid-Driven Cavity Navier-Stokes
    Project Converted to Python.

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
import sys
from subprocess import call
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Time stepping subroutine                                                #
#-------------------------------------------------------------------------#
def ddt_prov(M, N, dx, dy, nu, u, v, ul, ur, ub, ut, vl, vr, vb, vt):

    """ Calculates time rate of change of provisional velocity field
        u*,v* (without pressure term) """
    #=====================================================================#
    # Domain variables                                                    #
    #=====================================================================#
    dx2 = dx*dx                         # spatial step size
    dy2 = dy*dy                         # spatial step size
    Lu  = np.zeros((N+1,M))             # linear u-velocity         i.e., laplacian(u)
    Lv  = np.zeros((N,M+1))             # linear v-velocity         i.e., laplacian(v)
    Nx  = np.zeros((N+1,M))             # non-linear u-velocity     i.e., u dot grad(u)
    Ny  = np.zeros((N,M+1))             # non-linear v-velocity     i.e., u dot grad(v)
    #=====================================================================#
    # Linear terms for u-velocity                                         #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Lower left wall                                                     #
    #---------------------------------------------------------------------#
    Lu[1][1]    = nu*(u[1][2]-2.*u[1][1]+ul[1])/dx2 
    Lu[1][1]    += nu*(u[2][1]-3.*u[1][1]+2.*ub[1])/dy2     
    #---------------------------------------------------------------------#
    # Bottom wall                                                         #
    #---------------------------------------------------------------------#
    for i in range(2, M-1):
        Lu[1][i]    = nu*(u[1][i+1]-2.*u[1][i]+u[1][i-1])/dx2
        Lu[1][i]    += nu*(u[2][i]-3.*u[1][i]+2.*ub[i])/dy2
    #---------------------------------------------------------------------#
    # Lower right wall                                                    #
    #---------------------------------------------------------------------#
    Lu[1][M-1]  = nu*(ur[1]-2.*u[1][M-1]+u[1][M-2])/dx2
    Lu[1][M-1]  += nu*(u[2][M-1]-3.*u[1][M-1]+2.*ub[M-1])/dy2
    #---------------------------------------------------------------------#
    # Interior domain                                                     #
    #---------------------------------------------------------------------#
    for j in range(2, N):
        #-----------------------------------------------------------------#
        # Left wall                                                       #
        #-----------------------------------------------------------------#
        Lu[j][1]    = nu*(u[j][2]-2.*u[j][1]+ul[j])/dx2
        Lu[j][1]    += nu*(u[j+1][1]-2.*u[j][1]+u[j-1][1])/dy2
        #-----------------------------------------------------------------#
        # Interior                                                        #
        #-----------------------------------------------------------------#
        for i in range(2, M-1):
            Lu[j][i]    = nu*(u[j][i+1]-2.*u[j][i]+u[j][i-1])/dx2
            Lu[j][i]    += nu*(u[j+1][i]-2.*u[j][i]+u[j-1][i])/dy2
        #-----------------------------------------------------------------#
        # Right wall                                                      #
        #-----------------------------------------------------------------#
        Lu[j][M-1]  = nu*(ur[j]-2.*u[j][M-1]+u[j][M-2])/dx2
        Lu[j][M-1]  += nu*(u[j+1][M-1]-2.*u[j][M-1]+u[j-1][M-1])/dy2
    #---------------------------------------------------------------------#
    # Upper left wall                                                     #
    #---------------------------------------------------------------------#
    Lu[N][1] = nu*(u[N][2]-2.*u[N][1]+ul[N])/dx2
    Lu[N][1] += nu*(2.*ut[1]-3.*u[N][1]+u[N-1][1])/dy2
    #---------------------------------------------------------------------#
    # Top wall                                                            #
    #---------------------------------------------------------------------#
    for i in range(2, M-1):
        Lu[N][i]    = nu*(u[N][i+1]-2.*u[N][i]+u[N][i-1])/dx2
        Lu[N][i]    += nu*(2.*ut[i]-3.*u[N][i]+u[N-1][i])/dy2
    #---------------------------------------------------------------------#
    # Upper right wall                                                    #
    #---------------------------------------------------------------------#
    Lu[N][M-1]  = nu*(ur[N]-2.*u[N][M-1]+u[N][M-2])/dx2
    Lu[N][M-1]  += nu*(2.*ut[M-1]-3.*u[N][M-1]+u[N-1][M-1])/dy2
    #=====================================================================#
    # Linear terms for v-velocity                                         #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Lower left                                                          #
    #---------------------------------------------------------------------#
    Lv[1][1] = nu*(v[1][2]-3.*v[1][1]+2.*vl[1])/dx2
    Lv[1][1] += nu*(v[2][1]-2.*v[1][1]+vb[1])/dy2
    #---------------------------------------------------------------------#
    # Bottom wall                                                         #
    #---------------------------------------------------------------------#
    for i in range(2, M):
        Lv[1][i] = nu*(v[1][i+1]-2.*v[1][i]+v[1][i-1])/dx2
        Lv[1][i] += nu*(v[2][i]-2.*v[1][i]+vb[i])/dy2
    #---------------------------------------------------------------------#
    # Lower right                                                         #
    #---------------------------------------------------------------------#
    Lv[1][M] = nu*(2.0*vr[1]-3.*v[1][M]+v[1][M-1])/dx2
    Lv[1][M] += nu*(v[2][M]-2.*v[1][M]+vb[M])/dy2
    #---------------------------------------------------------------------#
    # Interior points                                                     #
    #---------------------------------------------------------------------#
    for j in range(2, N-1):
        #-----------------------------------------------------------------#
        # Left wall                                                       #
        #-----------------------------------------------------------------#
        Lv[j][1] = nu*(v[j][2]-3.*v[j][1]+2.*vl[j])/dx2
        Lv[j][1] += nu*(v[j+1][1]-2.*v[j][1]+v[j-1][1])/dy2
        #-----------------------------------------------------------------#
        # interior points                                                 #
        #-----------------------------------------------------------------#
        for i in range(2, M):
            Lv[j][i] = nu*(v[j][i+1]-2.*v[j][i]+v[j][i-1])/dx2
            Lv[j][i] += nu*(v[j+1][i]-2.*v[j][i]+v[j-1][i])/dy2
        #-----------------------------------------------------------------#
        # Right wall                                                      #
        #-----------------------------------------------------------------#
        Lv[j][M] = nu*(2.*vr[j]-3.*v[j][M]+v[j][M-1])/dx2
        Lv[j][M] += nu*(v[j+1][M]-2.*v[j][M]+v[j-1][M])/dy2
    #---------------------------------------------------------------------#
    # Upper left wall                                                     #
    #---------------------------------------------------------------------#
    Lv[N-1][1] = nu*(v[N-1][2]-3.*v[N-1][1]+2.*vl[N-1])/dx2
    Lv[N-1][1] += nu*(vt[1]-2.*v[N-1][1]+v[N-2][1])/dy2
    #---------------------------------------------------------------------#
    # Top wall                                                            #
    #---------------------------------------------------------------------#
    for i in range(2, M):
        Lv[N-1][i] = nu*(v[N-1][i+1]-2.*v[N-1][i]+v[N-1][i-1])/dx2
        Lv[N-1][i] += nu*(vt[i]-2.*v[N-1][i]+v[N-2][i])/dy2
    #---------------------------------------------------------------------#
    # Upper right wall                                                    #
    #---------------------------------------------------------------------#
    Lv[N-1][M] = nu*(2.*vr[N-1]-3.*v[N-1][M]+v[N-1][M-1])/dx2
    Lv[N-1][M] += nu*(vt[M]-2.*v[N-1][M]+v[N-2][M])/dy2
    #=====================================================================#
    # Non-linear terms for u-velocity                                     #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Looping over domain                                                 #
    #---------------------------------------------------------------------#
    for j in range(1, N+1):
        #-----------------------------------------------------------------#
        # Left wall                                                       #
        #-----------------------------------------------------------------#
        Nx[j][1] = (0.25*(u[j][2]+u[j][1])**2 - 0.25*(u[j][1]+ul[j])**2)/dx
        #-----------------------------------------------------------------#
        # Interior domain                                                 #
        #-----------------------------------------------------------------#
        for i in range(2, M-1):
            Nx[j][i] = (0.25*(u[j][i+1] + u[j][i])**2 - 0.25*(u[j][i] +\
                        u[j][i-1])**2)/dx
        #-----------------------------------------------------------------#
        # Right wall                                                      #
        #-----------------------------------------------------------------#
        Nx[j][M-1] = (0.25*(ur[j]+u[j][M-1])**2 - 0.25*(u[j][M-1]+ \
                        u[j][M-2])**2)/dx
    #---------------------------------------------------------------------#
    # Bottom wall                                                         #
    #---------------------------------------------------------------------#
    for i in range(1, M):
        Nx[1][i] += (0.25*(u[2][i]+u[1][i])*(v[1][i+1]+v[1][i]) - \
                    0.5*(ub[i])*(vb[i]+vb[i+1]))/dy
    #---------------------------------------------------------------------#
    # Interior domain                                                     #
    #---------------------------------------------------------------------#
    for j in range(2, N):
        for i in range(1, M):
            Nx[j][i] += (0.25*(u[j+1][i]+u[j][i])*(v[j][i+1]+v[j][i]) -\
                        0.25*(u[j][i]+u[j-1][i])*(v[j-1][i]+v[j-1][i+1]))/dy
    #---------------------------------------------------------------------#
    # Top wall                                                            #
    #---------------------------------------------------------------------#
    for i in range(1, M):
        Nx[N][i] += (0.5*(ut[i])*(vt[i+1]+vt[i]) - \
                    0.25*(u[N][i]+u[N-1][i])*(v[N-1][i]+v[N-1][i+1]))/dy
    #=====================================================================#
    # Non-linear terms for v-velocity                                     #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Left wall                                                           #
    #---------------------------------------------------------------------#
    for i in range(1, M+1):
        Ny[1][i] = (0.25*(v[2][i]+v[1][i])**2 - \
                        0.25*(v[1][i]+vb[i])**2)/dy
    #---------------------------------------------------------------------#
    # Interior points dv/dy                                               #
    #---------------------------------------------------------------------#
    for j in range(2, N-1):
        for i in range(1, M+1):
            Ny[j][i] = (0.25*(v[j+1][i] + v[j][i])**2 - \
                            0.25*(v[j][i] + v[j-1][i])**2)/dy
    #---------------------------------------------------------------------#
    # Top wall                                                            # 
    #---------------------------------------------------------------------#
    for i in range(1, M+1):
        Ny[N-1][i] = (0.25*(vt[i]+v[N-1][i])**2 - \
                        0.25*(v[N-1][i]+v[N-2][i])**2)/dy
    #---------------------------------------------------------------------#
    # Interior points dv/dx                                               # 
    #---------------------------------------------------------------------#
    for j in range(1, N):
        #-----------------------------------------------------------------#
        # Left domain                                                     # 
        #-----------------------------------------------------------------#
        Ny[j][1] += (0.25*(u[j][1]+u[j+1][1])*(v[j][2]+v[j][1]) -\
                        0.5*(ul[j]+ul[j+1])*(vl[j]))/dx
        #-----------------------------------------------------------------#
        # Interior domain                                                 # 
        #-----------------------------------------------------------------#
        for i in range(2, M):
            Ny[j][i] += (0.25*(u[j][i]+u[j+1][i])*(v[j][i+1]+v[j][i]) -\
                        0.25*(u[j][i-1]+u[j+1][i-1])*(v[j][i]+v[j][i-1]))/dx
        #-----------------------------------------------------------------#
        # Right wall domain                                               # 
        #-----------------------------------------------------------------#
        Ny[j][M] += (0.5*(ur[j]+ur[j+1])*(vr[j]) - \
                        0.25*(u[j][M-1]+u[j+1][M-1])*(v[j][M]+v[j][M-1]))/dx

    return [Lu, Lv, Nx, Ny]
#-------------------------------------------------------------------------#
# Pressure implementation                                                 #
#-------------------------------------------------------------------------#
def calcpress(M, N, dx, dy, maxerror, ustar, vstar, ul, ur, ub, ut, vl,\
                vr, vb, vt):
    
    """ Calculating the pressure; Gauss-Seidel part """
    #=====================================================================#
    # Domain variables                                                    #
    #=====================================================================#
    gsiter      = 0
    normnew     = 0.
    diffnorm    = 0.
    normerror   = 2.*maxerror     # normerror is relative error in norm of diagonal
    rhs         = 0.                         
    dx2         = dx*dx
    dy2         = dy*dy
    #=====================================================================#
    # Pressure variables                                                  #
    #=====================================================================#
    P       = np.zeros((N+1,M+1))
    Pold    = np.zeros((N+1,M+1))
    #=====================================================================#
    # Gauss-Seidel iteration                                              #
    #=====================================================================#
    while((normerror > maxerror and gsiter < 5000) or (gsiter < 30)):
        #-----------------------------------------------------------------#
        # Norm variables                                                  #
        #-----------------------------------------------------------------#
        normnew     = 0.
        diffnorm    = 0.        
        gsiter      += 1
        #-----------------------------------------------------------------#
        # Lower left                                                      #
        #-----------------------------------------------------------------#
        rhs     = (ustar[1][1]-ul[1])/dx + (vstar[1][1]-vb[1])/dy
        P[1][1] = 1./(1./dx2+1./dy2)*(P[1][2]/dx2+P[2][1]/dy2-rhs)
        #-----------------------------------------------------------------#
        # lower                                                           #
        #-----------------------------------------------------------------#
        for i in range(2, M):
            rhs     = (ustar[1][i]-ustar[1][i-1])/dx + (vstar[1][i]-vb[i])/dy
            P[1][i] = 1./(2./dx2+1./dy2)*((P[1][i+1]+P[1][i-1])/dx2+P[2][i]/dy2-rhs)
        #-----------------------------------------------------------------#
        # lower right                                                     #
        #-----------------------------------------------------------------#
        rhs     = (ur[1]-ustar[1][M-1])/dx + (vstar[1][M]-vb[M])/dy
        P[1][M] = 1./(1./dx2+1./dy2)*(P[1][M-1]/dx2+P[2][M]/dy2-rhs)
        #-----------------------------------------------------------------#
        # left                                                            #
        #-----------------------------------------------------------------#
        for j in range(2, N):
            rhs     = (ustar[j][1]-ul[j])/dx + (vstar[j][1]-vstar[j-1][1])/dy
            P[j][1] = 1./(1./dx2+2./dy2)*(P[j][2]/dx2+(P[j+1][1]+P[j-1][1])/dy2-rhs)
            #-------------------------------------------------------------#
            # away from the boundaries                                    #
            #-------------------------------------------------------------#
            for i in range(2, M):
                rhs     = (ustar[j][i]-ustar[j][i-1])/dx + (vstar[j][i]-vstar[j-1][i])/dy
                P[j][i] = 1./(2./dx2+2./dy2)*((P[j][i+1]+P[j][i-1])/dx2+(P[j+1][i]+P[j-1][i])/dy2-rhs)
            #-------------------------------------------------------------#
            # right                                                       #
            #-------------------------------------------------------------#
            rhs     = (ur[j]-ustar[j][M-1])/dx + (vstar[j][M]-vstar[j-1][M])/dy
            P[j][M] = 1./(1./dx2+2./dy2)*(P[j][M-1]/dx2+(P[j+1][M]+P[j-1][M])/dy2-rhs)
        #-----------------------------------------------------------------#
        # upper left                                                      #
        #-----------------------------------------------------------------#
        rhs     = (ustar[N][1]-ul[N])/dx + (vt[1]-vstar[N-1][1])/dy
        P[N][1] = 1./(1./dx2+1./dy2)*(P[N][2]/dx2+P[N-1][1]/dy2-rhs)
        #-----------------------------------------------------------------#
        # upper                                                           #
        #-----------------------------------------------------------------#
        for i in range(2, M):
            rhs     = (ustar[N][i]-ustar[N][i-1])/dx + (vt[i]-vstar[N-1][i])/dy
            P[N][i] = 1./(2./dx2+1./dy2)*((P[N][i+1]+P[N][i-1])/dx2+P[N-1][i]/dy2-rhs)
        #-----------------------------------------------------------------#
        # upper right                                                     #
        #-----------------------------------------------------------------#
        rhs     = (ur[N]-ustar[N][M-1])/dx + (vt[M]-vstar[N-1][M])/dy
        P[N][M] = 1./(1./dx2+1./dy2)*(P[N][M-1]/dx2+P[N-1][M]/dy2-rhs)
        #-----------------------------------------------------------------#
        # updating norms                                                  #
        #-----------------------------------------------------------------#
        for j in range(1,N+1):
            for i in range(1,M+1):
                normnew     += P[j][i]**2
                diffnorm    += (P[j][i]-Pold[j][i])**2
        #-----------------------------------------------------------------#
        # new norm                                                        #
        #-----------------------------------------------------------------#
        normerror = diffnorm/normnew
        #-----------------------------------------------------------------#
        # storing old pressure                                            #
        #-----------------------------------------------------------------#
        Pold = np.copy(P)    
    #---------------------------------------------------------------------#
    # Print statement                                                     #
    #---------------------------------------------------------------------#
    print("%16.5e: %5i iterations" % (normerror, gsiter))

    return P
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    #  path variables                                                     #
    #---------------------------------------------------------------------#
    call(['clear'])
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    M   = 32        # nx
    N   = 32       # ny
    lx  = 1.0       # x length 
    ly  = 1.0       # y length
    dx  = lx/M      # x spatial step size
    dy  = ly/N      # y spatial step size 
    #---------------------------------------------------------------------#
    # Writing flags                                                       #
    #---------------------------------------------------------------------#
    write_vel_field             = True
    write_terms                 = False
    write_P                     = False
    write_divergence_field      = True
    plot_vel_field              = False
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    u       = np.zeros((N+1,M))
    v       = np.zeros((N,M+1))
    ustar   = np.zeros((N+1,M))
    vstar   = np.zeros((N,M+1))
    ul      = np.zeros((N+1))
    ur      = np.zeros((N+1))
    vl      = np.zeros((N+1))
    vr      = np.zeros((N+1))
    ub      = np.zeros((M+1))
    ut      = np.zeros((M+1))
    vb      = np.zeros((M+1))
    vt      = np.zeros((M+1))
    
    
    #maxerror = 1.e-12
    maxerror = 1.e-8
    
    diverg = 0.
    maxdiverg = 0.
    
    ## set up BC's, velocity components at boundaries are all set to 0 until now -- lid-driven cavity
    #for i in range(1,M+1):
    #    ut[i] = 1.
    
    # set up BC's, velocity components at boundaries are all set to 0 until now -- unidirectional flow
    for i in range(1,M+1):
        ut[i] = 1.
        ub[i] = 0.0
    
    for j in range(1,N+1):
        ur[j] = 0.0
        ul[j] = 0.0
    
    # initial condition has already been set to 0 everywhere
    
    nu = 0.01    # viscosity
    
    dt = 0.25*dx*dx/nu/2.
    print(dt)
    tmax = 2.5                         # final time
    
    # file output
    if write_vel_field: fvelfld = open('vel_field.txt', 'w')
    if write_terms: fterms = open('terms.txt', 'w')
    if write_P: fP = open('P.txt', 'w')
    if write_divergence_field: fdiverg = open('divg_field.txt', 'w')
    
    Xgrid,Ygrid = np.meshgrid(np.linspace(dx/2,lx-dx/2,M),np.linspace(dy/2,ly-dy/2,N))
    
    t = 0.0
    # main time stepping loop
    for nstep in range(1, int(1.0000001*tmax/dt)+1):
        t += dt
        [Lu, Lv, Nx, Ny] = ddt_prov(M, N, dx, dy, nu, u, v, ul, ur, ub, ut, vl, vr, vb, vt)
        
        
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
    
        for j in range(1, N+1):
            for i in range(1, M):
                ustar[j][i] = u[j][i] + dt*(-Nx[j][i] + Lu[j][i])
    
        for j in range(1, N):
            for i in range(1, M+1):
                vstar[j][i] = v[j][i] + dt*(-Ny[j][i] + Lv[j][i])
    
        P = calcpress(M, N, dx, dy, maxerror, ustar, vstar, ul, ur, ub, ut, vl, vr, vb, vt)
    
        if write_P:
            for j in range(1, N+1):
                for i in range(1, M+1):
                    fP.write("%16.5e" % P[j][i])
                fP.write("\n")
            fP.write("\n")
    
        for j in range(1, N+1):
            for i in range(1, M):
                u[j][i] = ustar[j][i] - (P[j][i+1]-P[j][i])/dx
    
        for j in range(1, N):
            for i in range(1, M+1):
                v[j][i] = vstar[j][i] - (P[j+1][i]-P[j][i])/dy
    
        # output cell-centered (interpolated) field values
        if write_vel_field or plot_vel_field:
            uinterp = np.zeros((N,M))
            vinterp = np.zeros((N,M))
            for j in range(1, N+1):
                # left
                velcent = (ul[j]+u[j][1])/2.
                uinterp[j-1][0] = velcent
                if write_vel_field: fvelfld.write("%16.5e" % velcent)
                for i in range(2, M):
                    velcent = (u[j][i]+u[j][i-1])/2.
                    uinterp[j-1][i-1] = velcent
                    if write_vel_field: fvelfld.write("%16.5e" % velcent)
                # right
                velcent = (ur[j]+u[j][M-1])/2.
                uinterp[j-1][M-1] = velcent
                if write_vel_field: fvelfld.write("%16.5e" % velcent)
                if write_vel_field: fvelfld.write("\n")
            if write_vel_field: fvelfld.write("\n")
    
            for i in range(1, M+1):
                velcent = (v[1][i]+vb[i])/2.
                vinterp[0][i-1] = velcent
                if write_vel_field: fvelfld.write("%16.5e" % velcent)
            if write_vel_field: fvelfld.write("\n")
            for j in range(2, N):
                for i in range(1, M+1):
                    velcent = (v[j][i]+v[j-1][i])/2.
                    vinterp[j-1][i-1] = velcent
                    if write_vel_field: fvelfld.write("%16.5e" % velcent)
                if write_vel_field: fvelfld.write("\n")
            for i in range(1, M+1):
                velcent = (v[N-1][i]+vt[i])/2.
                vinterp[N-1][i-1] = velcent
                if write_vel_field: fvelfld.write("%16.5e" % velcent)
            if write_vel_field: fvelfld.write("\n")
            if write_vel_field: fvelfld.write("\n")
    
        if plot_vel_field:
            fig, ax = plt.subplots()
            q = ax.quiver(Xgrid,Ygrid,uinterp,vinterp)
            ax.set_aspect('equal')
            plt.show()
    
        # compute cell divergence
    
        # bottom left
        diverg = (u[1][1]-ul[1])/dx + (v[1][1]-vb[1])/dy
        maxdiverg = diverg
        if write_divergence_field: fdiverg.write("%16.5e" % diverg)
        # across the bottom
        for i in range(2, M):
            diverg = (u[1][i]-u[1][i-1])/dx + (v[1][i]-vb[i])/dy
            if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
            if write_divergence_field: fdiverg.write("%16.5e" % diverg)
        # bottom right
        diverg = (ur[1]-u[1][M-1])/dx + (v[1][M]-vb[M])/dy
        if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
        if write_divergence_field: fdiverg.write("%16.5e" % diverg)
        if write_divergence_field: fdiverg.write("\n")
    
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
            if abs(diverg)>abs(maxdiverg): maxdiverg = diverg
            if write_divergence_field: fdiverg.write("%16.5e" % diverg)
            if write_divergence_field: fdiverg.write("\n")
    
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
    
        print("%5i/%5i | Maximum velocity divergence: %16.5e" % (nstep, int(1.0000001*tmax/dt), maxdiverg))
    
        if nstep ==6:
            print(t)
            for j in range(0, N+1):
                for i in range(0,M+1):
                    print('%i %i %.5f'  %(j , i, P[j][i]))
            #sys.exit(0)
    if write_vel_field: fvelfld.close()
    if write_terms: fterms.close()
    if write_P: fP.close()
    if write_divergence_field: fdiverg.close()
    