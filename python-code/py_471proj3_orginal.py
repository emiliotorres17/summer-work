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
import os
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Time stepping subroutine                                                #
#-------------------------------------------------------------------------#
def ddt_prov(Min, Nin, DX, DY, NU, U, V, Ul, Ur, Ub, Ut, Vl, Vr, Vb, Vt):

    """ Calculates time rate of change of provisional velocity field
        U*,V* (without pressure term) """
    #=====================================================================#
    # Domain variables                                                    #
    #=====================================================================#
    dx2 = DX*DX                             # spatial step size
    dy2 = DY*DY                             # spatial step size
    LU  = np.zeros((Nin+1,Min))             # linear U-velocity         I.e., laplacian(U)
    LV  = np.zeros((Nin,Min+1))             # linear V-velocity         I.e., laplacian(V)
    NX  = np.zeros((Nin+1,Min))             # non-linear U-velocity     I.e., U dot grad(U)
    NY  = np.zeros((Nin,Min+1))             # non-linear V-velocity     I.e., U dot grad(V)
    #=====================================================================#
    # Linear terms for U-velocity                                         #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Lower left wall                                                     #
    #---------------------------------------------------------------------#
    LU[1][1]    = NU*(U[1][2]-2.*U[1][1]+Ul[1])/dx2
    LU[1][1]    += NU*(U[2][1]-3.*U[1][1]+2.*Ub[1])/dy2
    #---------------------------------------------------------------------#
    # Bottom wall                                                         #
    #---------------------------------------------------------------------#
    for I in range(2, Min-1):
        LU[1][I]    = NU*(U[1][I+1]-2.*U[1][I]+U[1][I-1])/dx2
        LU[1][I]    += NU*(U[2][I]-3.*U[1][I]+2.*Ub[I])/dy2
    #---------------------------------------------------------------------#
    # Lower right wall                                                    #
    #---------------------------------------------------------------------#
    LU[1][Min-1]  = NU*(Ur[1]-2.*U[1][Min-1]+U[1][Min-2])/dx2
    LU[1][Min-1]  += NU*(U[2][Min-1]-3.*U[1][Min-1]+2.*Ub[Min-1])/dy2
    #---------------------------------------------------------------------#
    # Interior domain                                                     #
    #---------------------------------------------------------------------#
    for J in range(2, Nin):
        #-----------------------------------------------------------------#
        # Left wall                                                       #
        #-----------------------------------------------------------------#
        LU[J][1]    = NU*(U[J][2]-2.*U[J][1]+Ul[J])/dx2
        LU[J][1]    += NU*(U[J+1][1]-2.*U[J][1]+U[J-1][1])/dy2
        #-----------------------------------------------------------------#
        # Interior                                                        #
        #-----------------------------------------------------------------#
        for I in range(2, Min-1):
            LU[J][I]    = NU*(U[J][I+1]-2.*U[J][I]+U[J][I-1])/dx2
            LU[J][I]    += NU*(U[J+1][I]-2.*U[J][I]+U[J-1][I])/dy2
        #-----------------------------------------------------------------#
        # Right wall                                                      #
        #-----------------------------------------------------------------#
        LU[J][Min-1]  = NU*(Ur[J]-2.*U[J][Min-1]+U[J][Min-2])/dx2
        LU[J][Min-1]  += NU*(U[J+1][Min-1]-2.*U[J][Min-1]+U[J-1][Min-1])/dy2
    #---------------------------------------------------------------------#
    # Upper left wall                                                     #
    #---------------------------------------------------------------------#
    LU[Nin][1] = NU*(U[Nin][2]-2.*U[Nin][1]+Ul[Nin])/dx2
    LU[Nin][1] += NU*(2.*Ut[1]-3.*U[Nin][1]+U[Nin-1][1])/dy2
    #---------------------------------------------------------------------#
    # Top wall                                                            #
    #---------------------------------------------------------------------#
    for I in range(2, Min-1):
        LU[Nin][I]    = NU*(U[Nin][I+1]-2.*U[Nin][I]+U[Nin][I-1])/dx2
        LU[Nin][I]    += NU*(2.*Ut[I]-3.*U[Nin][I]+U[Nin-1][I])/dy2
    #---------------------------------------------------------------------#
    # Upper right wall                                                    #
    #---------------------------------------------------------------------#
    LU[Nin][Min-1]  = NU*(Ur[Nin]-2.*U[Nin][Min-1]+U[Nin][Min-2])/dx2
    LU[Nin][Min-1]  += NU*(2.*Ut[Min-1]-3.*U[Nin][Min-1]+U[Nin-1][Min-1])/dy2
    #=====================================================================#
    # Linear terms for V-velocity                                         #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Lower left                                                          #
    #---------------------------------------------------------------------#
    LV[1][1] = NU*(V[1][2]-3.*V[1][1]+2.*Vl[1])/dx2
    LV[1][1] += NU*(V[2][1]-2.*V[1][1]+Vb[1])/dy2
    #---------------------------------------------------------------------#
    # Bottom wall                                                         #
    #---------------------------------------------------------------------#
    for I in range(2, Min):
        LV[1][I] = NU*(V[1][I+1]-2.*V[1][I]+V[1][I-1])/dx2
        LV[1][I] += NU*(V[2][I]-2.*V[1][I]+Vb[I])/dy2
    #---------------------------------------------------------------------#
    # Lower right                                                         #
    #---------------------------------------------------------------------#
    LV[1][Min] = NU*(2.0*Vr[1]-3.*V[1][Min]+V[1][Min-1])/dx2
    LV[1][Min] += NU*(V[2][Min]-2.*V[1][Min]+Vb[Min])/dy2
    #---------------------------------------------------------------------#
    # Interior points                                                     #
    #---------------------------------------------------------------------#
    for J in range(2, Nin-1):
        #-----------------------------------------------------------------#
        # Left wall                                                       #
        #-----------------------------------------------------------------#
        LV[J][1] = NU*(V[J][2]-3.*V[J][1]+2.*Vl[J])/dx2
        LV[J][1] += NU*(V[J+1][1]-2.*V[J][1]+V[J-1][1])/dy2
        #-----------------------------------------------------------------#
        # interior points                                                 #
        #-----------------------------------------------------------------#
        for I in range(2, Min):
            LV[J][I] = NU*(V[J][I+1]-2.*V[J][I]+V[J][I-1])/dx2
            LV[J][I] += NU*(V[J+1][I]-2.*V[J][I]+V[J-1][I])/dy2
        #-----------------------------------------------------------------#
        # Right wall                                                      #
        #-----------------------------------------------------------------#
        LV[J][Min] = NU*(2.*Vr[J]-3.*V[J][Min]+V[J][Min-1])/dx2
        LV[J][Min] += NU*(V[J+1][Min]-2.*V[J][Min]+V[J-1][Min])/dy2
    #---------------------------------------------------------------------#
    # Upper left wall                                                     #
    #---------------------------------------------------------------------#
    LV[Nin-1][1] = NU*(V[Nin-1][2]-3.*V[Nin-1][1]+2.*Vl[Nin-1])/dx2
    LV[Nin-1][1] += NU*(Vt[1]-2.*V[Nin-1][1]+V[Nin-2][1])/dy2
    #---------------------------------------------------------------------#
    # Top wall                                                            #
    #---------------------------------------------------------------------#
    for I in range(2, Min):
        LV[Nin-1][I] = NU*(V[Nin-1][I+1]-2.*V[Nin-1][I]+V[Nin-1][I-1])/dx2
        LV[Nin-1][I] += NU*(Vt[I]-2.*V[Nin-1][I]+V[Nin-2][I])/dy2
    #---------------------------------------------------------------------#
    # Upper right wall                                                    #
    #---------------------------------------------------------------------#
    LV[Nin-1][Min] = NU*(2.*Vr[Nin-1]-3.*V[Nin-1][Min]+V[Nin-1][Min-1])/dx2
    LV[Nin-1][Min] += NU*(Vt[Min]-2.*V[Nin-1][Min]+V[Nin-2][Min])/dy2
    #=====================================================================#
    # Non-linear terms for U-velocity                                     #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Looping over domain                                                 #
    #---------------------------------------------------------------------#
    for J in range(1, Nin+1):
        #-----------------------------------------------------------------#
        # Left wall                                                       #
        #-----------------------------------------------------------------#
        NX[J][1] = (0.25*(U[J][2]+U[J][1])**2 - 0.25*(U[J][1]+Ul[J])**2)/DX
        #-----------------------------------------------------------------#
        # Interior domain                                                 #
        #-----------------------------------------------------------------#
        for I in range(2, Min-1):
            NX[J][I] = (0.25*(U[J][I+1] + U[J][I])**2 - 0.25*(U[J][I] +\
                        U[J][I-1])**2)/DX
        #-----------------------------------------------------------------#
        # Right wall                                                      #
        #-----------------------------------------------------------------#
        NX[J][Min-1] = (0.25*(Ur[J]+U[J][Min-1])**2 - 0.25*(U[J][Min-1]+ \
                        U[J][Min-2])**2)/DX
    #---------------------------------------------------------------------#
    # Bottom wall                                                         #
    #---------------------------------------------------------------------#
    for I in range(1, Min):
        NX[1][I] += (0.25*(U[2][I]+U[1][I])*(V[1][I+1]+V[1][I]) - \
                    0.5*(Ub[I])*(Vb[I]+Vb[I+1]))/DY
    #---------------------------------------------------------------------#
    # Interior domain                                                     #
    #---------------------------------------------------------------------#
    for J in range(2, Nin):
        for I in range(1, Min):
            NX[J][I] += (0.25*(U[J+1][I]+U[J][I])*(V[J][I+1]+V[J][I]) -\
                        0.25*(U[J][I]+U[J-1][I])*(V[J-1][I]+V[J-1][I+1]))/DY
    #---------------------------------------------------------------------#
    # Top wall                                                            #
    #---------------------------------------------------------------------#
    for I in range(1, Min):
        NX[Nin][I] += (0.5*(Ut[I])*(Vt[I+1]+Vt[I]) - \
                    0.25*(U[Nin][I]+U[Nin-1][I])*(V[Nin-1][I]+V[Nin-1][I+1]))/DY
    #=====================================================================#
    # Non-linear terms for V-velocity                                     #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Left wall                                                           #
    #---------------------------------------------------------------------#
    for I in range(1, Min+1):
        NY[1][I] = (0.25*(V[2][I]+V[1][I])**2 - \
                        0.25*(V[1][I]+Vb[I])**2)/DY
    #---------------------------------------------------------------------#
    # Interior points dv/DY                                               #
    #---------------------------------------------------------------------#
    for J in range(2, Nin-1):
        for I in range(1, Min+1):
            NY[J][I] = (0.25*(V[J+1][I] + V[J][I])**2 - \
                            0.25*(V[J][I] + V[J-1][I])**2)/DY
    #---------------------------------------------------------------------#
    # Top wall                                                            #
    #---------------------------------------------------------------------#
    for I in range(1, Min+1):
        NY[Nin-1][I] = (0.25*(Vt[I]+V[Nin-1][I])**2 - \
                        0.25*(V[Nin-1][I]+V[Nin-2][I])**2)/DY
    #---------------------------------------------------------------------#
    # Interior points dv/DX                                               #
    #---------------------------------------------------------------------#
    for J in range(1, Nin):
        #-----------------------------------------------------------------#
        # Left domain                                                     #
        #-----------------------------------------------------------------#
        NY[J][1] += (0.25*(U[J][1]+U[J+1][1])*(V[J][2]+V[J][1]) -\
                        0.5*(Ul[J]+Ul[J+1])*(Vl[J]))/DX
        #-----------------------------------------------------------------#
        # Interior domain                                                 #
        #-----------------------------------------------------------------#
        for I in range(2, Min):
            NY[J][I] += (0.25*(U[J][I]+U[J+1][I])*(V[J][I+1]+V[J][I]) -\
                        0.25*(U[J][I-1]+U[J+1][I-1])*(V[J][I]+V[J][I-1]))/DX
        #-----------------------------------------------------------------#
        # Right wall domain                                               #
        #-----------------------------------------------------------------#
        NY[J][Min] += (0.5*(Ur[J]+Ur[J+1])*(Vr[J]) - \
                        0.25*(U[J][Min-1]+U[J+1][Min-1])*(V[J][Min]+V[J][Min-1]))/DX

    return [LU, LV, NX, NY]
#-------------------------------------------------------------------------#
# Pressure implementation                                                 #
#-------------------------------------------------------------------------#
def calcpress(Min, Nin, DX, DY, Maxerror, Ustar, Vstar, Ul, Ur, Vb, Vt):

    """ Calculating the pressure; Gauss-Seidel part """
    #=====================================================================#
    # Domain variables                                                    #
    #=====================================================================#
    gsiter      = 0
    normnew     = 0.
    diffnorm    = 0.
    normerror   = 2.*Maxerror     # normerror is relative error in norm of diagonal
    rhs         = 0.
    dx2         = DX*DX
    dy2         = DY*DY
    #=====================================================================#
    # Pressure variables                                                  #
    #=====================================================================#
    Press   = np.zeros((Nin+1,Min+1))
    Pold    = np.zeros((Nin+1,Min+1))
    #=====================================================================#
    # Gauss-Seidel iteration                                              #
    #=====================================================================#
    while((normerror > Maxerror and gsiter < 5000) or (gsiter < 30)):
        #-----------------------------------------------------------------#
        # Norm variables                                                  #
        #-----------------------------------------------------------------#
        normnew     = 0.
        diffnorm    = 0.
        gsiter      += 1
        #-----------------------------------------------------------------#
        # Lower left                                                      #
        #-----------------------------------------------------------------#
        rhs     = (Ustar[1][1]-Ul[1])/DX + (Vstar[1][1]-Vb[1])/DY
        Press[1][1] = 1./(1./dx2+1./dy2)*(Press[1][2]/dx2+Press[2][1]/dy2-rhs)
        #-----------------------------------------------------------------#
        # lower                                                           #
        #-----------------------------------------------------------------#
        for I in range(2, Min):
            rhs     = (Ustar[1][I]-Ustar[1][I-1])/DX + (Vstar[1][I]-Vb[I])/DY
            Press[1][I] = 1./(2./dx2+1./dy2)*((Press[1][I+1]+Press[1][I-1])/dx2+Press[2][I]/dy2-rhs)
        #-----------------------------------------------------------------#
        # lower right                                                     #
        #-----------------------------------------------------------------#
        rhs     = (Ur[1]-Ustar[1][Min-1])/DX + (Vstar[1][Min]-Vb[Min])/DY
        Press[1][Min] = 1./(1./dx2+1./dy2)*(Press[1][Min-1]/dx2+Press[2][Min]/dy2-rhs)
        #-----------------------------------------------------------------#
        # left                                                            #
        #-----------------------------------------------------------------#
        for J in range(2, Nin):
            rhs     = (Ustar[J][1]-Ul[J])/DX + (Vstar[J][1]-Vstar[J-1][1])/DY
            Press[J][1] = 1./(1./dx2+2./dy2)*(Press[J][2]/dx2+(Press[J+1][1]+Press[J-1][1])/dy2-rhs)
            #-------------------------------------------------------------#
            # away from the boundaries                                    #
            #-------------------------------------------------------------#
            for I in range(2, Min):
                rhs     = (Ustar[J][I]-Ustar[J][I-1])/DX + (Vstar[J][I]-Vstar[J-1][I])/DY
                Press[J][I] = 1./(2./dx2+2./dy2)*((Press[J][I+1]+Press[J][I-1])/dx2+(Press[J+1][I]+Press[J-1][I])/dy2-rhs)
            #-------------------------------------------------------------#
            # right                                                       #
            #-------------------------------------------------------------#
            rhs     = (Ur[J]-Ustar[J][Min-1])/DX + (Vstar[J][Min]-Vstar[J-1][Min])/DY
            Press[J][Min] = 1./(1./dx2+2./dy2)*(Press[J][Min-1]/dx2+(Press[J+1][Min]+Press[J-1][Min])/dy2-rhs)
        #-----------------------------------------------------------------#
        # upper left                                                      #
        #-----------------------------------------------------------------#
        rhs     = (Ustar[Nin][1]-Ul[Nin])/DX + (Vt[1]-Vstar[Nin-1][1])/DY
        Press[Nin][1] = 1./(1./dx2+1./dy2)*(Press[Nin][2]/dx2+Press[Nin-1][1]/dy2-rhs)
        #-----------------------------------------------------------------#
        # upper                                                           #
        #-----------------------------------------------------------------#
        for I in range(2, Min):
            rhs         = (Ustar[Nin][I]-Ustar[Nin][I-1])/DX +\
                            (Vt[I]-Vstar[Nin-1][I])/DY
            Press[Nin][I] = 1./(2./dx2+1./dy2)*((Press[Nin][I+1]+\
                            Press[Nin][I-1])/dx2+Press[Nin-1][I]/dy2-rhs)
        #-----------------------------------------------------------------#
        # upper right                                                     #
        #-----------------------------------------------------------------#
        rhs     = (Ur[Nin]-Ustar[Nin][Min-1])/DX + (Vt[Min]-Vstar[Nin-1][Min])/DY
        Press[Nin][Min] = 1./(1./dx2+1./dy2)*(Press[Nin][Min-1]/dx2+Press[Nin-1][Min]/dy2-rhs)
        #-----------------------------------------------------------------#
        # updating norms                                                  #
        #-----------------------------------------------------------------#
        for J in range(1,Nin+1):
            for I in range(1,Min+1):
                normnew     += Press[J][I]**2
                diffnorm    += (Press[J][I]-Pold[J][I])**2
        #-----------------------------------------------------------------#
        # new norm                                                        #
        #-----------------------------------------------------------------#
        normerror = diffnorm/normnew
        #-----------------------------------------------------------------#
        # storing old pressure                                            #
        #-----------------------------------------------------------------#
        Pold = np.copy(Press)
    #---------------------------------------------------------------------#
    # Print statement                                                     #
    #---------------------------------------------------------------------#
    print("%16.5e: %5i iterations" % (normerror, gsiter))
    String = "%16.5e: %5i iterations\n" % (normerror, gsiter)

    return Press, String
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
    N   = 32        # ny
    lx  = 1.0       # x length
    ly  = 1.0       # y length
    dx  = lx/M      # x spatial step size
    dy  = ly/N      # y spatial step size
    #---------------------------------------------------------------------#
    # Path variables                                                      #
    #---------------------------------------------------------------------#
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + "%cpython-data%c"               %(sep, sep)
    if os.path.exists(data_path) is False:
        os.mkdir(data_path)
    #---------------------------------------------------------------------#
    # Writing flags                                                       #
    #---------------------------------------------------------------------#
    write_vel_field             = True
    write_terms                 = True
    write_P                     = True
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
    #---------------------------------------------------------------------#
    # Setting  the error                                                  #
    #---------------------------------------------------------------------#
    maxerror    = 1.e-12
    diverg      = 0.
    maxdiverg   = 0.
    #---------------------------------------------------------------------#
    # Set up BC's, velocity components at boundaries are all set to 0     #
    # until now -- lid-driven cavity                                      #
    #---------------------------------------------------------------------#
    #---------------------------------------------------------------------#
    # Setting the u-velocity at the top and bottom walls                  #
    #---------------------------------------------------------------------#
    for i in range(1,M+1):
        ut[i] = 1.
        ub[i] = 0.0
    #---------------------------------------------------------------------#
    # Setting the u-velocity at the left and right walls                  #
    #---------------------------------------------------------------------#
    for j in range(1,N+1):
        ur[j] = 0.0
        ul[j] = 0.0
    #---------------------------------------------------------------------#
    # transport properties                                                #
    #---------------------------------------------------------------------#
    nu      = 0.01                  # viscosity
    dt      = 0.25*dx*dx/nu/2.      # time step
    tmax    = 30.0                   # final time
    print(dt)
    #---------------------------------------------------------------------#
    # Opening file output                                                 #
    #---------------------------------------------------------------------#
    if write_vel_field:
        u_field     = open(data_path + 'u-velocity.txt', 'w')
        v_field     = open(data_path + 'v-velocity.txt', 'w')
        uout        = ""
        vout        = ""
    if write_terms:
        Luterms     = open(data_path + 'Lu-terms.txt', 'w')
        Lvterms     = open(data_path + 'Lv-terms.txt', 'w')
        Nxterms     = open(data_path + 'Nx-terms.txt', 'w')
        Nyterms     = open(data_path + 'Ny-terms.txt', 'w')
        Luout       = ""
        Lvout       = ""
        Nxout       = ""
        Nyout       = ""
    if write_P:
        fp          = open(data_path + 'P.txt', 'w')
        pout        = ""
    if write_divergence_field:
        fdiverg     = open(data_path + 'divergence_field.txt', 'w')
        divout      = ""
    #---------------------------------------------------------------------#
    # grid for plotting                                                   #
    #---------------------------------------------------------------------#
    [Xgrid,Ygrid]   = np.meshgrid(np.linspace(dx/2,lx-dx/2,M),\
                                    np.linspace(dy/2,ly-dy/2,N))
    #---------------------------------------------------------------------#
    # Main time stepping loop                                             #
    #---------------------------------------------------------------------#
    string  = ''
    t       = 0.0
    for nstep in range(1, int(1.0000001*tmax/dt)+1):
        t += dt
        #-----------------------------------------------------------------#
        # Calculating linear and non linear terms                         #
        #-----------------------------------------------------------------#
        [Lu, Lv, Nx, Ny]    = ddt_prov(M, N, dx, dy, nu, u, v, ul, ur,\
                                        ub, ut, vl, vr, vb, vt)
        #-----------------------------------------------------------------#
        # Updating u and v star                                           #
        #-----------------------------------------------------------------#
        for j in range(1, N+1):
            for i in range(1, M):
                ustar[j][i] = u[j][i] + dt*(-Nx[j][i] + Lu[j][i])
        for j in range(1, N):
            for i in range(1, M+1):
                vstar[j][i] = v[j][i] + dt*(-Ny[j][i] + Lv[j][i])
        #-----------------------------------------------------------------#
        # Updating pressure                                               #
        #-----------------------------------------------------------------#
        [P, press_string]   = calcpress(M, N, dx, dy, maxerror, ustar, \
                                        vstar, ul, ur, vb, vt)
        string              += press_string
        #-----------------------------------------------------------------#
        # Updating velocity                                               #
        #-----------------------------------------------------------------#
        for j in range(1, N+1):
            for i in range(1, M):
                u[j][i] = ustar[j][i] - (P[j][i+1]-P[j][i])/dx
        for j in range(1, N):
            for i in range(1, M+1):
                v[j][i] = vstar[j][i] - (P[j+1][i]-P[j][i])/dy
        #-----------------------------------------------------------------#
        # Computing cell divergence                                       #
        #-----------------------------------------------------------------#
        # Bottom left                                                     #
        #-----------------------------------------------------------------#
        diverg      = (u[1][1]-ul[1])/dx + (v[1][1]-vb[1])/dy
        maxdiverg   = diverg
        if write_divergence_field:
            divout  += "%16.5e"      %(diverg)
        #-----------------------------------------------------------------#
        # Across the bottom                                               #
        #-----------------------------------------------------------------#
        for i in range(2, M):
            diverg = (u[1][i]-u[1][i-1])/dx + (v[1][i]-vb[i])/dy
            if abs(diverg)>abs(maxdiverg):
                maxdiverg = diverg
            if write_divergence_field:
                divout  += "%16.5e"      %(diverg)
        #-----------------------------------------------------------------#
        # Bottom right                                                    #
        #-----------------------------------------------------------------#
        diverg  = (ur[1]-u[1][M-1])/dx + (v[1][M]-vb[M])/dy
        if abs(diverg)>abs(maxdiverg):
            maxdiverg = diverg
        if write_divergence_field:
            divout  += "%16.5e"      %(diverg)
        if write_divergence_field:
            divout  += "\n"
        #-----------------------------------------------------------------#
        # Middle left                                                     #
        #-----------------------------------------------------------------#
        for j in range(2, N):
            diverg  = (u[j][1]-ul[j])/dx + (v[j][1]-v[j-1][1])/dy
            if abs(diverg)>abs(maxdiverg):
                maxdiverg = diverg
            if write_divergence_field:
                divout  += "%16.5e"      %(diverg)
            #-------------------------------------------------------------#
            # Across the middle                                           #
            #-------------------------------------------------------------#
            for i in range(2, M):
                diverg  = (u[j][i]-u[j][i-1])/dx + (v[j][i]-v[j-1][i])/dy
                if abs(diverg)>abs(maxdiverg):
                    maxdiverg = diverg
                if write_divergence_field:
                    divout  += "%16.5e"      %(diverg)
            #-------------------------------------------------------------#
            # Middle right                                                #
            #-------------------------------------------------------------#
            diverg  = (ur[j]-u[j][M-1])/dx + (v[j][M]-v[j-1][M])/dy
            if abs(diverg)>abs(maxdiverg):
                maxdiverg = diverg
            if write_divergence_field:
                divout  += "%16.5e"      %(diverg)
            if write_divergence_field:
                divout  += "\n"
        #-----------------------------------------------------------------#
        # Top left                                                        #
        #-----------------------------------------------------------------#
        # top left
        diverg  = (u[N][1]-ul[N])/dx + (vt[1]-v[N-1][1])/dy
        if abs(diverg)>abs(maxdiverg):
            maxdiverg = diverg
        if write_divergence_field:
            divout  += "%16.5e"                  %(diverg)
        #-----------------------------------------------------------------#
        # Across the top                                                  #
        #-----------------------------------------------------------------#
        for i in range(2, M):
            diverg  = (u[N][i]-u[N][i-1])/dx + (vt[i]-v[N-1][i])/dy
            if abs(diverg)>abs(maxdiverg):
                maxdiverg = diverg
            if write_divergence_field:
                divout  += "%16.5e"                  %(diverg)
        #-----------------------------------------------------------------#
        # Top right                                                       #
        #-----------------------------------------------------------------#
        diverg = (ur[M]-u[N][M-1])/dx + (vt[M]-v[N-1][M])/dy
        if abs(diverg)>abs(maxdiverg):
            maxdiverg = diverg
        if write_divergence_field:
            divout  += "%16.5e"                  %(diverg)
        if write_divergence_field:
            divout  += "\n"
        if write_divergence_field:
            divout  += "\n"
        #-----------------------------------------------------------------#
        # Writing output files                                            #
        #-----------------------------------------------------------------#
        #-----------------------------------------------------------------#
        # velocity                                                        #
        #-----------------------------------------------------------------#
        if write_terms:
            for j in range(1, N+1):
                for i in range(1, M):
                    Luout   += "%16.5e"     %(Lu[j][i])
                Luout   += "\n"
            Luout   += "\n"
            for j in range(1, N):
                for i in range(1, M+1):
                    Lvout   += "%16.5e"     %(Lv[j][i])
                Lvout   += "\n"
            Lvout   += "\n"
            for j in range(1, N+1):
                for i in range(1, M):
                    Nxout   += "%16.5e"     %(Nx[j][i])
                Nxout   += "\n"
            Nxout   += "\n"
            for j in range(1, N):
                for i in range(1, M+1):
                    Nyout   += "%16.5e"     %(Ny[j][i])
                Nyout   += "\n"
            Nyout   += "\n"
        #-----------------------------------------------------------------#
        # Pressure                                                        #
        #-----------------------------------------------------------------#
        if write_P:
            for j in range(1, N+1):
                for i in range(1, M+1):
                    pout    += "%16.5e"     %(P[j][i])
                pout    += "\n"
            pout    += "\n"
        #-----------------------------------------------------------------#
        # Cell centered u-velocity                                        #
        #-----------------------------------------------------------------#
        if write_vel_field or plot_vel_field:
            uinterp = np.zeros((N,M))
            vinterp = np.zeros((N,M))
            for j in range(1, N+1):
                # left
                velcent = (ul[j]+u[j][1])/2.
                uinterp[j-1][0] = velcent
                if write_vel_field:
                    uout    += "%16.5e"  %(velcent)
                for i in range(2, M):
                    velcent = (u[j][i]+u[j][i-1])/2.
                    uinterp[j-1][i-1] = velcent
                    if write_vel_field:
                        uout    += "%16.5e"     %(velcent)
                # right
                velcent = (ur[j]+u[j][M-1])/2.
                uinterp[j-1][M-1] = velcent
                if write_vel_field:
                    uout   += "%16.5e"         %(velcent)
                if write_vel_field:
                    uout    += "\n"
            if write_vel_field:
                uout    += "\n"
        #-----------------------------------------------------------------#
        # Cell centered v-velocity                                        #
        #-----------------------------------------------------------------#
            for i in range(1, M+1):
                velcent         = (v[1][i]+vb[i])/2.
                vinterp[0][i-1] = velcent
                if write_vel_field:
                    vout    += "%16.5e"     %(velcent)
            if write_vel_field:
                vout    += "\n"
            for j in range(2, N):
                for i in range(1, M+1):
                    velcent = (v[j][i]+v[j-1][i])/2.
                    vinterp[j-1][i-1] = velcent
                    if write_vel_field:
                        vout    += "%16.5e"     %(velcent)
                if write_vel_field:
                    vout    += "\n"
            for i in range(1, M+1):
                velcent = (v[N-1][i]+vt[i])/2.
                vinterp[N-1][i-1] = velcent
                if write_vel_field:
                    vout    += "%16.5e"     %(velcent)
            if write_vel_field:
                vout    += "\n"
            if write_vel_field:
                vout    += "\n"
        #-----------------------------------------------------------------#
        # Plotting                                                        #
        #-----------------------------------------------------------------#
        if plot_vel_field:
            fig, ax = plt.subplots()
            q = ax.quiver(Xgrid,Ygrid,uinterp,vinterp)
            ax.set_aspect('equal')
            plt.show()
        print("%5i/%5i | Maximum velocity divergence: %16.5e"\
                        % (nstep, int(1.0000001*tmax/dt), maxdiverg))
        string += "%5i/%5i | Maximum velocity divergence: %16.5e\n"\
                        % (nstep, int(1.0000001*tmax/dt), maxdiverg)
    #---------------------------------------------------------------------#
    # Writing output                                                      #
    #---------------------------------------------------------------------#
    f = open(data_path + 'output.out', 'w')
    f.write(string)
    f.close()
    #---------------------------------------------------------------------#
    # Closing output files                                                #
    #---------------------------------------------------------------------#
    if write_vel_field:
        u_field.write(uout)
        u_field.close()
        v_field.write(vout)
        v_field.close()
    if write_terms:
        Luterms.write(Luout)
        Luterms.close()
        Lvterms.write(Lvout)
        Lvterms.close()
        Nxterms.write(Nxout)
        Nxterms.close()
        Nyterms.write(Nyout)
        Nyterms.close()
    if write_P:
        fp.write(pout)
        fp.close()
    if write_divergence_field:
        fdiverg.write(divout)
        fdiverg.close()
