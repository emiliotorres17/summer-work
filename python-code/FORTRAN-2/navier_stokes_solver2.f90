program navier_stokes
    use precision_m    
    use ns_lib 
    implicit none
    integer,    parameter                           :: M = 32, N = 32
    real(WP),   dimension(0:N)                      :: x, y
    real(WP)                                        :: dx, dy, dt, nu
    real(WP),   dimension(0:N)                      :: ul, ur, vl, vr
    real(WP),   dimension(0:M)                      :: ub, ut, vb, vt
    real(WP),   dimension(0:N, 0:M-1)               :: u, ustar, uold
    real(WP),   dimension(0:N-1,0:M)                :: v, vstar, vold
    real(WP),   dimension(0:N, 0:M-1)               :: Lu
    real(WP),   dimension(0:N, 0:M-1)               :: Nx
    real(WP),   dimension(0:N-1, 0:M)               :: Lv, Ny
    real(WP),   dimension(0:N, 0:M)                 :: P
    real(WP)                                        :: maxerror
    integer                                         :: i, j, nstep, tstep_F
    real(WP)                                        :: tmax, t
    integer                                         :: counter
    real(WP)                                        :: max_div
    real(WP)                                        :: div
    real(WP)                                        :: res_u, res_v
    real(WP)                                        :: ss_check, ss_max
    real(WP)                                        :: dist
    !---------------------------------------------------------------------!
    ! Preallocating variables                                             !
    !---------------------------------------------------------------------!
    counter     = 1
    u           = 0.0_WP
    uold        = 0.0_WP
    v           = 0.0_WP
    vold        = 0.0_WP
    ustar       = 0.0_WP
    vstar       = 0.0_WP
    !---------------------------------------------------------------------!
    ! Boundary conditions                                                 !
    !---------------------------------------------------------------------!
    ut          = 1.0_WP
    ub          = 0.0_WP
    ul          = 0.0_WP
    ur          = 0.0_WP
    vt          = 0.0_WP
    vb          = 0.0_WP
    vl          = 0.0_WP
    vr          = 0.0_WP
    !---------------------------------------------------------------------!
    ! Transport properties                                                !
    !---------------------------------------------------------------------!
    nu          = 0.01_WP
    dx          = 1.0_WP/real(M)
    dy          = 1.0_WP/real(N)
    dt          = 0.25_WP*dx*dy/nu/2.0_WP
    !---------------------------------------------------------------------!
    ! Simulation running criteria                                         !
    !---------------------------------------------------------------------!
    tmax        = 2.0_WP
    maxerror    = (10.0_WP)**(-12.0_WP)
    tstep_F     = int(1.0000001_WP*tmax/dt)
    ss_check    = 1.00
    ss_max      = 1.0e-7
    print '(A,/,4X,A,f16.10,/,4X,A,f16.10,/,4X,A,f16.10, &
            & /,4X,A,f16.10,/,4X,A,ES16.5)', &
            'Simulation running criteria', &    
            'time step -->', dt, &
            'x-spatial step size -->', dx, &
            'y-spatial step size -->', dy, &
            'steady state criteria -->', ss_check, &
            'maximum error criteria -->', maxerror
    !---------------------------------------------------------------------!
    ! Print variables                                                     !
    !---------------------------------------------------------------------!
    open(unit=1, file='FORTRAN-data/FORTRAN-32-data/output.out')
    10 format(I8, A, I8, 4X, A, f10.6, 4X, A, ES16.5, /, 4X, A, &
                ES16.5,/,4X, A, ES16.5)
    !---------------------------------------------------------------------!
    ! Time stepping loop                                                  !
    !---------------------------------------------------------------------!
    nstep   = 0
    do while (nstep < tstep_F .and. ss_check > ss_max) 
        !-----------------------------------------------------------------!
        ! Updating time step and simulation time                          !
        !-----------------------------------------------------------------!
        nstep   = nstep + 1
        t       = t + dt   
        !-----------------------------------------------------------------!
        ! Calculating linear and non-linear derivatives                   !
        !-----------------------------------------------------------------!
        call time_derv(Lu, Lv, Nx, Ny, M, N, dx, dy, nu, u, v, ul, ur, &
                            ub, ut, vl, vr, vb, vt)
        !-----------------------------------------------------------------!
        ! Updating ustar and vstar                                        !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i  = 1, M-1
                ustar(j,i) = u(j,i) + dt*(-Nx(j,i) + Lu(j,i))
            end do
        end do
        do j = 1, N-1
            do i = 1, M
                vstar(j,i) = v(j,i) + dt*(-Ny(j,i) + Lv(j,i))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Updating pressure                                               !
        !-----------------------------------------------------------------!
        call calcpress(P, M, N, dx, dy, maxerror, ustar, vstar, ul, &
                                ur, vb, vt)
        !-----------------------------------------------------------------!
        ! Updating u and v                                                !
        !-----------------------------------------------------------------!
        do j  = 1, N
            do i = 1, M-1
                u(j,i) = ustar(j,i) - (P(j,i+1)-P(j,i))/dx
            end do
        end do 
        do j  = 1, N-1
            do i  = 1, M
                v(j,i) = vstar(j,i) - (P(j+1,i)-P(j,i))/dy
            end do 
        end do
        !!-----------------------------------------------------------------!
        !! Update boundary conditions                                      !
        !!-----------------------------------------------------------------!
        !ul          = u(:,1)
        !ur          = ul
        !vl(0:N-1)   = v(:,1)
        !vr          = vl
        !-----------------------------------------------------------------!
        ! Calculating velocity divergence                                 !
        !-----------------------------------------------------------------!
        call divergence(max_div, div, M, N, u, v, ul, ur, vt, vb, dx, dy)
        !-----------------------------------------------------------------!
        ! Steady state check                                              !
        !-----------------------------------------------------------------!
        res_u       = maxval(abs((u-uold)/dt))
        res_v       = maxval(abs((u-uold)/dt))
        ss_check    = max(res_u, res_v)
        uold        = u
        vold        = v
        !-----------------------------------------------------------------!
        ! Print statement                                                 !
        !-----------------------------------------------------------------!
        print '(I8, A, I8, 4X, A, f10.6, 4X, A, ES10.3)', &
                nstep, '/', tstep_F, &
                'time -->', t, &
                'max velocity -->', max(maxval(u),maxval(v))
        !-----------------------------------------------------------------!
        ! Divergence print out                                            !
        !-----------------------------------------------------------------!
        print '(4X, A, ES16.5)', &
            '| Maximum velocity diveregence -->', max_div
        !-----------------------------------------------------------------!
        ! Steady state print out                                          !
        !-----------------------------------------------------------------!
        print '(4X, A, ES16.5)', &
                '| Steady state check -->', ss_check
        !-----------------------------------------------------------------!
        ! Writing output                                                  !
        !-----------------------------------------------------------------!
        write(1,  10)&
                nstep, '/', tstep_F, &
                'time -->', t, &
                'max velocity -->', max(maxval(u),maxval(v)), &
                '| Maximum velocity diveregence -->', max_div, &
                '| Steady state check -->', ss_check
        !-----------------------------------------------------------------!
        ! Writing u-velocity                                              !
        !-----------------------------------------------------------------!
        !-----------------------------------------------------------------!
        ! Cell centered u and v velocities                                !
        !-----------------------------------------------------------------!
        call u_write(u,M,N,ul,ur)
        call v_write(v,M,N,vt,vb)
        !-----------------------------------------------------------------!
        ! Linear derivative terms                                         !
        !-----------------------------------------------------------------!
        call Lv_write(Lv,M,N) 
        call Lu_write(Lu,M,N) 
        !-----------------------------------------------------------------!
        ! Non-linear derivatives                                          !
        !-----------------------------------------------------------------!
        call Ny_write(Ny,M,N)
        call Nx_write(Nx,M,N)
        !-----------------------------------------------------------------!
        ! Pressure                                                        !
        !-----------------------------------------------------------------!
        call P_write(P,M,N) 
    end do
    close(unit=1)
    close(unit=129)
end program navier_stokes
