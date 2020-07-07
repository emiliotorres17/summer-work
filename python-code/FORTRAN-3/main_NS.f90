program main_NS
    !=====================================================================!
    ! Main preamble                                                       !
    !=====================================================================!
    use precision_m                                       
    use navier_stokes_library   
    implicit none
    integer, parameter                      :: M = 128, N = 128
    real(WP)                                :: dx, dy, nu
    real(WP)                                :: dt, dt1, dt2, a, b
    real(WP), dimension(0:N+1, 0:M)         :: u, ustar, uold, Lx, Nx, a_grid 
    real(WP), dimension(0:N, 0:M+1)         :: v, vstar, vold, Ly, Ny, b_grid
    real(WP), dimension(0:N+1, 0:M+1)       :: P, Pold
    real(WP)                                :: ul, ur, ut, ub, vl, vr, vt, vb
    real(WP)                                :: max_error
    real(WP)                                :: tfinal, t, dt_temp
    real(WP)                                :: gs_err
    real(WP)                                :: ss_val, res_u, res_v
    integer                                 :: gs_iter
    integer                                 :: i, j, k
    integer                                 :: iter, kiter
    integer                                 :: gs_iter_max
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    nu          = 0.005_WP
    dx          = 1.0_WP/dble(M)
    dy          = 1.0_WP/dble(N)
    dt_temp     = 0.25_WP*(dx)**(2.0_WP)/nu
    tfinal      = 50.0_WP
    k           = 0
    Pold        = 0.0_WP
    !---------------------------------------------------------------------!
    ! Boundary conditions                                                 !
    !---------------------------------------------------------------------!
    ul          = 0.0_WP
    ur          = 0.0_WP
    ut          = 1.0_WP
    ub          = 0.5_WP
    vl          = 0.0_WP
    vr          = 0.0_WP
    vt          = 0.0_WP
    vb          = 0.0_WP
    !=====================================================================!
    ! Preallocating velocities and initializing velocities                !
    !=====================================================================!
    u           = 0.0_WP
    v           = 0.0_WP
    ustar       = 0.0_WP
    vstar       = 0.0_WP
    !---------------------------------------------------------------------!
    ! Applying BCs                                                        !
    !---------------------------------------------------------------------!
    !-----------------------------------------------------------------!
    ! Setting u-x boundary conditions                                 !
    !-----------------------------------------------------------------!
    u(0,:)   = 2.0_WP*ub - u(1,:)
    u(N+1,:) = 2.0_WP*ut - u(N,:)
    u(:,0)   = ul  
    u(:,M)   = ur 
    !-----------------------------------------------------------------!
    ! Setting u-y boundary conditions                                 !
    !-----------------------------------------------------------------!
    v(:,0)   = 2.0_WP*vl - v(:,1)  
    v(:,M+1) = 2.0_WP*vr - v(:,M)  
    v(0,:)   = vb
    v(N,:)   = vt
    !call velocity_boundary_condition(u, v, M, N, u, v, ul, ur, &
    !                                        vl, vr, ut, ub, vt, vb) 
    !call velocity_boundary_condition(ustar, vstar, M, N, ustar, vstar, ul,&
    !                                        ur, vl, vr, ut, ub, vt, vb) 
    !=====================================================================!
    ! Time loop                                                           !
    !=====================================================================!
    t           = 0.0_WP
    iter        = 0
    uold        = u
    vold        = v
    gs_iter_max = 17e03
    !=====================================================================!
    ! Writing variables                                                   !
    !=====================================================================!
    open(unit=1, file='data/u-temp.dat')
    open(unit=2, file='data/v-temp.dat')
    open(unit=3, file='data/u-star-temp.dat')
    open(unit=4, file='data/v-star-temp.dat')
    open(unit=5, file='data/p-temp.dat')
    open(unit=120, file='data/output.txt')
    10 format(300ES25.10)
    12 format(A, ES25.5, 4X, A, I5, A, I10, /, 4x, A, I10, & 
                /, 4x, A, ES25.5, /, 4X, A, ES25.5, /, 4X, A, ES27.14, &
                /, 4X, A, ES27.14)
    13 format(/)
    do while (t < tfinal)
        k = k + 1
        !-----------------------------------------------------------------!
        ! time step calculation                                           !
        !-----------------------------------------------------------------!
        a_grid  = u                             ! a value for the grid
        b_grid  = v                             ! b value for the grid
        a       = maxval(abs(a_grid))           ! maximum a value
        b       = maxval(abs(b_grid))           ! maximum b value
        dt1     = 0.25_WP*(dx)**(2.0_WP)/nu     ! parabolic time step constraint
        dt2     = dx/(a+b)                      ! hyperbolic time step constraint
        dt      = min(dt1, dt2)                 ! calculating time step
        kiter   = int(tfinal/dt + 1)
        t       = t + dt
        !-----------------------------------------------------------------!
        ! GS criteria                                                     !
        !-----------------------------------------------------------------!
        if (t < 5.0_WP) then
            max_error   = (10.0_WP)**(-2.0_WP)
        elseif (t > 5.0_WP .and. t < 10.0_WP ) then
            max_error   = (10.0_WP)**(-4.0_WP)
        elseif (t > 10.0_WP .and. t < 15.0_WP ) then
            max_error   = (10.0_WP)**(-5.0_WP)
        elseif (t > 15.0_WP .and. t < 20.0_WP ) then
            max_error   = (10.0_WP)**(-6.0_WP)
        elseif (t > 20.0_WP .and. t < 25.0_WP ) then
            max_error   = (10.0_WP)**(-7.0_WP)
        else 
            max_error   = (10.0_WP)**(-8.0_WP)
        end if



        !if (t < 1.0_WP) then
        !    max_error   = (10.0_WP)**(-3.0_WP)
        !    gs_iter_max = 200
        !elseif (t > 1.0_WP .and. t < 2.0_WP) then
        !    max_error   = (10.0_WP)**(-3.0_WP)
        !    gs_iter_max = 200
        !elseif (t > 2.0_WP .and. t < 3.0_WP) then
        !    gs_iter_max = 200
        !    max_error   = (10.0_WP)**(-3.0_WP)
        !elseif ( t > 3.0_WP .and. t < 8.0_WP) then 
        !    max_error   = (10.0_WP)**(-5.0_WP)
        !    gs_iter_max = 5000
        !elseif ( t > 8.0_WP .and. t < 10.0_WP) then 
        !    gs_iter_max = 5000
        !    max_error   = (10.0_WP)**(-7.0_WP)
        !elseif ( t > 10.0_WP .and. t < 12.0_WP) then 
        !    max_error   = (10.0_WP)**(-9.0_WP)
        !    gs_iter_max = 5000
        !elseif ( t > 12.0_WP .and. t < 14.0_WP) then 
        !    max_error   = (10.0_WP)**(-10.0_WP)
        !    gs_iter_max = 5000
        !else
        !    max_error   = (10.0_WP)**(-12.0_WP)
        !    gs_iter_max = 5000
        !end if
        !-----------------------------------------------------------------!
        ! Star time derivative                                            !
        !-----------------------------------------------------------------!
        call time_derv_calc(Nx, Lx, Ny, Ly, M, N, u, v, dx, dy, nu)
        !-----------------------------------------------------------------!
        ! Updating star velocities                                        !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i = 1, M-1
                ustar(j,i)   = u(j,i) + dt*(-Nx(j,i) + Lx(j,i))
            end do
        end do
        do j = 1, N-1
            do i = 1, M
                vstar(j,i)  = v(j,i) + dt*(-Ny(j,i) + Ly(j,i))
            end do
        end do
        !call velocity_boundary_condition(ustar, vstar, M, N, ustar, vstar, &
        !                                    ul, ur, vl, vr, ut, ub, vt, vb) 
        !-----------------------------------------------------------------!
        ! Setting u-x boundary conditions                                 !
        !-----------------------------------------------------------------!
        ustar(0,:)   = 2.0_WP*ub - ustar(1,:)
        ustar(N+1,:) = 2.0_WP*ut - ustar(N,:)
        ustar(:,0)   = ul  
        ustar(:,M)   = ur 
        !-----------------------------------------------------------------!
        ! Setting u-y boundary conditions                                 !
        !-----------------------------------------------------------------!
        vstar(0,:)   = vb
        vstar(N,:)   = vt
        vstar(:,0)   = 2.0_WP*vl - vstar(:,1)  
        vstar(:,M+1) = 2.0_WP*vr - vstar(:,M)  
        !-----------------------------------------------------------------!
        !-----------------------------------------------------------------!
        ! Pressure                                                        !
        !-----------------------------------------------------------------!
        call pressure_calc(P, gs_iter, gs_err, M, N, Pold, ustar, vstar, &
                            dx, dy, dt, max_error, gs_iter_max)
        Pold = P
        !-----------------------------------------------------------------!
        ! Velocities                                                      !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i = 1, M-1
                u(j,i)  = ustar(j,i) - dt/dx*(P(j,i+1) - P(j,i))   
            end do
        end do
        do j = 1, N-1
            do i = 1, M
                v(j,i)  = vstar(j,i) - dt/dy*(P(j+1,i) - P(j,i))   
            end do
        end do
        !-----------------------------------------------------------------!
        ! Setting u-x boundary conditions                                 !
        !-----------------------------------------------------------------!
        u(0,:)   = 2.0_WP*ub - u(1,:)
        u(N+1,:) = 2.0_WP*ut - u(N,:)
        u(:,0)   = ul  
        u(:,M)   = ur 
        !-----------------------------------------------------------------!
        ! Setting u and v boundary conditions                             !
        !-----------------------------------------------------------------!
        v(:,0)   = 2.0_WP*vl - v(:,1)  
        v(:,M+1) = 2.0_WP*vr - v(:,M)  
        v(0,:)   = vb
        v(N,:)   = vt
        !-----------------------------------------------------------------!
        ! Steady state check                                              !
        !-----------------------------------------------------------------!
        res_u   = maxval(abs(u-uold)/dt)
        res_v   = maxval(abs(v-vold)/dt)
        ss_val  = max(res_u, res_v)
        uold    = u
        vold    = v
        !-----------------------------------------------------------------!
        ! Writing output                                                  !
        !-----------------------------------------------------------------!
        write(120, 12)&
                'time --> ', t, &
                'time step -->', k, '/', kiter, &
                'GS iters -->', gs_iter, &
                'convergence error -->', gs_err, &
                'steady state check -->', ss_val, &
                'mid u-velocity -->', 0.5*(u(N/2, M/2) + u((N+2)/2,M/2)), &
                'mid v-velocity -->', 0.5*(v(N/2, M/2) + v(N/2,(M+2)/2))
        !-----------------------------------------------------------------!
        ! Printout                                                        !
        !-----------------------------------------------------------------!
        print '(A, ES25.5, 4X, A, I5, A, I10)', &
                'time --> ', t, &
                'time step -->', k, '/', kiter
        if ((a+b)*dx/nu > 2.0_WP) then
            print *, 'Unstable - Cell Reynolds Number'
        end if
        print '(4x, A, I10)', 'GS iters -->', gs_iter
        print '(4x, A, ES25.5)', 'convergence error -->', gs_err
        print '(4X, A, ES25.5)', 'steady state check -->', ss_val
        print '(4X, A, ES27.14)', 'mid u-velocity -->', 0.5*(u(N/2, M/2) + u((N+2)/2,M/2))
        print '(4X, A, ES27.14)', 'mid v-velocity -->', 0.5*(v(N/2, M/2) + v(N/2,(M+2)/2))
    end do
    !---------------------------------------------------------------------!
    ! Writing velocities                                                  !
    !---------------------------------------------------------------------!
    do j = 0, N+1
        write(1, 10) (u(j,i), i=0, M)
    end do
    write(1, 13)
    do j = 0, N
        write(2, 10) (v(j,i), i=0, M+1)
    end do
    write(2, 13)
    !---------------------------------------------------------------------!
    ! Writing u and v star velocities                                     !
    !---------------------------------------------------------------------!
    do j = 0, N+1
        write(3, 10) (ustar(j,i), i=0, M)
    end do
    do j = 0, N
        write(4, 10) (vstar(j,i), i=0, M+1)
    end do
    !---------------------------------------------------------------------!
    ! Writing pressure                                                    !
    !---------------------------------------------------------------------!
    do j = 0, N+1
        write(5, 10) (P(j,i), i=0, M+1)
    end do
    write(5, 13)
end program main_NS
