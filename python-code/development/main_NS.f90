program main_NS
    !=====================================================================!
    ! Main preamble                                                       !
    !=====================================================================!
    use precision_m                                       
    use navier_stokes_library   
    implicit none
    integer, parameter                      :: M = 64, N = 64
    real(WP)                                :: dx, dy, nu
    real(WP)                                :: dt, dt1, dt2, a, b
    real(WP), dimension(0:N+1, 0:M)         :: u, ustar, uold, Lx, Nx, a_grid 
    real(WP), dimension(0:N, 0:M+1)         :: v, vstar, vold, Ly, Ny, b_grid
    real(WP), dimension(0:N+1, 0:M+1)       :: P, Pold, grad_p, dp_const
    real(WP)                                :: ut, ul, ur, ub, vl, vr, vt, vb
    real(WP)                                :: max_error
    real(WP)                                :: tfinal, t, dt_temp
    real(WP)                                :: gs_err
    real(WP)                                :: ss_val, res_u, res_v
    real(WP)                                :: rho, mu
    integer                                 :: gs_iter
    integer                                 :: i, j, k
    integer                                 :: iter, kiter
    integer                                 :: gs_iter_max
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    nu          = 0.1
    rho         = 1.0_WP
    mu          = rho*nu
    dx          = 1.0_WP/dble(M)
    dy          = 1.0_WP/dble(N)
    dt_temp     = 0.25_WP*(dx)**(2.0_WP)/nu
    tfinal      = 11.0_WP
    k           = 0
    Pold        = 0.0_WP
    dp_const    = 1.0_WP
    !---------------------------------------------------------------------!
    ! Boundary conditions                                                 !
    !---------------------------------------------------------------------!
    ul          = 0.0_WP
    ur          = 0.0_WP
    ub          = 0.0_WP
    ut          = 0.0_WP
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
    u(0,:)      = 2.0_WP*ub - u(1,:)
    u(N+1,:)    = 2.0_WP*ut - u(N,:)
    !-----------------------------------------------------------------!
    ! Setting u-y boundary conditions                                 !
    !-----------------------------------------------------------------!
    v(:,0)      = v(:,M+1)  
    v(0,:)      = vb
    v(N,:)      = vt
    !=====================================================================!
    ! Time loop                                                           !
    !=====================================================================!
    t           = 0.0_WP
    iter        = 0
    gs_iter_max = 5e03
    max_error   = 1e-08
    uold        = u
    vold        = v
    !=====================================================================!
    ! Writing variables                                                   !
    !=====================================================================!
    open(unit=1, file='data/data-64-dev/u-temp.dat')
    open(unit=2, file='data/data-64-dev/v-temp.dat')
    open(unit=3, file='data/data-64-dev/u-star-temp.dat')
    open(unit=4, file='data/data-64-dev/v-star-temp.dat')
    open(unit=5, file='data/data-64-dev/p-temp.dat')
    open(unit=120, file='data/data-64-dev/output.txt')
    10 format(300ES25.10)
    12 format(A, ES25.5, 4X, A, I5, A, I10, /, 4x, A, I10, & 
                /, 4x, A, ES25.5, & 
                /, A, ES25.16, &
                /, 4X, A, ES25.5, &
                /, 4X, A, ES27.14, &
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
        ! Star time derivative                                            !
        !-----------------------------------------------------------------!
        call time_derv_calc(Nx, Lx, Ny, Ly, M, N, u, v, dx, dy, nu)
        !-----------------------------------------------------------------!
        ! Updating star velocities                                        !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i = 0, M
                ustar(j,i)   = u(j,i) + dt*(-Nx(j,i) + Lx(j,i)) + &
                                   dt*dp_const(j,i)
            end do                  
        end do
        do j = 1, N-1
            do i = 1, M
                vstar(j,i)  = v(j,i) + dt*(-Ny(j,i) + Ly(j,i))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Setting u-x boundary conditions                                 !
        !-----------------------------------------------------------------!
        ustar(0,:)   = 2.0_WP*ub - ustar(1,:)
        ustar(N+1,:) = 2.0_WP*ut - ustar(N,:)
        !-----------------------------------------------------------------!
        ! Setting u-y boundary conditions                                 !
        !-----------------------------------------------------------------!
        vstar(0,:)   = vb
        vstar(N,:)   = vt
        vstar(:,0)   = vstar(:,M+1)  
        !-----------------------------------------------------------------!
        !-----------------------------------------------------------------!
        ! Pressure                                                        !
        !-----------------------------------------------------------------!
        call pressure_calc(P, gs_iter, gs_err, M, N, Pold, ustar, vstar, &
                            dx, dy, dt, max_error, gs_iter_max)
        Pold = P
        !-----------------------------------------------------------------!
        ! Pressure gradient                                               !
        !-----------------------------------------------------------------!
        grad_p = 0.0_WP
        do j = 1, N
            do i = 1, M
                grad_p(j,i) = (P(j,i+1) - P(j,i))/dx
            end do
        end do
        !-----------------------------------------------------------------!
        ! Velocities                                                      !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i = 0, M
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
        !-----------------------------------------------------------------!
        ! Setting u and v boundary conditions                             !
        !-----------------------------------------------------------------!
        v(0,:)  = vb
        v(N,:)  = vt
        v(:,0)  = v(:,M+1)  
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
                'maximum pressure -->', maxval(abs(P)), &
                'convergence error -->', gs_err, &
                'steady state check -->', ss_val, &
                'mid u-velocity -->', 0.5*(u(N/2, M/2) + u((N+2)/2,M/2)), &
                'mid v-velocity -->', 0.5*(v(N/2, M/2) + v(N/2,(M+2)/2))
        !-----------------------------------------------------------------!
        ! Printout                                                        !
        !-----------------------------------------------------------------!
        print '(A, ES25.5, 4X, A, I10, A, I10)', &
                'time --> ', t, &
                'time step -->', k, '/', kiter
        if ((a+b)*dx/nu > 2.0_WP) then
            print *, 'Unstable - Cell Reynolds Number'
        end if
        print '(4x, A, I10)', 'GS iters -->', gs_iter
        print '(4x, A, ES25.16)', 'Maximum pressure -->', maxval(abs(P))
        print '(4x, A, ES25.16)', 'Maximum grad(P)-->', dt*maxval(abs(grad_p(1:N, 1:M)))
        print '(4x, A, ES25.16)', 'Maximum u -->', maxval(abs(ustar))
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
