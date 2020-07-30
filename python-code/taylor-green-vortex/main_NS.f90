program main_NS
    !=====================================================================!
    ! Main preamble                                                       !
    !=====================================================================!
    use precision_m                                       
    use navier_stokes_library   
    implicit none
    integer, parameter                      :: M = 128, N = 128, kmax = 20000
    real(WP), parameter                     :: pi = 4.0_WP*atan(1.0_WP)
    real(WP)                                :: dx, dy, nu
    real(WP)                                :: dt, dt1, dt2, a, b
    real(WP), dimension(0:N+1, 0:M)         :: u, ustar, uold, Lx, Nx, a_grid 
    real(WP), dimension(0:N, 0:M+1)         :: v, vstar, vold, Ly, Ny, b_grid
    real(WP), dimension(0:N+1, 0:M+1)       :: P, Pold, grad_p, dp_const
    real(WP)                                :: max_error
    real(WP)                                :: tfinal, t, dt_temp
    real(WP)                                :: gs_err
    real(WP)                                :: ss_val, res_u, res_v
    real(WP)                                :: rho, mu
    integer                                 :: gs_iter
    integer                                 :: i, j, k
    integer                                 :: iter, kiter
    integer                                 :: gs_iter_max
    real(WP), dimension(0:M)                :: x1
    real(WP), dimension(0:N+1)              :: y1
    real(WP), dimension(0:M+1)              :: x2
    real(WP), dimension(0:N)                :: y2
    real(WP), dimension(0:kmax)             :: u_max, v_max, tval
    !=====================================================================!
    ! Files and formats                                                   !
    !=====================================================================!
    print *, 'here'
    !---------------------------------------------------------------------!
    ! Initial conditions                                                  !
    !---------------------------------------------------------------------!
    24 format(ES25.16)
    23 format(256ES25.16)
    open(unit=101, file='data/decay-data-2/u-x.dat')
    open(unit=102, file='data/decay-data-2/u-y.dat')
    open(unit=103, file='data/decay-data-2/u-IC.dat')
    open(unit=104, file='data/decay-data-2/v-x.dat')
    open(unit=105, file='data/decay-data-2/v-y.dat')
    open(unit=106, file='data/decay-data-2/v-IC.dat')
    !---------------------------------------------------------------------!
    ! Time varying files                                                  !
    !---------------------------------------------------------------------!
    open(unit=1, file='data/decay-data-2/u-temp.dat')
    open(unit=2, file='data/decay-data-2/v-temp.dat')
    open(unit=3, file='data/decay-data-2/u-star-temp.dat')
    open(unit=4, file='data/decay-data-2/v-star-temp.dat')
    open(unit=5, file='data/decay-data-2/p-temp.dat')
    open(unit=120, file='data/decay-data-2/output.txt')
    10 format(300ES25.10)
    12 format(A, ES25.5, 4X, A, I5, A, I10, /, 4x, A, I10, & 
                /, 4x, A, ES25.5, & 
                /, A, ES25.16, &
                /, 4X, A, ES25.5, &
                /, 4X, A, ES27.14, &
                /, 4X, A, ES27.14)
    13 format(/)
    !---------------------------------------------------------------------!
    ! Maximum velocities and corresponding times                          !
    !---------------------------------------------------------------------!
    print *, 'here-->2'
    11 format(3ES35.18)
    print *, 'here-->3'
    open(unit=206, file='data/decay-data-2/maximum-velocities.dat')
    print *, 'here-->4'
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    nu          = 0.1
    rho         = 1.0_WP
    mu          = rho*nu
    dx          = 2.0_WP*pi/dble(M)
    dy          = 2.0_WP*pi/dble(N)
    dt_temp     = 0.25_WP*(dx)**(2.0_WP)/nu
    tfinal      = 20.0_WP
    k           = 0
    Pold        = 0.0_WP
    dp_const    = 1.0_WP
    !=====================================================================!
    ! Preallocating velocities and initializing velocities                !
    !=====================================================================!
    !---------------------------------------------------------------------!
    ! Maximum velocities                                                  !
    !---------------------------------------------------------------------!
    u_max       = 0.0
    v_max       = 0.0
    tval        = 0.0
    !---------------------------------------------------------------------!
    ! Velocities                                                          !
    !---------------------------------------------------------------------!
    u           = 0.0_WP
    v           = 0.0_WP
    ustar       = 0.0_WP
    vstar       = 0.0_WP
    !---------------------------------------------------------------------!
    ! Setting the x and y domain for (u-x-y stagger grid)                 !
    !---------------------------------------------------------------------!
    x1(0)   = 0.0_WP
    do i = 1, M
        x1(i)   = x1(i-1) + dx
    end do 
    y1(0)   = -0.5_WP*dx
    do j = 1, N+1
        y1(j)   = y1(j-1) +  dx
    end do
    !---------------------------------------------------------------------!
    ! Initial conditions for u(x,y,t=0)                                   !
    !---------------------------------------------------------------------!
    do j = 1,N
        do i = 0,M
            u(j,i)  = cos(x1(i))*sin(y1(j))
        end do
    end do
    !---------------------------------------------------------------------!
    ! Setting the x and y domain for (v-x-y stagger grid)                 !
    !---------------------------------------------------------------------!
    x2(0)   = -0.5_WP*dx
    do i = 1, M+1
        x2(i)   = x2(i-1) + dx
    end do 
    y2(0)   = 0.0
    do j = 1, N
        y2(j)   = y2(j-1) +  dx
    end do
    !---------------------------------------------------------------------!
    ! Initial conditions for v(x,y,t=0)                                   !
    !---------------------------------------------------------------------!
    do j = 1,N-1
        do i = 0,M+1
            v(j,i)  = -sin(x2(i))*cos(y2(j))
        end do
    end do
    !=====================================================================!
    ! Applying BCs                                                        !
    !=====================================================================!
    !---------------------------------------------------------------------!
    ! Setting u-x boundary conditions                                     !
    !---------------------------------------------------------------------!
    u(0,:)      = -u(1,:)
    u(N+1,:)    = -u(N,:)
    !---------------------------------------------------------------------!
    ! Setting u-y boundary conditions                                     !
    !---------------------------------------------------------------------!
    !v(:,M+1)    = v(:,0)  
    v(0,:)      = -sin(x2)*cos(0.0_WP)
    v(N,:)      = -sin(x2)*cos(2.0_WP*pi)
    !v(:,M+1)    = v(:,0)
    !---------------------------------------------------------------------!
    ! Initial condition for u                                             !
    !---------------------------------------------------------------------!
    do i = 0, M
        write (101, 24) x1(i)
        print '(f25.16)', x1(i)
    end do
    do j = 0, N+1
        write (102, 24) y1(j)
    end do
    do j = 0, N+1
        write (103, 23) (u(j,i), i=0, M)
    end do
    !---------------------------------------------------------------------!
    ! Initial condition for v                                             !
    !---------------------------------------------------------------------!
    do i = 0, M+1
        write (104, 24) x2(i)
    end do
    do j = 0, N
        write (105, 24) y2(j)
    end do
    do j = 0, N
        write (106, 23) (v(j,i), i=0, M+1)
    end do
    !=====================================================================!
    ! Time loop                                                           !
    !=====================================================================!
    t           = 0.0_WP
    iter        = 0
    gs_iter_max = 3e05
    max_error   = 1e-5
    uold        = u
    vold        = v
    !---------------------------------------------------------------------!
    ! Storing the maximum velocity output                                 !
    !---------------------------------------------------------------------!
    u_max(0)    = maxval(u)
    v_max(0)    = maxval(v)
    tval(0)     = t
    !=====================================================================!
    ! Writing variables                                                   !
    !=====================================================================!
    do while (t < tfinal)
        k       = k + 1
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
                ustar(j,i)   = u(j,i) + dt*(-Nx(j,i) + Lx(j,i))
            end do                  
        end do
        do j = 1, N-1
            do i = 0, M
                vstar(j,i)  = v(j,i) + dt*(-Ny(j,i) + Ly(j,i))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Setting u-x boundary conditions                                 !
        !-----------------------------------------------------------------!
        ustar(0,:)   = -ustar(1,:)
        ustar(N+1,:) = -ustar(N,:) 
        !-----------------------------------------------------------------!
        ! Setting u-y boundary conditions                                 !
        !-----------------------------------------------------------------!
        vstar(0,:)      = -sin(x2)*exp(-2.0_WP*nu*t)
        vstar(N,:)      = -sin(x2)*exp(-2.0_WP*nu*t)
        vstar(:,M+1)    = vstar(:,0)  
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
            do i = 0, M
                u(j,i)  = ustar(j,i) - dt/dx*(P(j,i+1) - P(j,i))
            end do
        end do
        do j = 1, N-1
            do i = 0, M
                v(j,i)  = vstar(j,i) - dt/dy*(P(j+1,i) - P(j,i))   
            end do
        end do
        !-----------------------------------------------------------------!
        ! Setting u-x boundary conditions                                 !
        !-----------------------------------------------------------------!
        u(0,:)   = -u(1,:) 
        u(N+1,:) = -u(N,:)
        !-----------------------------------------------------------------!
        ! Setting u and v boundary conditions                             !
        !-----------------------------------------------------------------!
        v(0,:)  = -sin(x2)*exp(-2.0_WP*nu*t)
        v(N,:)  = -sin(x2)*exp(-2.0_WP*nu*t)
        v(:,M+1)= v(:,0)  
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
        !-----------------------------------------------------------------!
        ! Updating iterations and getting the maximum velocities          !
        !-----------------------------------------------------------------!
        print   *, maxval(u)
        print   *, maxval(V)
        u_max(k)    = maxval(u)
        v_max(k)    = maxval(v)
        tval(k)     = t
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
    !---------------------------------------------------------------------!
    ! Writing maximum velocities and corresponding times                  !
    !---------------------------------------------------------------------!
    print *, k
    do i = 0,k
        write(206, 11) u_max(i), v_max(i), tval(i)
    end do

end program main_NS
