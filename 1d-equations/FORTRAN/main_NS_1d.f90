program navier_stokes_1D
    !=====================================================================!
    ! Main preamble                                                       !
    !=====================================================================!
    use precision_m                                       
    use navier_stokes_library_1d   
    implicit none 
    integer, parameter                              :: M = 4096 
    real(WP)                                        :: dy, dt, nu, dif
    real(WP)                                        :: sigma
    real(WP), dimension(0:M+1)                      :: y
    real(WP), dimension(0:M)                        :: v, vstar
    real(WP), dimension(0:M+1)                      :: u , unew
    real(WP), dimension(0:M+1)                      :: spec1, spec1_new
    real(WP)                                        :: ut, ub, vt, vb, v1, v2
    real(WP), dimension(0:M+1)                      :: p, pold
    real(WP)                                        :: dp_dx
    real(WP)                                        :: t, tfinal
    real(WP)                                        :: rdy, rdy2
    real(WP)                                        :: gs_error, gs_max
    real(WP), dimension(0:M+1)                      :: phinew
    real(WP)                                        :: phi1, phi2
    real(WP)                                        :: dt1, dt2
    integer                                         :: counter, iter_max, n_gs
    integer                                         :: j, i, q
    integer                                         :: print_count = 0
    real(WP)                                        :: phi_bot = 0.0_WP
    real(WP)                                        :: phi_top = 1.0_WP
    !---------------------------------------------------------------------!
    ! File and formats                                                    !
    !---------------------------------------------------------------------!
    10 format(f25.16)
    11 format(2f25.16)
    15 format(I8)
    14 format(/)
    open(unit=20, file='data/data-scalar-test/u-data.dat')
    open(unit=21, file='data/data-scalar-test/v-data.dat')
    open(unit=22, file='data/data-scalar-test/phi-data.dat')
    open(unit=23, file='data/data-scalar-test/pressure-data.dat')
    open(unit=24, file='data/data-scalar-test/constants.dat')
    !---------------------------------------------------------------------!
    ! Domain variables                                                    !
    !---------------------------------------------------------------------!
    dy      = 1.0_WP/dble(M)
    rdy     = 1.0_WP/dy
    rdy2    = (1.0_WP/dy)**2.0_WP
    nu      = 1.0_WP
    dif     = 0.00001_WP
    dt1     = 0.9_WP*(dy)**2.0_WP/(2.0_WP*dif)
    dt2     = (10.0_WP)**(-7.0_WP) 
    dt      = min(dt1, dt2)
    print *, dy
    print *, dt
    dp_dx   = 0.0_WP
    t       = 0.0_WP
    tfinal  = 1.0001_WP 
    !---------------------------------------------------------------------!
    ! Wall velocities                                                     !
    !---------------------------------------------------------------------!
    ut      = 0.0_WP
    ub      = 0.0_WP
    vt      = 0.1_WP
    vb      = 0.1_WP
    !---------------------------------------------------------------------!
    ! Preallocating variables                                             !
    !---------------------------------------------------------------------!
    u       = 0.0_WP
    unew    = 0.0_WP
    v       = 0.0_WP
    vstar   = vt
    p       = 0.0_WP
    pold    = 0.0_WP
    spec1   = 0.0_WP
    phinew  = 0.0_WP
    !---------------------------------------------------------------------!
    ! Setting boundary conditions                                         !
    !---------------------------------------------------------------------!
    u(0)        = 2.0_WP*ub - u(1)
    u(M+1)      = 2.0_WP*ut - u(M)
    v(0)        = vb
    v(M)        = vt
    v(1:M-1)    = vt 
    !---------------------------------------------------------------------!
    ! Convergence criteria                                                !
    !---------------------------------------------------------------------!
    counter     = 0
    gs_max      = 1e-08
    iter_max    = 100
    !!---------------------------------------------------------------------!
    !! Initial condition for phi                                           !
    !!---------------------------------------------------------------------!
    y(0)        = -0.5_WP*dy
    do i = 1, M+1
        y(i)    = y(i-1)+dy 
    end do
    !do i = 1, M
    !    spec1(i)  = -8.0_WP*(y(i)-0.5_WP)**2.0_WP + 2.0_WP
    !end do
    !spec1(0)    = 2.0_WP*phi_bot-spec1(1)
    !spec1(M+1)  = 2.0_WP*phi_top-spec1(M)
    !---------------------------------------------------------------------!
    ! Initial condition for phi                                           !
    !---------------------------------------------------------------------!
    do i = 0, M+1
        if (y(i) < 0.5) then
            spec1(i) = 0.0_WP
        else if (i >= 0.5) then
            spec1(i) = 1.0_WP
        end if
    end do
    spec1(0)    = 2.0_WP*phi_bot-spec1(1)
    spec1(M+1)  = 2.0_WP*phi_top-spec1(M)
    !---------------------------------------------------------------------!
    ! Writing initial condition                                           !
    !---------------------------------------------------------------------!
    do i = 0, M+1
        write(22, 11)    y(i), spec1(i)
    end do
    write(22,14)
    !---------------------------------------------------------------------!
    ! Time loop                                                           !
    !---------------------------------------------------------------------!
    do while (t < tfinal)
        q = q+1
        !-----------------------------------------------------------------!
        ! Updating time values                                            !
        !-----------------------------------------------------------------!
        t       = t + dt
        counter = counter + 1
        !-----------------------------------------------------------------!
        ! Calculating u                                                   !
        !-----------------------------------------------------------------!
        call u_time_derv(unew, M, u, v, dp_dx, dt, dy, nu)
        u       = unew
        u(0)    = 2.0_WP*ub - u(1)
        u(M+1)  = 2.0_WP*ut - u(M)
        !-----------------------------------------------------------------!
        ! Calculating scalar transport                                    !
        !-----------------------------------------------------------------!
        !call scalar_transport(spec1_new, M, spec1, v, dt, dy, dif)
        call rk4_density(spec1_new, M, spec1, v, dt, dy, dif)
        !do i = 1, M
        !    phi1        = 0.5_WP*(spec1(i) + spec1(i+1))    
        !    phi2        = 0.5_WP*(spec1(i) + spec1(i-1))    
        !    phinew(i)   = spec1(i) - dt*rdy*(phi1*v(i) - phi2*v(i-1)) + &
        !                    dif*dt*rdy2*(spec1(i+1) - 2.0_WP*spec1(i) +  spec1(i-1))
        !end do
        !spec1       = phinew
        spec1       = spec1_new
        spec1(0)    = 2.0_WP*phi_bot-spec1(1)
        spec1(M+1)  = 2.0_WP*phi_top-spec1(M)
        !-----------------------------------------------------------------!
        ! Calculating vstar                                               !
        !-----------------------------------------------------------------!
        do j = 1, M-1
            v1          = 0.5_WP*(v(j) + v(j+1))
            v2          = 0.5_WP*(v(j) + v(j-1))
            vstar(j)    = v(j) - dt*rdy*(v1**2.0_WP - v2**2.0_WP) &
                                + dt*rdy2*nu*(v(j+1)-2.0_WP*v(j)+v(j-1))
        end do
        !-----------------------------------------------------------------!
        ! vstar boundary conditions                                       !
        !-----------------------------------------------------------------!
        vstar(0)    = vb
        vstar(M)    = vt
        !-----------------------------------------------------------------!
        ! Calculating the pressure                                        !
        !-----------------------------------------------------------------!
        call pressure_calc(p, n_gs, gs_error, M, pold, vstar, &
                                dy, dt, gs_max, iter_max)
        pold    = p
        !-----------------------------------------------------------------!
        ! Calculating v^{n+1}                                             !
        !-----------------------------------------------------------------!
        do j = 1, M-1
            v(j) = vstar(j) - dt*rdy*(P(j+1)-P(j))
        end do
        v(0)    = vb
        v(M)    = vt
        !-----------------------------------------------------------------!
        ! Print statement                                                 !
        !-----------------------------------------------------------------!
        if (counter > 10) then
            print '(A, I8)', 'q -->', q
            print '(A, ES25.16)', 'time -->', t
            print '(4x, A, I10)', 'GS iters -->', n_gs
            print '(4x, A, ES25.16)', 'convergence error -->', gs_error
            print '(4x, A, ES25.16)', 'Maximum u -->', maxval(abs(u))
            print '(4X, A, ES27.14)', 'mid u-velocity -->', 0.5*(u(M/2) + u(M/2+1))
            print '(4X, A, ES27.14)', 'mid v-velocity -->', 0.5*(v(M/2) + v(M/2+1))
            print '(4X, A, ES27.14)', 'phi average -->', sum(spec1)/M
            print *, '/'
            counter = 0
        end if
        !-----------------------------------------------------------------!
        ! Time dependent printout                                         !
        !-----------------------------------------------------------------!
        if  (t > 0.01_WP .and. print_count == 0) then
            !-------------------------------------------------------------!
            ! Print phi field at t=1                                      !
            !-------------------------------------------------------------!
            do i = 0, M+1
                write(22, 11)    y(i), spec1(i)
            end do
            write(22,14)
            print_count = print_count + 1
            stop
        else if  (t > 2.0_WP .and. print_count == 1) then
            !-------------------------------------------------------------!
            ! Print phi field at t=3                                      !
            !-------------------------------------------------------------!
            do i = 0, M+1
                write(22, 11)    y(i), spec1(i)
            end do
            write(22,14)
            print_count = print_count + 1
        else if  (t > 3.0_WP .and. print_count == 2) then
            !-------------------------------------------------------------!
            ! Print phi field at t=3                                      !
            !-------------------------------------------------------------!
            do i = 0, M+1
                write(22, 11)    y(i), spec1(i)
            end do
            write(22,14)
            print_count = print_count + 1
        end if 
    end do
    !---------------------------------------------------------------------!
    ! Print u velocity                                                    !
    !---------------------------------------------------------------------!
    do i = 0, M+1
        write(20, 10)    u(i)
    end do
    !---------------------------------------------------------------------!
    ! Print v velocity                                                    !
    !---------------------------------------------------------------------!
    do i = 0, M
        write(21, 10)    v(i)
    end do
    !---------------------------------------------------------------------!
    ! Print pressure field                                                !
    !---------------------------------------------------------------------!
    do i = 0, M+1
        write(23, 11)    y(i), P(i)
    end do
    print *, maxval(spec1)
    !---------------------------------------------------------------------!
    ! Printing constants                                                  !
    !---------------------------------------------------------------------!
    write(24, 10)       t
    write(24, 10)       nu
    write(24, 10)       dif
    write(24, 15)       M
end program
