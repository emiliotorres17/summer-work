program navier_stokes_1D
    !=====================================================================!
    ! Main preamble                                                       !
    !=====================================================================!
    use precision_m                                       
    use navier_stokes_library_1d   
    implicit none 
    integer, parameter                              :: M = 128
    real(WP)                                        :: dy, dt, nu, dif
    real(WP)                                        :: sigma
    real(WP), dimension(0:M+1)                      :: y
    real(WP), dimension(0:M)                        :: v, vstar
    real(WP), dimension(0:M+1)                      :: u , unew
    real(WP), dimension(M+1)                        :: spec1, spec1_new
    real(WP)                                        :: ut, ub, vt, vb, v1, v2
    real(WP), dimension(0:M+1)                      :: p, pold
    real(WP)                                        :: dp_dx
    real(WP)                                        :: t, tfinal
    real(WP)                                        :: rdy, rdy2
    real(WP)                                        :: gs_error, gs_max
    integer                                         :: counter, iter_max, n_gs
    integer                                         :: j, i, q
    !---------------------------------------------------------------------!
    ! File and formats                                                    !
    !---------------------------------------------------------------------!
    10 format(f25.16)
    11 format(2f25.16)
    open(unit=20, file='data/data-scalar-transport/u-data.dat')
    open(unit=21, file='data/data-scalar-transport/v-data.dat')
    open(unit=22, file='data/data-scalar-transport/phi-data.dat')
    open(unit=23, file='data/data-scalar-transport/pressure-data.dat')
    !---------------------------------------------------------------------!
    ! Domain variables                                                    !
    !---------------------------------------------------------------------!
    dy      = 1.0_WP/dble(M) 
    rdy     = 1.0_WP/dy
    rdy2    = (1.0_WP/dy)**2.0_WP
    nu      = 0.1_WP
    dif     = 0.01_WP
    !dt      = 0.25_WP*dy**2.0_WP/nu
    dt      = 0.0001_WP
    dp_dx   = 0.0_WP
    t       = 0.0_WP
    tfinal  = 10.0_WP 
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
    vstar   = 0.0_WP
    p       = 0.0_WP
    pold    = 0.0_WP
    spec1   = 0.0_WP
    !---------------------------------------------------------------------!
    ! Setting boundary conditions                                         !
    !---------------------------------------------------------------------!
    u(0)    = 2.0_WP*ub - u(1)
    u(M+1)  = 2.0_WP*ut - u(M)
    v(0)    = vb
    v(M)    = vt
    !---------------------------------------------------------------------!
    ! Convergence criteria                                                !
    !---------------------------------------------------------------------!
    counter     = 0
    gs_max      = 1e-08
    iter_max    = 100
    !---------------------------------------------------------------------!
    ! Initial condition for phi                                           !
    !---------------------------------------------------------------------!
    sigma       = 1.0_WP/6.0_WP
    y(0)        = -0.5_WP*dy
    do i = 1, M+1
        y(i)    = y(i-1)+dy 
    end do
    do i = 1, M
        spec1(i)  = -8.0_WP*(y(i)-0.5_WP)**2.0_WP + 2.0_WP
    end do
    spec1(0)    = -spec1(1)
    spec1(M+1)  = -spec1(M)
    !---------------------------------------------------------------------!
    ! Time loop                                                           !
    !---------------------------------------------------------------------!
    do while (t < tfinal)
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
        call scalar_transport(spec1_new, M, spec1, v, dt, dy, dif)
        spec1       = spec1_new
        spec1(0)    = -spec1(1)
        spec1(M+1)  = -spec1(M)
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
        !call pressure_calc(p, n_gs, gs_error, M, pold, vstar, &
        !                        dy, dt, gs_max, iter_max)
        !pold    = p
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
        if (counter > 100) then
            print '(A, I8)', 'q -->', q
            print '(A, ES25.16)', 'time -->', t
            print '(4x, A, I10)', 'GS iters -->', n_gs
            print '(4x, A, ES25.16)', 'convergence error -->', gs_error
            print '(4x, A, ES25.16)', 'Maximum u -->', maxval(abs(u))
            print '(4X, A, ES27.14)', 'mid u-velocity -->', 0.5*(u(M/2) + u(M/2+1))
            print '(4X, A, ES27.14)', 'mid v-velocity -->', 0.5*(v(M/2) + v(M/2+1))
            print *, '/'
            counter = 0
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
    ! Print phi field                                                     !
    !---------------------------------------------------------------------!
    do i = 0, M+1
        write(22, 11)    y(i), spec1(i)
    end do
    !---------------------------------------------------------------------!
    ! Print pressure field                                                !
    !---------------------------------------------------------------------!
    do i = 0, M+1
        write(23, 11)    y(i), P(i)
    end do
    print *, maxval(spec1)
end program
