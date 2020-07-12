program btcs
    use precision_m
    use library
    !=====================================================================!
    ! Preamble                                                            !
    !=====================================================================!
    implicit none
    integer, parameter                      :: M = 64
    real(WP), dimension(0:M)                :: u, u18, u36, u54, u72, u90
    real(WP), dimension(1:M-1)              :: b, d, sol
    real(WP), dimension(2:M-1)              :: a 
    real(WP), dimension(1:M-2)              :: c
    real(WP)                                :: Uo
    real(WP)                                :: dt, dx, nu, alpha, L
    real(WP)                                :: t, tfinal
    integer                                 :: i, k
    !=====================================================================!
    ! File formats                                                        !
    !=====================================================================!
    10 format(ES25.10)
    open(unit=1, file='data/u-IC.dat')
    open(unit=2, file='data/u-final.dat')
    open(unit=3, file='data/u-18.dat')
    open(unit=4, file='data/u-36.dat')
    open(unit=5, file='data/u-54.dat')
    open(unit=6, file='data/u-72.dat')
    open(unit=7, file='data/u-90.dat')
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    L       = 0.04_WP
    tfinal  = 1.0_WP
    t       = 0.0_WP
    dx      = L/dble(M)
    Uo      = 40.0_WP
    nu      = 0.000217_WP
    dt      = 0.25_WP*dx
    alpha   = nu*dt/(dx**2.0_WP)
    !=====================================================================!
    ! Preallocating storing variables                                     !
    !=====================================================================!
    u18     = 0.0_WP
    u36     = 0.0_WP
    u54     = 0.0_WP
    u72     = 0.0_WP
    u90     = 0.0_WP
    !=====================================================================!
    ! Initial conditions                                                  !
    !=====================================================================!
    u(0)    = Uo
    u(1:M)  = 0.0_WP
    sol     = 0.0_WP
    d       = 0.0_WP
    b       = 0.0_WP
    c       = 0.0_WP
    a       = 0.0_WP
    !=====================================================================!
    ! Write initial condition                                             !
    !=====================================================================!
    do i = 0, M
        write(1, 10) u(i) 
    end do
    close(unit=1)
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    k   = 1
    do while (t < tfinal)
        t           = t + dt
        a           = alpha
        b           = -(2.0_WP*alpha + 1)
        c           = alpha
        d(1)        = -u(1)-alpha*u(0)
        d(2:M-2)    = -u(2:M-2)
        d(M-1)      = -u(M-1) - alpha*u(M)
        call tri_solve(sol, M, a, b, c, d) 
        u(1:M-1)    = sol
        u(0)        = Uo
        u(M)        = 0.0_WP
        !-----------------------------------------------------------------!
        ! Writing different solutions                                     !
        !-----------------------------------------------------------------!
        if ( t > 0.18 .and. k == 1) then
            u18 = u
            k   = k + 1
        end if
        if ( t > 0.36 .and. k == 2) then
            u36 = u
            k   = k + 1
        end if
        if ( t > 0.54 .and. k == 3) then
            u54 = u
            k   = k + 1
        end if
        if ( t > 0.72 .and. k == 4) then
            u72 = u
            k   = k + 1
        end if
        if ( t > 0.90 .and. k == 5) then
            u90 = u
            k   = k + 1
        end if
    end do
    !=====================================================================!
    ! Writing solutions                                                   !
    !=====================================================================!
    do i = 0, M
        write(2, 10) u(i) 
        write(3, 10) u18(i) 
        write(4, 10) u36(i) 
        write(5, 10) u54(i) 
        write(6, 10) u72(i) 
        write(7, 10) u90(i) 
    end do
    close(unit=2)
    close(unit=3)
    close(unit=4)

end program btcs
