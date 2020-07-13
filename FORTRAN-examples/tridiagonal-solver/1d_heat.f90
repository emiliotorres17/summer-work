program heat_equation
    use precision_m
    use library
    !=====================================================================!
    ! Preamble                                                            !
    !=====================================================================!
    implicit none
    real(WP), parameter                     :: pi = 4.0_WP*atan(1.0_WP)
    integer, parameter                      :: M = 16, Mex = 128
    real(WP), dimension(0:M)                :: u, usol, ucn, x
    real(WP), dimension(0:Mex)              :: uex, xex
    real(WP), dimension(1:M-1)              :: b, d, sol, sol_cn
    real(WP), dimension(2:M-1)              :: a 
    real(WP), dimension(1:M-2)              :: c
    real(WP)                                :: dt, dx, L, alpha, dx_e
    real(WP)                                :: t, tfinal
    integer                                 :: i, n
    real(WP)                                :: k
    !=====================================================================!
    ! Defining files                                                      !
    !=====================================================================!
    10 format(2ES25.12)
    open(unit=1, file='data/u-IC.dat') 
    open(unit=2, file='data/u-exact.dat') 
    open(unit=3, file='data/u-approx.dat') 
    open(unit=4, file='data/u-cn.dat') 
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    L       = 1.0_WP
    dx      = L/dble(M)
    dx_e    = L/dble(Mex)
    dt      = 0.25 
    k       = 0.01_WP
    alpha   = dt*k/dx**2.0_WP 
    t       = 0.0_WP
    tfinal  = 1.0_WP
    !=====================================================================!
    ! Preallocating variables                                             !
    !=====================================================================!
    u       = 0.0_WP
    usol    = 0.0_WP
    ucn     = 0.0_WP
    x(0)    = 0.0_WP
    xex(0)  = 0.0_WP
    do i = 1, M
        x(i) = x(i-1) + dx
    end do
    !=====================================================================!
    ! Initial condition                                                   !
    !=====================================================================!
    do i = 0, M 
        u(i)    = 4.0_WP*sin(3.0_WP*pi*x(i)/L)
        ucn(i)  = 4.0_WP*sin(3.0_WP*pi*x(i)/L)
        write(1, 10) x(i), u(i)    
    end do
    close(unit=1)
    !=====================================================================!
    ! Exact solution at t = tfinal                                        !
    !=====================================================================!
    do i = 1, Mex 
        xex(i)  = xex(i-1) + dx_e
    end do
    do i = 0, Mex 
        uex(i)  = 4.0_WP*sin(3.0_WP*pi*xex(i)/L)*&
                    exp(-k*(3.0_WP*pi/L)**2.0_WP*tfinal)
        write(2, 10) xex(i), uex(i)    
    end do
    close(unit=2)
    !=====================================================================!
    ! BTCS method                                                         !
    !=====================================================================!
    do while (t < 1.0)
        t = t + dt  
        print '(f25.12)',  t 
        a           = alpha
        b           = -(2.0_WP*alpha + 1)
        c           = alpha
        d(1)        = -u(1)
        d(2:M-2)    = -u(2:M-2)
        d(M-1)      = -u(M-1)
        call tri_solve(sol, M, a, b, c, d) 
        u(1:M-1)    = sol
        u(0)        = 0.0_WP
        u(M)        = 0.0_WP
    end do
    !=====================================================================!
    ! Crank-Nickelson method                                              !
    !=====================================================================!
    t   = 0
    do while (t < 1.0)
        t = t + dt  
        print '(f25.12)',  t 
        a           = 0.5_WP*alpha
        b           = -(alpha + 1)
        c           = 0.5_WP*alpha
        do i = 1, M-1
            d(i) = -0.5_WP*alpha*ucn(i+1) + (alpha - 1)*ucn(i) &
                    - 0.5_WP*alpha*ucn(i-1) 
        end do
        call tri_solve(sol_cn, M, a, b, c, d) 
        ucn(1:M-1)    = sol_cn
        ucn(0)        = 0.0_WP
        ucn(M)        = 0.0_WP
    end do
    !=====================================================================!
    ! Writing BTCS solution                                               !
    !=====================================================================!
    do i = 0, M
        write(3, 10) x(i), u(i)
    end do
    close(unit=3)
    !=====================================================================!
    ! Writing CN solution                                                 !
    !=====================================================================!
    do i = 0, M
        write(4, 10) x(i), ucn(i)
    end do
    close(unit=4)
end program heat_equation
