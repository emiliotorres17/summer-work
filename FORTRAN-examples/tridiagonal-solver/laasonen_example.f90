program btcs
    use precision_m
    use library
    !=====================================================================!
    ! Preamble                                                            !
    !=====================================================================!
    implicit none
    integer, parameter                      :: M = 16
    real(WP), dimension(0:M)                :: u
    real(WP), dimension(0:M)                :: b, d
    real(WP), dimension(1:M)                :: a 
    real(WP), dimension(0:M-1)              :: c
    real(WP)                                :: Uo
    real(WP)                                :: dt, dx, nu, alpha, L
    integer                                 :: i, k
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    L       = 0.04_WP
    dx      = L/dble(M)
    Uo      = 40.0_WP
    nu      = 0.000217_WP
    dt      = 0.002_WP
    alpha   = nu*dt/(dx**2.0_WP)
    !=====================================================================!
    ! Initial conditions                                                  !
    !=====================================================================!
    u(0)    = Uo
    u(1:M)  = 0.0_WP
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    do k = 1, 8
        a           = alpha
        b           = -(alpha + 2)
        c           = alpha
        d(0)        = Uo  
        d(1)        = -u(1)-alpha*u(1)
        d(2:M-2)    = -u(2:M-2)
        d(M-1)      = -u(M-1) - alpha*0.0_WP
        call tri_solve(u(1:M-1), M, a, b, c, d) 
    end do
end program btcs
