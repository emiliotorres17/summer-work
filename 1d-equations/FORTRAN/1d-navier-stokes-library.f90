module navier_stokes_library_1d
    !=====================================================================!
    ! Dependencies                                                        !
    !=====================================================================!
    use precision_m
    implicit none
    public :: time_derv_calc, pressure_calc, scalar_transport, u_time_derv
    contains
    !=====================================================================!
    ! Time derivate                                                       !
    !=====================================================================!
    subroutine time_derv_calc(Ny, Ly, M, v, dy, nu)
        !-----------------------------------------------------------------!
        ! Time derivative preamble                                        !
        !-----------------------------------------------------------------!
        integer,  intent(in)                                :: M
        real(WP), intent(in)                                :: nu
        real(WP), dimension(0:M), intent(in)                :: v 
        real(WP), intent(in)                                :: dy
        real(WP)                                            :: rdy, rdy2
        real(WP), dimension(0:M), intent(out)               :: Ny, Ly 
        integer                                             :: j
        real(WP)                                            :: v1, v2
        !-----------------------------------------------------------------!
        ! Domain variables                                                !
        !-----------------------------------------------------------------!
        rdy     = 1.0_WP/dy
        rdy2    = 1.0_WP/(dy)**2.0_WP
        Ly      = 0.0_WP
        Ny      = 0.0_WP
        !-----------------------------------------------------------------!
        ! Looping over domain                                             !
        !-----------------------------------------------------------------!
        do j = 1, M-1
            !-------------------------------------------------------------!
            ! Linear derivative                                           !
            !-------------------------------------------------------------!
            Ly(j)   = nu*rdy2*(v(j+1)-2.0_WP*v(j)+v(j-1))
            !-------------------------------------------------------------!
            ! Non-linear time derivative                                  !
            !-------------------------------------------------------------!
            v1      = 0.5*(v(j) + v(j+1))
            v2      = 0.5*(v(j) + v(j-1))
            Ny(j)   = rdy*(v1**2.0-v2**2.0)
        end do

    end subroutine 
    !=====================================================================!
    ! Pressure calculation                                                !
    !=====================================================================!
    subroutine pressure_calc(P, n_gs, GS_error, M, Pold, vstar, &
                                dy, dt, GS_max, iter_max)
        !-----------------------------------------------------------------!
        ! Pressure calculation preamble                                   !
        !-----------------------------------------------------------------!
        integer,  intent(in)                                :: M
        integer, intent(in)                                 :: iter_max
        real(WP), dimension(0:M), intent(in)                :: vstar 
        real(WP), intent(in)                                :: dy, dt
        real(WP), intent(in)                                :: GS_max
        real(WP), dimension(0:M+1), intent(in)              :: Pold
        real(WP)                                            :: rdy, rdy2, rdt
        real(WP)                                            :: rhsC, rhsC2
        real(WP), dimension(0:M+1), intent(out)             :: P
        integer,  intent(out)                               :: n_gs
        real(WP), intent(out)                               :: GS_error
        real(WP), dimension(0:M+1)                          :: rhs, rhs2
        integer                                             :: j
        real(WP), dimension(0:M+1)                          :: res
        !-----------------------------------------------------------------!
        ! Domain variables                                                !
        !-----------------------------------------------------------------!
        p       = pold
        rdy     = 1.0_WP/(dy)
        rdy2    = 1.0_WP/(dy**2.0_WP)
        rdt     = 1.0_WP/(dt)
        rhsC    = 0.5_WP*dy*rdt
        rhsC2   = rdt*rdy
        !-----------------------------------------------------------------!
        ! Convergence criteria                                            !
        !-----------------------------------------------------------------!
        GS_error    = 1.0_WP
        n_gs        = 1
        !-----------------------------------------------------------------!
        ! Calculating RHS                                                 !
        !-----------------------------------------------------------------!
        do j = 1, M
            rhs(j)  = rhsC*(vstar(j) - vstar(j-1))
            rhs2(j) = rhsC2*(vstar(j)-vstar(j-1))
        end do
        !-----------------------------------------------------------------!
        ! GS loop                                                         !
        !-----------------------------------------------------------------!
        do while (GS_error > GS_max .and. n_gs < iter_max)
            !-------------------------------------------------------------!
            ! Gs calculation                                              !
            !-------------------------------------------------------------!
            do j = 1, M
                P(j)    = 0.5*(P(j+1)-P(j)) - rhs(j)
            end do 
            !-------------------------------------------------------------!
            ! Pressure Neumann boundary conditions                        !
            !-------------------------------------------------------------!
            P(0)    = P(1)
            P(M+1)  = P(M)
            !-------------------------------------------------------------!
            ! Residual calculation                                        !
            !-------------------------------------------------------------!
            do j = 1, M
                res(j)  = rdy2*(P(j+1) - 2.0_WP*P(j) + P(j-1)) - rhs2(j)
            end do
            !-------------------------------------------------------------!
            ! Residual check and iteration update                         !
            !-------------------------------------------------------------!
            GS_error    = maxval(abs(res)) 
            n_gs    = n_gs + 1
        end do
    end subroutine
    !=====================================================================!
    ! u time derivative                                                   !
    !=====================================================================!
    subroutine u_time_derv(unew, M, u, v, dp_dx, dt, dy, nu)
        !-----------------------------------------------------------------!
        ! u time derivative preamble                                      !
        !-----------------------------------------------------------------!
        integer, intent(in)                             :: M
        real(WP), dimension(0:M), intent(in)            :: v
        real(WP), dimension(0:M+1), intent(in)          :: u
        real(WP), intent(in)                            :: dp_dx
        real(WP), intent(in)                            :: dt, dy, nu
        real(WP), dimension(0:M+1), intent(out)         :: unew
        real(WP)                                        :: u1, u2
        real(WP)                                        :: rdy, rdy2
        integer                                         :: j
        !-----------------------------------------------------------------!
        ! Domain variables                                                !
        !-----------------------------------------------------------------!
        rdy     = 1.0_WP/dy
        rdy2    = 1.0_WP/(dy**2.0_WP)
        !-----------------------------------------------------------------!
        ! Looping over domain                                             !
        !-----------------------------------------------------------------!
        do j = 1, M  
            u1      = 0.5_WP*(u(j) + u(j+1)) 
            u2      = 0.5_WP*(u(j) + u(j-1)) 
            unew(j) = u(j)-rdy*dt*(u1*v(j) - u2*v(j-1))  &
                            + rdy2*dt*nu*(u(j+1)-2.0_WP*u(j) + u(j-1)) &
                            + dt*dp_dx
        end do
    end subroutine
    !=====================================================================!
    ! Scalar transport                                                    !
    !=====================================================================!
    subroutine scalar_transport(phinew, M, phiold, v, dt, dy, dif)
        integer, intent(in)                             :: M
        real(WP), dimension(0:M), intent(in)            :: v
        real(WP), dimension(0:M+1), intent(in)          :: phiold
        real(WP), intent(in)                            :: dt, dy, dif
        real(WP), dimension(0:M+1), intent(out)         :: phinew
        real(WP)                                        :: phi1, phi2
        real(WP)                                        :: rdy, rdy2
        integer                                         :: i
        !-----------------------------------------------------------------!
        ! Domain variables                                                !
        !-----------------------------------------------------------------!
        phinew  = 0.0_WP
        rdy     = 1.0_WP/dy
        rdy2    = (1.0_WP/dy)**2.0_WP
        !-----------------------------------------------------------------!
        ! time advancing                                                  !
        !-----------------------------------------------------------------!
        do i = 1, M
            phi1        = 0.5_WP*(phiold(i) + phiold(i+1))    
            phi2        = 0.5_WP*(phiold(i) + phiold(i-1))    
            phinew(i)   = phiold(i) - dt*rdy*(phi1*v(i) - phi2*v(i-1)) + &
                            dif*dt*rdy2*(phiold(i+1) - 2.0_WP*phiold(i) +  phiold(i-1))
            !phinew(i)   = phiold(i) + dif*dt*rdy2*(phiold(i+1) - 2.0_WP*phiold(i) +  phiold(i-1))
        end do
        !-----------------------------------------------------------------!
        ! Neumann boundary conditions                                     !
        !-----------------------------------------------------------------!
        phinew(0)       = 0.0_WP
        phinew(M+1)     = 0.0_WP
    end subroutine
    !=====================================================================!
    ! Scalar transport                                                    !
    !=====================================================================!
    subroutine rk4_density(phinew, M, phiold, v, dt, dy, dif)
        integer, intent(in)                             :: M
        real(WP), dimension(0:M), intent(in)            :: v
        real(WP), dimension(0:M+1), intent(in)          :: phiold
        real(WP), intent(in)                            :: dif
        real(WP), intent(in)                            :: dt, dy
        real(WP), dimension(0:M+1), intent(out)         :: phinew
        real(WP), dimension(1:M)                        :: k1, k2, k3, k4
        real(WP), dimension(0:M+1)                      :: s1, s2, s3
        real(WP)                                        :: phi1, phi2
        real(WP)                                        :: rdy, rdy2
        integer                                         :: i
        !-----------------------------------------------------------------!
        ! Domain variables                                                !
        !-----------------------------------------------------------------!
        phinew  = 0.0_WP
        rdy     = 1.0_WP/dy
        rdy2    = (1.0_WP/dy)**2.0_WP
        !-----------------------------------------------------------------!
        ! Preallocating k's                                               !
        !-----------------------------------------------------------------!
        k1      = 0.0
        k2      = 0.0
        k3      = 0.0
        k4      = 0.0
        !-----------------------------------------------------------------!
        ! Finding K1                                                      !
        !-----------------------------------------------------------------!
        do i = 1, M
            !-------------------------------------------------------------!
            ! K1                                                          !
            !-------------------------------------------------------------!
            phi1    = 0.5_WP*(phiold(i+1) + phiold(i))
            phi2    = 0.5_WP*(phiold(i-1) + phiold(i))
            k1(i)   = -rdy*((phi1*v(i) - phi2*v(i-1))) &
                        + rdy2*dif*(phiold(i+1) -2.0_WP*phiold(i) + phiold(i-1))
        end do
        s1(1:M)     = phiold(1:M) + 0.5_WP*dt*k1  
        s1(0)       = phiold(0)
        s1(M+1)     = phiold(M+1)
        !-----------------------------------------------------------------!
        ! Finding K2                                                      !
        !-----------------------------------------------------------------!
        do i = 1, M
            !-------------------------------------------------------------!
            ! K2                                                          !
            !-------------------------------------------------------------!
            phi1    = (s1(i+1) + s1(i)) 
            phi2    = (s1(i-1) + s1(i)) 
            k2(i)   = -rdy*((phi1*v(i) - phi2*v(i-1))) &
                        + rdy2*dif*(s1(i+1) -2.0_WP*s1(i) + s1(i-1))
        end do
        s2(1:M)     = phiold(1:M) + 0.5_WP*dt*k2  
        s2(0)       = phiold(0)
        s2(M+1)     = phiold(M+1)
        !-----------------------------------------------------------------!
        ! Finding K3                                                      !
        !-----------------------------------------------------------------!
        do i = 1, M
            !-------------------------------------------------------------!
            ! K2                                                          !
            !-------------------------------------------------------------!
            phi1    = (s2(i+1) + s2(i)) 
            phi2    = (s2(i-1) + s2(i)) 
            k3(i)   = -rdy*((phi1*v(i) - phi2*v(i-1)))  &
                        + rdy2*dif*(s2(i+1) -2.0_WP*s2(i) + s2(i-1))
        end do
        s3(1:M)     = phiold(1:M) + 0.5_WP*dt*k3  
        s3(0)       = phiold(0)
        s3(M+1)     = phiold(M+1)
        !-----------------------------------------------------------------!
        ! Finding K4                                                      !
        !-----------------------------------------------------------------!
        do i = 1, M
            !-------------------------------------------------------------!
            ! K4                                                          !
            !-------------------------------------------------------------!
            phi1    = (s3(i+1) + s3(i)) 
            phi2    = (s3(i-1) + s3(i)) 
            k4(i)   = -rdy*((phi1*v(i) - phi2*v(i-1)))  & 
                        + rdy2*dif*(s3(i+1) -2.0_WP*s3(i) + s3(i-1))
        end do 
        !-----------------------------------------------------------------!
        ! Finding phi_new                                                 !
        !-----------------------------------------------------------------!
        phinew(1:M) = phiold(1:M) + (dt/6.0_WP)*(k1 + 2.0_WP*k2 + 2.0_WP*k3 &
                                + k4)
        phinew(0)   = phiold(0)
        phinew(M+1) = phiold(M+1)
    end subroutine
end module


