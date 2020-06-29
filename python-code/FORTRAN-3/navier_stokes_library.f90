module navier_stokes_library
    !=====================================================================!
    ! Dependencies                                                        !
    !=====================================================================!
    use precision_m
    implicit none
    public  :: time_derv_calc, pressure_calc
    contains    
    !=====================================================================!
    ! Time derivate calculation                                           !
    !=====================================================================!
    subroutine time_derv_calc(Nx, Lx, Ny, Ly, M, N, u, v, dx, dy, nu)
        !-----------------------------------------------------------------!
        ! Time derivative preamble                                        !
        !-----------------------------------------------------------------!
        integer,  intent(in)                                :: M, N
        real(WP), intent(in)                                :: nu
        real(WP), dimension(0:N+1, 0:M), intent(in)         :: u 
        real(WP), dimension(0:N, 0:M+1), intent(in)         :: v 
        real(WP), intent(in)                                :: dx, dy
        real(WP)                                            :: rdx, rdy, dx2, dy2
        real(WP), dimension(0:N+1,0:M), intent(out)         :: Nx, Lx 
        real(WP), dimension(0:N,0:M+1), intent(out)         :: Ny, Ly 
        integer                                             :: i ,j
        real(WP)                                            :: u1, u2, u3, u4, v1, v2, v3, v4
        !-----------------------------------------------------------------!
        ! Differencing variables                                          !
        !-----------------------------------------------------------------!
        rdx = 1.0_WP/dx
        rdy = 1.0_WP/dy
        dx2 = 1.0_WP/(dx*dx)
        dy2 = 1.0_WP/(dy*dy)
        !-----------------------------------------------------------------!
        ! Calculating Lx and Nx                                           !
        !-----------------------------------------------------------------!
        do j = 1, N 
            do i = 1, M
                !---------------------------------------------------------!
                ! Calculating Ly                                          !
                !---------------------------------------------------------!
                Lx(j,i) = nu*dx2*(u(j,i+1) + 2.0_WP*u(i,j) + u(i-1,j)) 
                Lx(j,i) = Lx(j,i) + nu*dy2*(u(j+1,i) + 2.0_WP*u(i,j) + u(i,j-1)) 
                !---------------------------------------------------------!
                ! Calculating Nx                                          !
                !---------------------------------------------------------!
                u1      = 0.5_WP*(u(j,i) + u(j-1, i))
                u2      = 0.5_WP*(u(j+1,i) + u(j, i))
                u3      = 0.5_WP*(u(j+1,i) + u(j,i))
                u4      = 0.5_WP*(u(j,i) + u(j-1,i))
                v1      = 0.5_WP*(v(j,i+1) + v(j,i))
                v2      = 0.5_WP*(v(j,i) + v(j-1,i))
                Nx(j,i) = rdx*((u2)**2.0_WP + (u1)**2.0_WP)
                Nx(j,i) = Nx(j,i) + rdx*(v1*u3 + v2*u4)
            end do
        end do
        !-----------------------------------------------------------------!
        ! Calculating Ly and Ny                                           !
        !-----------------------------------------------------------------!
        do j = 1, N 
            do i = 1, M
                !---------------------------------------------------------!
                ! Calculating Ly                                          !
                !---------------------------------------------------------!
                Ly(j,i) = nu*dx2*(v(j,i+1) + 2.0_WP*v(i,j) + v(i-1,j)) 
                Ly(j,i) = Lx(j,i) + nu*dy2*(v(j+1,i) + 2.0_WP*v(i,j) + v(i,j-1)) 
                !---------------------------------------------------------!
                ! Calculating Ny                                          !
                !---------------------------------------------------------!
                u1      = 0.5_WP*(u(j+1,i) + u(j,i))
                v1      = 0.5_WP*(v(j,i+1) + v(j,i))
                u2      = 0.5_WP*(u(j-1,i) + u(j,i))
                v2      = 0.5_WP*(v(j,i-1) + v(j,i))
                v3      = 0.5_WP*(v(j+1,i) + v(j,i))
                v4      = 0.5_WP*(v(j-1,i) + v(j,i))
                Ny(j,i) = rdx*(u1*v1 - u2*v2)
                Ny(j,i) = Ny(j,i) + rdy*(v3**2.0_WP - v4**2.0_WP)
            end do
        end do
    end subroutine
    !=====================================================================!
    ! Pressure calculation                                                !
    !=====================================================================!
    subroutine pressure_calc(P, M, N, ustar, vstar, dx, dy, dt, GS_max)
        !-----------------------------------------------------------------!
        ! Pressure calculation preamble                                   !
        !-----------------------------------------------------------------!
        integer,  intent(in)                                :: M, N
        real(WP), dimension(0:N+1, 0:M), intent(in)         :: ustar 
        real(WP), dimension(0:N, 0:M+1), intent(in)         :: vstar 
        real(WP), intent(in)                                :: dx, dy, dt
        real(WP), intent(in)                                :: GS_max
        real(WP)                                            :: rdx, rdy, rdt, dx2, dy2
        real(WP), dimension(0:N+1, 0:M+1), intent(out)      :: P
        real(WP), dimension(0:N+1, 0:M+1)                   :: f
        integer                                             :: i ,j
        real(WP)                                            :: GS_error
        real(WP), dimension(1:M, 1:N)                       :: res
        integer                                             :: n_gs
        !-----------------------------------------------------------------!
        ! Differentiation variables                                       !
        !-----------------------------------------------------------------!
        dx2     = 1.0_WP/(dx*dx)
        dy2     = 1.0_WP/(dy*dy)
        rdx     = 1.0_WP/dx
        rdy     = 1.0_WP/dy
        rdt     = 1.0_WP/dt
        !-----------------------------------------------------------------!
        ! Calculate RHS                                                   !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i = 1, M
                f(j,i)  = rdx*(ustar(j,i) - ustar(j,i-1)) + &
                          rdy*(vstar(j,i) - ustar(j-1,i))
            end do
        end do
        !-----------------------------------------------------------------!
        ! GS iteration                                                    !
        !-----------------------------------------------------------------!
        do while (GS_error > GS_max)
            n_gs    = n_gs + 1
            !-------------------------------------------------------------!
            ! Calculating Pressure                                        !
            !-------------------------------------------------------------!
            do j = 1, N
                do i = 1, M
                    P(j,i) = 0.25_WP*(P(j,i+1) + P(j,i-1) + P(j+1,i) + P(j-1, i)) - &
                                0.25_WP*dx2*rdt*f(j,i)
                end do
            end do
            !-------------------------------------------------------------!
            ! Applying periodic boundary conditions                       !
            !-------------------------------------------------------------!
            P(:,0)      = P(:,1)
            P(:,M+1)    = P(:,M)
            P(0,:)      = P(1,:)
            P(M+1,:)    = P(M,:)
            !-------------------------------------------------------------!
            ! Calculating Pressure residual                               !
            !-------------------------------------------------------------!
            do j = 1, N
                do i = 1, M
                    res(j,i) = dx2*(P(j,i+1) - 2.0_WP*P(j,i) + P(j,i-1)) + &
                            dy2*(P(j+1,i) - 2.0_WP*P(j,i) + P(j-1,i)) - &
                            rdt*f(j,i)
                end do
            end do
            !-------------------------------------------------------------!
            ! Checking convergence                                        !
            !-------------------------------------------------------------!
            GS_error    = maxval(abs(res))
        end do
    end subroutine
end module