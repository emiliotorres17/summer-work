module navier_stokes_library
    !=====================================================================!
    ! Dependencies                                                        !
    !=====================================================================!
    use precision_m
    implicit none
    public  :: time_derv_calc, pressure_calc, velocity_boundary_condition, &
                velocity_write
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
        ! Preallocating output variables                                  !
        !-----------------------------------------------------------------!
        Lx  = 0.0_WP
        Nx  = 0.0_WP
        Ly  = 0.0_WP
        Ny  = 0.0_WP
        !-----------------------------------------------------------------!
        ! Calculating Lx and Nx                                           !
        !-----------------------------------------------------------------!
        do j = 1, N 
            do i = 0, M
                !---------------------------------------------------------!
                ! Left boundary                                           !
                !---------------------------------------------------------!
                if (i == 0) then
                    !-----------------------------------------------------!
                    ! Calculating Ly                                      !
                    !-----------------------------------------------------!
                    Lx(j,i) = nu*dx2*(u(j,i+1) - 2.0_WP*u(j,i) + u(j,M)) 
                    Lx(j,i) = Lx(j,i) + nu*dy2*(u(j+1,i) - 2.0_WP*u(j,i) + u(j-1,i)) 
                    !-----------------------------------------------------!
                    ! Calculating Nx                                      !
                    !-----------------------------------------------------!
                    u1      = 0.5_WP*(u(j,i)    + u(j, M))
                    u2      = 0.5_WP*(u(j,i+1)  + u(j, i))
                    u3      = 0.5_WP*(u(j+1,i)  + u(j,i))
                    u4      = 0.5_WP*(u(j,i)    + u(j-1,i))
                    v1      = 0.5_WP*(v(j,i+1)      + v(j,i))
                    v2      = 0.5_WP*(v(j-1,i+1)    + v(j-1,i))
                    Nx(j,i) = rdx*((u2)**2.0_WP - (u1)**2.0_WP)
                    Nx(j,i) = Nx(j,i) + rdx*(v1*u3 - v2*u4)
                !---------------------------------------------------------!
                ! Right boundary                                          !
                !---------------------------------------------------------!
                elseif ( i == M) then
                    !-----------------------------------------------------!
                    ! Calculating Ly                                      !
                    !-----------------------------------------------------!
                    Lx(j,i) = nu*dx2*(u(j,0) - 2.0_WP*u(j,i) + u(j,i-1)) 
                    Lx(j,i) = Lx(j,i) + nu*dy2*(u(j+1,i) - 2.0_WP*u(j,i) + u(j-1,i)) 
                    !-----------------------------------------------------!
                    ! Calculating Nx                                      !
                    !-----------------------------------------------------!
                    u1      = 0.5_WP*(u(j,i)    + u(j, i-1))
                    u2      = 0.5_WP*(u(j,0)  + u(j, i))
                    u3      = 0.5_WP*(u(j+1,i)  + u(j,i))
                    u4      = 0.5_WP*(u(j,i)    + u(j-1,i))
                    v1      = 0.5_WP*(v(j,i+1)      + v(j,i))
                    v2      = 0.5_WP*(v(j-1,i+1)    + v(j-1,i))
                    Nx(j,i) = rdx*((u2)**2.0_WP - (u1)**2.0_WP)
                    Nx(j,i) = Nx(j,i) + rdx*(v1*u3 - v2*u4)
                !---------------------------------------------------------!
                ! Interior points                                         !
                !---------------------------------------------------------!
                else
                    !-----------------------------------------------------!
                    ! Calculating Ly                                      !
                    !-----------------------------------------------------!
                    Lx(j,i) = nu*dx2*(u(j,i+1) - 2.0_WP*u(j,i) + u(j,i-1)) 
                    Lx(j,i) = Lx(j,i) + nu*dy2*(u(j+1,i) - 2.0_WP*u(j,i) + u(j-1,i)) 
                    !-----------------------------------------------------!
                    ! Calculating Nx                                      !
                    !-----------------------------------------------------!
                    u1      = 0.5_WP*(u(j,i)    + u(j, i-1))
                    u2      = 0.5_WP*(u(j,i+1)  + u(j, i))
                    u3      = 0.5_WP*(u(j+1,i)  + u(j,i))
                    u4      = 0.5_WP*(u(j,i)    + u(j-1,i))
                    v1      = 0.5_WP*(v(j,i+1)      + v(j,i))
                    v2      = 0.5_WP*(v(j-1,i+1)    + v(j-1,i))
                    Nx(j,i) = rdx*((u2)**2.0_WP - (u1)**2.0_WP)
                    Nx(j,i) = Nx(j,i) + rdx*(v1*u3 - v2*u4)
                end if 
            end do
        end do
        !-----------------------------------------------------------------!
        ! Calculating Ly and Ny                                           !
        !-----------------------------------------------------------------!
        do j = 1, N 
            do i = 0, M
                !---------------------------------------------------------!
                ! Left boundary                                           !
                !---------------------------------------------------------!
                if ( i == 0) then
                    !-----------------------------------------------------!
                    ! Calculating Ly                                      !
                    !-----------------------------------------------------!
                    Ly(j,i) = nu*dx2*(v(j,i+1) - 2.0_WP*v(j,i) + v(j,M)) 
                    Ly(j,i) = Ly(j,i) + nu*dy2*(v(j+1,i) - 2.0_WP*v(j,i) + v(j-1,i)) 
                    !-----------------------------------------------------!
                    ! Calculating Ny                                      !
                    !-----------------------------------------------------!
                    u1      = 0.5_WP*(u(j+1,i) + u(j,i))
                    v1      = 0.5_WP*(v(j,i+1) + v(j,i))
                    u2      = 0.5_WP*(u(j+1,M) + u(j,M))
                    v2      = 0.5_WP*(v(j,M)   + v(j,i))
                    v3      = 0.5_WP*(v(j+1,i) + v(j,i))
                    v4      = 0.5_WP*(v(j-1,i) + v(j,i))
                    Ny(j,i) = rdx*(u1*v1 - u2*v2)
                    Ny(j,i) = Ny(j,i) + rdy*(v3**2.0_WP - v4**2.0_WP)
                !---------------------------------------------------------!
                ! Interior points                                         !
                !---------------------------------------------------------!
                else
                    !-----------------------------------------------------!
                    ! Calculating Ly                                      !
                    !-----------------------------------------------------!
                    Ly(j,i) = nu*dx2*(v(j,i+1) - 2.0_WP*v(j,i) + v(j,i-1)) 
                    Ly(j,i) = Ly(j,i) + nu*dy2*(v(j+1,i) - 2.0_WP*v(j,i) + v(j-1,i)) 
                    !-----------------------------------------------------!
                    ! Calculating Ny                                      !
                    !-----------------------------------------------------!
                    u1      = 0.5_WP*(u(j+1,i) + u(j,i))
                    v1      = 0.5_WP*(v(j,i+1) + v(j,i))
                    u2      = 0.5_WP*(u(j+1,i-1) + u(j,i-1))
                    v2      = 0.5_WP*(v(j,i-1) + v(j,i))
                    v3      = 0.5_WP*(v(j+1,i) + v(j,i))
                    v4      = 0.5_WP*(v(j-1,i) + v(j,i))
                    Ny(j,i) = rdx*(u1*v1 - u2*v2)
                    Ny(j,i) = Ny(j,i) + rdy*(v3**2.0_WP - v4**2.0_WP)
                end if 
            end do
        end do
    end subroutine
    !=====================================================================!
    ! Pressure calculation                                                !
    !=====================================================================!
    subroutine pressure_calc(P, n_gs, GS_error, M, N, Pold, ustar, vstar, &
                                dx, dy, dt, GS_max, iter_max)
        !-----------------------------------------------------------------!
        ! Pressure calculation preamble                                   !
        !-----------------------------------------------------------------!
        integer,  intent(in)                                :: M, N
        integer, intent(in)                                 :: iter_max
        real(WP), dimension(0:N+1, 0:M), intent(in)         :: ustar 
        real(WP), dimension(0:N, 0:M+1), intent(in)         :: vstar 
        real(WP), intent(in)                                :: dx, dy, dt
        real(WP), intent(in)                                :: GS_max
        real(WP), dimension(0:N+1, 0:M+1), intent(in)       :: Pold
        real(WP)                                            :: rdx, rdy, rdt, dx2, dy2
        real(WP), dimension(0:N+1, 0:M+1), intent(out)      :: P
        integer,  intent(out)                               :: n_gs
        real(WP), intent(out)                               :: GS_error
        real(WP), dimension(0:N+1, 0:M+1)                   :: f, f2
        integer                                             :: i ,j
        real(WP), dimension(0:M+1, 0:N+1)                   :: res
        integer                                             :: counter = 0
        real(WP)                                            :: rho, omega
        real(WP)                                            :: pi = 4.0_WP*atan(1.0_WP)
        !-----------------------------------------------------------------!
        ! Differentiation variables                                       !
        !-----------------------------------------------------------------!
        dx2     = 1.0_WP/(dx*dx)
        dy2     = 1.0_WP/(dy*dy)
        rdx     = 1.0_WP/dx
        rdy     = 1.0_WP/dy
        rdt     = 1.0_WP/dt
        f       = 0.0_WP
        f2      = 0.0_WP
        P       = Pold
        rho     = 0.5_WP*(cos(pi/dble(M)) + cos(pi/dble(N)))
        omega   = 2.0_WP/(1.0_WP + sqrt(1.0_WP -  rho**(2.0_WP)))
        !omega   = 1.0_WP
        res     = 0.0_WP
        !-----------------------------------------------------------------!
        ! Calculate RHS                                                   !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i = 0, M
                if ( i == 0) then  
                    f(j,i)  = 0.25_WP*(dx/dt)*(ustar(j,i) - ustar(j,M) + &
                                     vstar(j,i) - vstar(j-1,i))
                else
                    f(j,i)  = 0.25_WP*(dx/dt)*(ustar(j,i) - ustar(j,i-1) + &
                                     vstar(j,i) - vstar(j-1,i))
                end if 
                if ( i == 0) then
                    f2(j,i) = rdt*rdx*(ustar(j,i) - ustar(j,M) + &
                                     vstar(j,i) - vstar(j-1,i))
                else
                    f2(j,i) = rdt*rdx*(ustar(j,i) - ustar(j,i-1) + &
                                     vstar(j,i) - vstar(j-1,i))
                end if
            end do
        end do
        !-----------------------------------------------------------------!
        ! GS iteration                                                    !
        !-----------------------------------------------------------------!
        n_gs        = 0
        GS_error    = 100.0_WP
        do while (GS_error > GS_max .and. n_gs < iter_max)
            n_gs    = n_gs + 1
            !-------------------------------------------------------------!
            ! Calculating Pressure                                        !
            !-------------------------------------------------------------!
            do j = 1, N
                do i = 0, M
                    if ( i == 0) then
                        P(j,i) = P(j,i) + & 
                                    omega*(0.25_WP*(P(j,i+1) + P(j,M) + &
                                    P(j+1,i) + P(j-1, i)) - &
                                    f(j,i) - P(j,i))
                    else
                        P(j,i) = P(j,i) + & 
                                    omega*(0.25_WP*(P(j,i+1) + P(j,i-1) + &
                                    P(j+1,i) + P(j-1, i)) - &
                                    f(j,i) - P(j,i))
                    end if
                end do
            end do
            !-------------------------------------------------------------!
            ! Applying periodic boundary conditions                       !
            !-------------------------------------------------------------!
            !P(:,0)      = P(:,M+1)
            P(:,M+1)    = P(:,0)
            P(0,:)      = P(1,:)
            P(N+1,:)    = P(N,:)
            !-------------------------------------------------------------!
            ! Calculating Pressure residual                               !
            !-------------------------------------------------------------!
            do j = 1, N
                do i = 0, M
                    if (i==0) then
                        res(j,i) = dx2*(P(j,i+1) - 2.0_WP*P(j,i) + P(j,M)) + &
                                    dy2*(P(j+1,i) - 2.0_WP*P(j,i) + P(j-1,i)) - &
                                    f2(j,i)
                    else
                        res(j,i) = dx2*(P(j,i+1) - 2.0_WP*P(j,i) + P(j,i-1)) + &
                                dy2*(P(j+1,i) - 2.0_WP*P(j,i) + P(j-1,i)) - &
                                f2(j,i)
                    end if 
                end do
            end do
            !-------------------------------------------------------------!
            ! Checking convergence                                        !
            !-------------------------------------------------------------!
            GS_error    = maxval(abs(res))
            counter     = counter + 1 
        end do
    end subroutine
    !=====================================================================!
    ! Velocity Boundary conditions                                        !
    !=====================================================================!
    subroutine velocity_boundary_condition(velx, vely, M, N, u, v, ul, ur, &
                                            vl, vr, ut, ub, vt, vb) 
        !-----------------------------------------------------------------!
        ! Boundary condition preamble                                     !
        !-----------------------------------------------------------------!
        integer, intent(in)                                 :: M, N
        real(WP), intent(in)                                :: ul ,ur, ut, ub, vl, vr, vt, vb
        real(WP), dimension(0:N+1, 0:M), intent(in)         :: u
        real(WP), dimension(0:N, 0:M+1), intent(in)         :: v
        real(WP), dimension(0:N+1, 0:M), intent(out)        :: velx
        real(WP), dimension(0:N, 0:M+1), intent(out)        :: vely
        !-----------------------------------------------------------------!
        ! Setting u-x boundary conditions                                 !
        !-----------------------------------------------------------------!
        velx(:,0)   = ul  
        velx(:,M)   = ur 
        velx(0,:)   = 2.0_WP*ub - u(1,:)
        velx(N+1,:) = 2.0_WP*ut - u(N,:)
        !-----------------------------------------------------------------!
        ! Setting u-y boundary conditions                                 !
        !-----------------------------------------------------------------!
        vely(:,0)   = 2.0_WP*vl - v(:,1)  
        vely(:,M+1) = 2.0_WP*vr - v(:,M)  
        vely(0,:)   = vb
        vely(N,:)   = vt
    end subroutine
    !=====================================================================!
    ! Interpolated velocity fields                                        !
    !=====================================================================!
    subroutine velocity_write(M, N, u, v, flag)
        !-----------------------------------------------------------------!
        ! Velocity write preamble                                         !
        !-----------------------------------------------------------------!
        integer, intent(in)                         :: M, N
        real(WP), dimension(0:N+1, 0:M)             :: u
        real(WP), dimension(0:N, 0:M)               :: v
        integer, intent(in)                         :: flag
        real(WP), dimension(1:N, 1:M)               :: ucent
        integer                                     :: i, j
        !-----------------------------------------------------------------!
        ! Averaging u-velocity                                            !
        !-----------------------------------------------------------------!
        if (flag == 0) then
            open(unit=1, file='u-velocity.dat')
            10 format(64ES25.10)
            do j = 1, N
                do i = 1, M
                    ucent(j,i) = 0.5_WP*(u(j,i+1) + u(j,i))  
                end do
            end do
            do j = 1, N
                write(1, 10) (u(j,i), i=1, M)
            end do
            close(unit=1)
        !-----------------------------------------------------------------!
        ! Averaging v-velocity                                            !
        !-----------------------------------------------------------------!
        else if (flag == 1) then
            open(unit=2, file='v-velocity.dat')
            11 format(64ES25.10)
            do j = 1, N
                do i = 1, M
                    ucent(j,i) = 0.5_WP*(v(j+1,i) + v(j,i))  
                end do
            end do
            do j = 1, N
                write(2, 11) (v(j,i), i=1, M)
            end do
            close(unit=2)
        end if
        end subroutine 
end module
