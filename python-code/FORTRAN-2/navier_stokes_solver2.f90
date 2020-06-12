program navier_stokes
    use precision_m    
    use ns_lib 
    real(WP)                                        :: xs = 0.0_WP, xf = 1.0_WP
    integer, parameter                              :: xlen = 101
    integer, parameter                              :: M = 16, N = 16
    real(WP), dimension(1:xlen)                     :: xvec
    real(WP)                                        :: dx, dy, nu
    real(WP), dimension(0:N)                        :: ul, ur, vl, vr
    real(WP), dimension(0:M)                        :: ub, ut, vb, vt
    real(WP), dimension(0:N, 0:M-1)                 :: u, ustar
    real(WP), dimension(0:N-1,0:M)                  :: v, vstar
    real(WP), dimension(0:N, 0:M-1)                 :: Lu
    real(WP), dimension(0:N, 0:M-1)                 :: Nx
    real(WP), dimension(0:N-1, 0:M)                 :: Lv, Ny
    real(WP),dimension(0:N, 0:M)                    :: P
    real(WP)                                        :: maxerror
    real(WP), dimension(0:N-1,0:M-1)                :: uinterp
    real(WP)                                        :: velcent 
    integer                                         :: i, j, nstep, tstep_F
    real(WP)                                        :: tmax, t
    integer                                         :: counter
    real(WP)                                        :: max_div
    real(WP)                                        :: div
    !---------------------------------------------------------------------!
    ! Preallocating variables                                             !
    !---------------------------------------------------------------------!
    counter     = 1
    u           = 0.0_WP
    v           = 0.0_WP
    ustar       = 0.0_WP
    vstar       = 0.0_WP
    !---------------------------------------------------------------------!
    ! Boundary conditions                                                 !
    !---------------------------------------------------------------------!
    ut          = 1.0_WP
    ub          = 0.0_WP
    ur          = 0.0_WP
    ul          = 0.0_WP
    vt          = 0.0_WP
    vb          = 0.0_WP
    vr          = 0.0_WP
    vl          = 0.0_WP
    !---------------------------------------------------------------------!
    ! Transport properties                                                !
    !---------------------------------------------------------------------!
    nu          = 0.01_WP
    dx          = 1.0_WP/real(M)
    dy          = 1.0_WP/real(N)
    dt          = 0.25_WP*dx*dy/nu/2.0_WP
    tmax        = 5.0_WP
    maxerror    = (10.0_WP)**(-12.0_WP)
    tstep_F     = int(1.0000001_WP*tmax/dt)
    print *, tstep_F
    !---------------------------------------------------------------------!
    ! Print variables                                                     !
    !---------------------------------------------------------------------!
    open(unit=1, file='FORTRAN-data/output.out')
    10 format(I8, A, I8, 4X, A, f10.6, 4X, A, ES16.5, /, I8, A, I8,&
                4X, A, ES16.5)
    !---------------------------------------------------------------------!
    ! Time stepping loop                                                  !
    !---------------------------------------------------------------------!
    do nstep = 1, tstep_F
        t   = t + dt   
        call time_derv(Lu, Lv, Nx, Ny, M, N, dx, dy, nu, u, v, ul, ur, ub, ut, &
                                    vl, vr, vb, vt)
        !-----------------------------------------------------------------!
        ! Updating ustar and vstar                                        !
        !-----------------------------------------------------------------!
        do j = 1, N
            do i  = 1, M-1
                ustar(j,i) = u(j,i) + dt*(-Nx(j,i) + Lu(j,i))
            end do
        end do
        do j = 1, N-1
            do i = 1, M
                vstar(j,i) = v(j,i) + dt*(-Ny(j,i) + Lv(j,i))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Updating pressure                                               !
        !-----------------------------------------------------------------!
        call calcpress(P, M, N, dx, dy, maxerror, ustar, vstar, ul, &
                                ur, ub, ut, vl, vr, vb, vt)
        !-----------------------------------------------------------------!
        ! Updating u and v                                                !
        !-----------------------------------------------------------------!
        do j  = 1, N
            do i = 1, M-1
                u(j,i) = ustar(j,i) - (P(j,i+1)-P(j,i))/dx
            end do
        end do 
        do j  = 1, N-1
            do i  = 1, M
                v(j,i) = vstar(j,i) - (P(j+1,i)-P(j,i))/dy
            end do 
        end do
        !-----------------------------------------------------------------!
        ! Calculating velocity divergence                                 !
        !-----------------------------------------------------------------!
        call divergence(max_div, div, M, N, u, v, ul, ur, ut, ub, &
                            vl, vr, vt, vb, dx, dy)
        !-----------------------------------------------------------------!
        ! Print statement                                                 !
        !-----------------------------------------------------------------!
        print '(I8, A, I8, 4X, A, f10.6, 4X, A, ES10.3)', &
                nstep, '/', tstep_F, &
                'time -->', t, &
                'max velocity -->', max(maxval(u),maxval(v))
        !-----------------------------------------------------------------!
        ! Divergence print out                                            !
        !-----------------------------------------------------------------!
        print '(I8, A , I8, A, ES16.5)', &
            nstep, '/', tstep_F, &
            '| Maximum velocity diveregence:', max_div
        !-----------------------------------------------------------------!
        ! Writing output                                                  !
        !-----------------------------------------------------------------!
        write(1,  10)&
                nstep, '/', tstep_F, &
                'time -->', t, &
                'max velocity -->', max(maxval(u),maxval(v)), &
                nstep, '/', tstep_F, &
                '| Maximum velocity diveregence:', max_div
        !-----------------------------------------------------------------!
        ! Writing u-velocity                                              !
        !-----------------------------------------------------------------!
        !-------------------------------------------------------------!
        ! Cell centered u-velocity                                    !
        !-------------------------------------------------------------!
        call u_write(u,M,N,ul,ur,ut,ub)
        call v_write(v,M,N,vl,vr,vt,vb)
    end do
    close(unit=1)
end program navier_stokes
