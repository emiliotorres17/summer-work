program SOR_example
    !=====================================================================!
    ! Purpose:                                                            !
    !   The purpose of this script is to  practice with the SOR method.   !
    !                                                                     !
    ! Author:                                                             !
    !   Emilio Torres                                                     !
    !=====================================================================!
    !=====================================================================!
    ! Preamble                                                            !
    !=====================================================================!
    use precision_m                                       
    implicit none
    integer, parameter                      :: itermax = 5000
    real(WP), parameter                     :: pi = 4.0_WP*atan(1.0_WP) 
    integer, parameter                      :: M = 512, N = 512
    real(WP), dimension(0:M)                :: x
    real(WP), dimension(0:N)                :: y
    real(WP), dimension(0:N, 0:M)           :: phi, res, f
    real(WP), dimension(1:itermax)          :: res_vec
    real(WP)                                :: omega, dx, dy, rho
    integer                                 :: i, j, k, counter
    !=====================================================================!
    ! Output files and formats                                            !
    !=====================================================================!
    open(unit=1, file='data/SOR-IC.dat')
    open(unit=2, file='data/SOR-data.dat')
    open(unit=3, file='data/residual-data.dat')
    10 format(550ES25.12)
    11 format(ES25.12)
    !=====================================================================!
    ! Spatial and time stepping variables                                 !
    !=====================================================================!
    dx      = 2.0_WP/dble(M)
    dy      = 2.0_WP/dble(N)
    rho     = 0.5_WP*(cos(pi/M) + cos(pi/N))
    omega   = 2.0_WP/(1.0_WP + sqrt(1.0_WP -  rho**(2.0_WP)))
    print *, omega
    !=====================================================================!
    ! Preallocating variables                                             !
    !=====================================================================!
    f       = 0.0_WP
    phi     = 0.0_WP
    res     = 0.0_WP
    res_vec = 0.0_WP
    counter = 500
    !=====================================================================!
    ! Calculating x and y vectors                                         !
    !=====================================================================!
    x(0)    = -1.0_WP
    y(0)    = -1.0_WP
    do i = 1, M 
        x(i) = x(i-1) + dx
    end do
    do j = 1, N 
        y(j) = y(j-1) + dx
    end do
    !=====================================================================!
    ! Calculating RHS                                                     !
    !=====================================================================!
    do j = 0, N
        do i = 0, M
            f(j,i) = 10.0_WP - 10.0_WP*(cos(x(i)))**2.0_WP + 10.0_WP*sin(y(j))
        end do
    end do
    !=====================================================================!
    ! Setting the initial condition                                       !
    !=====================================================================!
    do j = 0, N
        do i = 0, M
            phi(j,i)    = 0.5_WP*sin(pi*x(i)) * sin(4.0_WP*pi*x(i))
        end do
    end do
    phi(:,0)    = 0.0_WP
    phi(:,M)    = 0.0_WP
    phi(0,:)    = 0.0_WP
    phi(N,:)    = 0.0_WP
    !=====================================================================!
    ! Writing initial condition                                           !
    !=====================================================================!
    do j = 0, N
        write(1, 10) (phi(j,i), i = 0, M)
    end do
    !=====================================================================!
    ! Iteration loop                                                      !
    !=====================================================================!
    do k = 1, itermax
        !-----------------------------------------------------------------!
        ! SOR method                                                      !
        !-----------------------------------------------------------------!
        do j = 1, N-1
            do i = 1, M-1
                phi(j,i)    = phi(j,i) + &
                                omega*(0.25_WP*(phi(j,i+1) + phi(j,i-1) + phi(j+1,i) + phi(j-1,i)) &
                                - 0.25_WP*dx**2.0_WP*f(j,i) - phi(j,i)) 
            end do
        end do
        !-----------------------------------------------------------------!
        ! Maximum residual for each iteration                             !
        !-----------------------------------------------------------------!
        do j = 1, N-1
            do i = 1, M-1
                res(j,i)    = (phi(j,i+1) -2.0_WP*phi(j,i) + phi(j,i-1))/(dx*dx) + &
                                (phi(j+1,i) -2.0_WP*phi(j,i) + phi(j-1,i))/(dy*dy) - & 
                                f(j,i)
            end do
        end do
        res_vec(k)  = maxval(abs(res))
        counter     = counter + 1 
        if (counter > 100) then
            print *, res_vec(k)
            counter = 0
        end if 
    end do
    !=====================================================================!
    ! Writing initial condition                                           !
    !=====================================================================!
    do j = 0, N
        write(2, 10) (phi(j,i), i = 0, M)
    end do
    !=====================================================================!
    ! Writing residual iteration data                                     !
    !=====================================================================!
    do j = 1, itermax
        write(3, 11)    res_vec(j)
    end do 
        

    end program SOR_example
