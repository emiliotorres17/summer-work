program diffusion
    !---------------------------------------------------------------------!
    ! Preamble                                                            !
    !---------------------------------------------------------------------!
    use precision_m
    implicit none
    real(WP), parameter             :: pi = 4.0_WP*atan(1.0_WP)
    integer, parameter              :: n = 256
    real(WP), dimension(0:n)        :: x, y
    real(WP), dimension(0:n,0:n)    :: u, unew
    real(WP)                        :: dx, dy, dt
    integer                         :: i, j, k
    real(WP)                        :: del, p
    integer                         :: tf = 24
    real(WP)                        :: nu
    integer                         :: pcount
    !---------------------------------------------------------------------!
    ! Defining domain variables                                           !
    !---------------------------------------------------------------------!
    nu      = 0.01_WP
    dx      = 5.0_WP/(real(n))
    dy      = 5.0_WP/(real(n))
    dt      = (0.25_WP*dx*dy)
    print *, dx, dy
    !---------------------------------------------------------------------!
    ! Defining x and y array                                              !
    !---------------------------------------------------------------------!
    x(0:n)  = [(i*dx, i=0, n)]
    y(0:n)  = [(i*dy, i=0, n)]
    !---------------------------------------------------------------------!
    ! printing x and y arrays                                             !
    !---------------------------------------------------------------------!
    100 format (1x, 2(f15.10))
    open(unit=2, file='x.dat')
    open(unit=3, file='y.dat')
    do i = 0, n
        write(2,*) x(i)
        write(3,*) y(i)
    end do
    close(unit=2)
    close(unit=3)
    !---------------------------------------------------------------------!
    ! calculating initial conditions                                      !
    !---------------------------------------------------------------------!
    do j = 0,n
        do i = 0,n
            del = sqrt((x(i)-1.5_WP)**2.0_WP + (y(j)-3.0_WP)**2.0_WP)
            if (del <= 1.0_WP) then
                p = cos(0.5_WP*pi*del)
            else
                p = 0.0_WP
            end if
            u(i,j) = p
            print *, u(i,j)
        end do
    end do
    print *, maxval(u)
    open(unit=1, file='IC.dat')
    !---------------------------------------------------------------------!
    ! writing initial conditions                                          !
    !---------------------------------------------------------------------!
    do, i=0,n
        write(1,*) ( u(i,j), j=0,n )
    end do
    close(unit=1)
    !---------------------------------------------------------------------!
    ! time loop                                                          !
    !---------------------------------------------------------------------!
    unew    = u
    pcount  = 0
    print *, tf/dt
    do k = 1, int(tf/dt)
        !-----------------------------------------------------------------!
        ! looping over the domain                                         !
        !-----------------------------------------------------------------!
        do j = 1, n-1
            do i = 1, n-1
                unew(i,j) = u(i,j) + dt*nu*(&
                        (u(i+1,j)-2.0_WP*u(i,j)+u(i-1,j))/dx**2.0&
                        + (u(i,j+1)-2.0_WP*u(i,j)+u(i,j-1))/dy**2.0_WP)
            end do
        end do
        !-----------------------------------------------------------------!
        ! BC at x = 0                                                     !
        !-----------------------------------------------------------------!
        do j = 0, n
            unew(0,j) = 0.0_WP
            unew(n,j) = 0.0_WP
        end do
        !-----------------------------------------------------------------!
        ! BC at y = 0  and y = n                                          !
        !-----------------------------------------------------------------!
        do i = 0, n
            unew(i,0) = 0.0_WP
            unew(i,n) = 0.0_WP
        end do
        u = unew        ! updating solution
        !-----------------------------------------------------------------!
        ! Print statement                                                 !
        !-----------------------------------------------------------------!
        if (pcount == int(1000)) then
            print *, 'time step --> ', k*dt
            pcount = 0
        end if
        pcount = pcount + 1
    end do
    !---------------------------------------------------------------------!
    ! writing final solution                                              !
    !---------------------------------------------------------------------!
    print *, maxval(u)
    open(unit=4, file='data.dat')
    do, i = 0, n
        write(4,*) ( u(i,j), j=0,n )
    end do
    close(unit=4)
end program diffusion
