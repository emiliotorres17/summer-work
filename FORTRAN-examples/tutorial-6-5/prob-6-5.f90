program wave_equation
    !---------------------------------------------------------------------!
    ! Preamble                                                            !
    !---------------------------------------------------------------------!
    use precision_m
    implicit none
    real(WP), parameter             :: pi = 4.0_WP*atan(1.0_WP)
    integer, parameter              :: n = int(400/5)
    real(WP), dimension(0:n)        :: x
    real(WP), dimension(0:n)        :: u, unew
    real(WP)                        :: dx, dt
    integer                         :: i, j, k
    real(WP)                        :: tf = 0.5_WP
    real(WP)                        :: a, c
    !---------------------------------------------------------------------!
    ! Domain variables                                                    !
    !---------------------------------------------------------------------!
    a       = 250.0_WP
    dt      = 0.00001_WP
    dx      = 5.0_WP
    c       = a*dt/dx
    print *, c
    do i = 0, n
        x(i) = i*dx
        print *, x(i)
    end do
    !---------------------------------------------------------------------!
    ! Initial conditions                                                  !
    !---------------------------------------------------------------------!
    do i = 0, n
        if (x(i)>= 0.0 .and. x(i)<= 50.0) then
            u(i) =  0.0_WP
        else if (x(i)>= 50.0 .and. x(i)<=110.0) then
            u(i) = 100.0_WP*(sin(pi*(x(i)-50.0_WP)/60.0_WP))
        else if (x(i) >= 110.0_WP .and. x(i)<=400.0) then
            u(i) = 0.0_WP
        end if
    end do
    !---------------------------------------------------------------------!
    ! Writing initial condition                                           !
    !---------------------------------------------------------------------!
    open(unit=1, file='IC.dat')
    do, i = 0, n
        write(1,*) u(i), x(i)
    end do
    close(unit=1)
    !---------------------------------------------------------------------!
    ! Time stepping                                                       !
    !---------------------------------------------------------------------!
    print *, int(tf/dt)
    do k = 1, int(tf/dt)
        do i = 1, n
            unew(i) =  u(i) -  c*(u(i) - u(i-1))
        end do
        u = unew
    end do
    !---------------------------------------------------------------------!
    ! Writing solution at t = 0.5                                         !
    !---------------------------------------------------------------------!
    open(unit=1, file='c-temp.dat')
    do, i = 0, n
        write(1,*) u(i)
    end do
    close(unit=1)

end program wave_equation
