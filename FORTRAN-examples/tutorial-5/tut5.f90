program array
    !---------------------------------------------------------------------!
    ! Preamble                                                            !
    !---------------------------------------------------------------------!
    implicit none
    real,       parameter       :: pi = 4d0*atan(1d0)
    integer,    parameter       :: n = 100
    real,       dimension(1:n)  :: x,y
    real                        :: a = 0d0, b = 2d0*pi
    real                        :: increment
    integer                     :: i
    !---------------------------------------------------------------------!
    ! Calculating the increment                                           !
    !---------------------------------------------------------------------!
    increment   = (b-a)/(real(n)-1)
    !---------------------------------------------------------------------!
    ! Defining x and y                                                    !
    !---------------------------------------------------------------------!
    x(1)    = 0.0
    do i = 2,n
        x(i) =  x(i-1) + increment
    end do 
    y   = sin(x)
    !---------------------------------------------------------------------!
    ! Writing output                                                      !
    !---------------------------------------------------------------------!
    open(unit=1, file='data.dat')
    do i = 1,n
        write(1,*)  x(i), y(i)
    end do
    close(unit=1)


end program 