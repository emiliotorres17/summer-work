program tutorial
    !=====================================================================!
    ! Preamble                                                            !
    !=====================================================================!
    use precision_m
    use my_lib
    implicit none
    integer, parameter              :: len_x = 10
    real(WP), dimension(1:len_x)    :: x,y
    integer                         :: i
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
    call linspace(x, 1.0_WP, 10.0_WP, len_x)
    y   = x**(2.0_WP)

    do i = 1, len_x
        write(*, '(f5.1, (a), f5.1)') x(i), ' ', y(i)
    end do



end program tutorial
