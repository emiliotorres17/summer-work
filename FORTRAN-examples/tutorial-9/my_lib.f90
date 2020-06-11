module my_lib
    !=====================================================================!
    ! Preamble                                                            !
    !=====================================================================!
    use precision_m
    implicit none
    public :: linspace
    !=====================================================================!
    ! Subroutines                                                         !
    !=====================================================================!
    contains
        !=================================================================!
        ! Creates a 1-D dimensional array with evenly spaced elements     !
        !=================================================================!
        subroutine linspace(x, xstart, xend, xlen)
            !-------------------------------------------------------------!
            ! linspace preamble                                           !
            !-------------------------------------------------------------!
            real(WP), dimension(:), intent(out) :: x
            real(WP)                            :: xstart, xend, dx
            integer                             :: xlen, i       
            !-------------------------------------------------------------!
            ! Generating x output                                         !
            !-------------------------------------------------------------!
            dx          = (xend - xstart)/(xlen-1)              ! spatial step size
            x(1:xlen)   = [(xstart + ((i-1)*dx), i=1, xlen)]    ! x array
        end subroutine

end module my_lib
