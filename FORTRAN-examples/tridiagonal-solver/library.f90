module library
    use precision_m
    implicit none
    public :: tri_solve
    contains 
    !=====================================================================!
    ! Tridiagnol solver                                                   !
    !=====================================================================!
    subroutine tri_solve(sol, M, vec3, vec1, vec2, rhs)
        !-----------------------------------------------------------------!
        ! Time derivative preamble                                        !
        !-----------------------------------------------------------------!
        integer,  intent(in)                                :: M
        real(WP), dimension(1:M-1), intent(in)              :: vec1, rhs
        real(WP), dimension(1:M-1), intent(out)             :: sol
        real(WP), dimension(2:M-1), intent(in)              :: vec3
        real(WP), dimension(1:M-2), intent(in)              :: vec2
        real(WP), dimension(1:M-1)                          :: b, d
        real(WP), dimension(2:M-1)                          :: a
        real(WP), dimension(1:M-2)                          :: c
        integer                                             :: i
        !-----------------------------------------------------------------!
        ! Defining variables                                              !
        !-----------------------------------------------------------------!
        a   = vec3
        b   = vec1
        c   = vec2
        d   = rhs
        !-----------------------------------------------------------------!
        ! Elimination                                                     !
        !-----------------------------------------------------------------!
        do i = 2, M-1
            b(i)    = b(i) - c(i-1)*a(i)/b(i-1)
            d(i)    = d(i) - d(i-1)*a(i)/b(i-1)
        end do
        !-----------------------------------------------------------------!
        ! Back substitution                                               !
        !-----------------------------------------------------------------!
        d(M-1)    = d(M-1)/b(M-1)
        do i = M-2, 1, -1
            d(i)    = (d(i)-c(i)*d(i+1))/b(i)
        end do
        sol = d

    end subroutine 
end module library
