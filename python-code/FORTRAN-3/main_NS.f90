program main_NS
    !=====================================================================!
    ! Main preamble                                                       !
    !=====================================================================!
    use precision_m
    use navier_stokes_library   
    implicit none
    integer, parameter                      :: M = 32, N = 32
    real(WP)                                :: dx, dy, dt
    real(WP), dimension(0:N+1, 0:M)         :: u, ustar 
    real(WP), dimension(0:N, 0:M+1)         :: v, vstar
    real(WP), dimension(0:N+1, 0:M+1)       :: P
    real(WP)                                :: max_error
    real 
    !=====================================================================!
    ! Domain variables                                                    !
    !=====================================================================!
