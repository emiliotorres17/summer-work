program ns_solve
    !=====================================================================!
    ! Preamble                                                            !
    !=====================================================================!
    use precision_m
    implicit none
    integer, parameter              :: M = 128, N = 128
    integer                         :: n_count, n_gs, n_gs_t, i, j
    real(WP)                        :: dx, dy, dt1, dt2, dt
    real(WP)                        :: x, y, t
    real(WP)                        :: tfinal
    real(WP)                        :: Re, L, nu, Uwall
    real(WP)                        :: Urx, Ulx, Utx, Ubx
    real(WP)                        :: Ury, Uly, Uty, Uby
    real(WP)                        :: conv, conv_crit, conv_gs, conv_gs_limit
    real(WP)                        :: conv_u, conv_v
    real(WP)                        :: rh, rh2, rRe
    real(WP)                        :: a, b, sf
    real(WP)                        :: ua, ub, uc, ud, ue, va, vb, vc, vd, ve, vf
    real(WP), dimension(200000)     :: conv_hist_u, conv_hist_v
    real(WP), dimension(M+1, N+2)   :: u, u_star, a_grid, r_u, res_u
    real(WP), dimension(M+2, N+1)   :: v, v_star, b_grid, r_v, res_v
    real(WP), dimension(M+2, N+2)   :: phi, res_gs, gs_RHS, gs_RHS2
    real(WP), dimension(M+1, N+1)   :: omega
    !---------------------------------------------------------------------!
    ! Domain variables                                                    !
    !---------------------------------------------------------------------!
    Re      = 100.0_WP              ! Reynolds number
    L       = 1.0_WP                ! length and height of cavity
    nu      = 0.01_WP               ! kinematic viscosity
    Uwall   = Re*nu/L               ! velocity at top wall
    !---------------------------------------------------------------------!
    ! Setting the wall velocities                                         !
    !---------------------------------------------------------------------!
    Urx     = 0.0_WP                ! right wall u-velocity               
    Ulx     = 0.0_WP                ! left wall u-velocity
    Utx     = 1.0_WP                ! top wall u-velocity
    Ubx     = -1.0_WP               ! bottom wall u-velocity
    Ury     = -1.0_WP               ! right wall v-velocity
    Uly     = 1.0_WP                ! left wall v-velocity
    Uty     = 0.0_WP                ! top wall v-velocity
    Uby     = 0.0_WP                ! bottom wall v-velocity
    !---------------------------------------------------------------------!
    ! Steady state convergence criteria                                   !
    !---------------------------------------------------------------------!
    conv_crit = (10.0_WP)**(-7.0_WP)  ! convergence criteria
    !---------------------------------------------------------------------!
    ! Discretization domain                                               !
    !---------------------------------------------------------------------!
    dx  = L/M                   ! x-step size
    dy  = L/N                   ! y-step size
    rh  = 1.0_WP/dx             ! reciprocal of h (dx = dy)
    rh2 = 1.0_WP/dx**2          ! reciprocal of h^2
    rRe = 1.0_WP/Re             ! reciprocal of Re
    !---------------------------------------------------------------------!
    ! Preallocating fields                                                !
    !---------------------------------------------------------------------!
    u       = 0.0_WP            ! u-velocity field
    v       = 0.0_WP            ! v-velocity field
    u_star  = 0.0_WP            ! u-star-velocity field (no pressure)
    v_star  = 0.0_WP            ! v-star-velocity field (no pressure)
    phi     = 0.0_WP            ! Lagrangian multiplier
    res_gs  = 0.0_WP            ! Guass Seidel residual
    gs_RHS  = 0.0_WP            ! right hand side of GS
    gs_RHS2 = 0.0_WP            ! right hand side of GS
    omega   = 0.0_WP            ! vorticity
    !---------------------------------------------------------------------!
    ! Time counter variables                                              !
    !---------------------------------------------------------------------!
    tfinal  = 30.0_WP           ! final time step
    n_count = 0                 ! counter
    conv    = 5                 ! initial convergence value
    n_gs_t  = 0                 ! gauss-seidel counter (total number of iterations for all time)
    t       = 0.0_WP            ! time
    !---------------------------------------------------------------------!
    ! Boundary conditions                                                 !
    !---------------------------------------------------------------------!
    u(:,N+2)    = 2.0_WP*Utx-u_star(:,N+1)      ! u-star top wall    
    u(:,1)      = 2.0*Ubx-u_star(:,2)           ! u-star bottom wall 
    u(1,:)      = Ulx                           ! u-star left wall   
    u(M+1,:)    = Urx                           ! u-star right wall  
    v(:,N+1)    = Uty                           ! v-star top wall    
    v(:,1)      = Uby                           ! v-star bottom wall 
    v(1,:)      = 2.0_WP*Uly-v_star(2,:)        ! v-star left wall   
    v(M+2,:)    = 2.0_WP*Ury-v_star(M+1,:)      ! v-star right wall  
    !=====================================================================!
    ! Time loop                                                           !
    !=====================================================================!
    do while (t < tfinal)
        !-----------------------------------------------------------------!
        ! time step calculation                                           !
        !-----------------------------------------------------------------!
        sf      = 0.9_WP                        ! time step safety factor for stability
        a_grid  = u                             ! a value for the grid
        b_grid  = v                             ! b value for the grid
        a       = maxval(abs(a_grid))           ! maximum a value
        b       = maxval(abs(b_grid))           ! maximum b value
        dt1     = 0.25_WP*dx**(2.0_WP)/nu       ! parabolic time step constraint
        dt2     = dx/(a+b)                      ! hyperbolic time step constraint
        dt      = sf*min(dt1, dt2)              ! calculating time step
        !-----------------------------------------------------------------!
        ! Checking cell Reynolds number stability                         !
        !-----------------------------------------------------------------!
        if ((a+b)*dx/nu > 2) then
            print *, 'Unstable - Cell Reynolds Number'
        end if
        !-----------------------------------------------------------------!
        ! Increasing counters and resetting counters                      !
        !-----------------------------------------------------------------!
        t           = t+dt      ! increase time by stable dt
        n_count     = n_count+1 ! counter
        n_gs        = 0         ! gauss-seidel counter for each set of loops
        !=================================================================!
        ! Step I of Fractional Step Method                                !
        !=================================================================!
        !-----------------------------------------------------------------!
        ! Storing solutions to time step t^{n}                            !
        !-----------------------------------------------------------------!
        r_u = u                ! u-vel. solution for t^{n} for convergence
        r_v = v                ! v-vel. solution for t^{n} for convergence
        !-----------------------------------------------------------------!
        ! Calculating u-star values                                       !
        !-----------------------------------------------------------------!
        do j = 2, N+1
            do i = 2, M
                !---------------------------------------------------------!
                ! Calculating u-vel. and v-vel. average values            !
                !---------------------------------------------------------!
                ua = 0.5_WP*(u(i+1,j)+u(i,j))          ! u_i+1,j
                ub = 0.5_WP*(u(i,j)+u(i-1,j))          ! u_i,j
                uc = 0.5_WP*(u(i,j+1)+u(i,j))          ! u_i+1/2,j+1/2
                ud = 0.5_WP*(u(i,j-1)+u(i,j))          ! u_i+1/2,j-1/2
                va = 0.5_WP*(v(i+1,j)+v(i,j))          ! v_i+1/2,j+1/2
                vb = 0.5_WP*(v(i+1,j-1)+v(i,j-1))      ! v_i+1/2,j-1/2
                !---------------------------------------------------------!
                ! u-star values                                           !
                !---------------------------------------------------------!
                u_star(i,j) = u(i,j)+ dt*(&
                                -((ua**2.0_WP-ub**2.0_WP)*rh &
                                +(uc*va-ud*vb)*rh) &
                                + rRe*((u(i+1,j)-2*u(i,j) + u(i-1,j))*rh2 &
                                +(u(i,j+1)-2*u(i,j) + u(i,j-1))*rh2))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Calculating v-star values                                       !
        !-----------------------------------------------------------------!
        do j = 2, N
            do i = 2, M+1
                !---------------------------------------------------------!
                ! Calculating u-vel. and v-vel. average values            !
                !---------------------------------------------------------!
                uc = 0.5_WP*(u(i,j+1)+u(i,j))      !u_i+1/2,j+1/2
                ue = 0.5_WP*(u(i-1,j)+u(i-1,j+1))  !u_i-1/2,j+1/2
                va = 0.5_WP*(v(i+1,j)+v(i,j))      !v_i+1/2,j+1/2
                vc = 0.5_WP*(v(i,j+1)+v(i,j))      !v_i,j+1
                vd = 0.5_WP*(v(i,j)+v(i,j-1))      !v_i,j
                vf = 0.5_WP*(v(i,j)+v(i-1,j))      !v_i-1/2,j+1/2
                !---------------------------------------------------------!
                ! Calculating v-star values                               !
                !---------------------------------------------------------!
                v_star(i,j) = v(i,j)+dt*(&
                                -((uc*va-ue*vf)*rh &
                                +(vc**2-vd**2)*rh) &
                                +rRe*((v(i+1,j)-2*v(i,j)+v(i-1,j))*rh2&
                                +(v(i,j+1)-2*v(i,j)+ v(i,j-1))*rh2))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Boundary conditions                                             !
        !-----------------------------------------------------------------!
        u_star(:,N+2)   = 2.0_WP*Utx-u_star(:,N+1)      ! u-star top wall
        u_star(:,1)     = 2.0_WP*Ubx-u_star(:,2)        ! u-star bottom wall
        u_star(1,:)     = Ulx                           ! u-star left wall
        u_star(M+1,:)   = Urx                           ! u-star right wall
        v_star(:,N+1)   = Uty                           ! v-star top wall
        v_star(:,1)     = Uby                           ! v-star bottom wall
        v_star(1,:)     = 2.0_WP*Uly-v_star(2,:)        ! v-star left wall
        v_star(M+2,:)   = 2.0*Ury-v_star(M+1,:)         ! v-star right wall
        !=================================================================!
        ! Step II of Fractional Step Method                               !
        !=================================================================!
        !-----------------------------------------------------------------!
        ! Dynamic GS convergence conditions                               !
        !-----------------------------------------------------------------!
        conv_gs         = 1                     ! resetting the GS convergence
        conv_gs_limit   = (10.0_WP)**(-9.0_WP)
        !print *, conv
        !if  (conv*(10.0_WP)**(-2.0_WP) >= (10.0_WP)**(-3.0_WP)) then
        !    conv_gs_limit = 10.0_WP**(-12.0_WP)      ! convergence criteria for early in time
        !elseif (conv*10.0_WP**(-2.0_WP) <= 10.0_WP**(-7.0_WP)) then
        !    conv_gs_limit = 10.0_WP**(-12.0_WP)     ! convergence criteria near steady state
        !end if
        !-----------------------------------------------------------------!
        ! Finding solution to Lagrange multiplier using Gauss-Seidel      !
        ! iterates until Gauss-Seidel convergence criteria is met         !
        !-----------------------------------------------------------------!
        do j = 2, N+1       ! pre-calculate RHS of GS and residual equations
            do i = 2, M+1
                gs_RHS(i,j) = 0.25_WP*(dx/dt)*(u_star(i,j)-u_star(i-1,j)+&
                                v_star(i,j)-v_star(i,j-1))
                gs_RHS2(i,j) = (1.0_WP/dt)*(((u_star(i,j)-u_star(i-1,j))*rh)+&
                                ((v_star(i,j)-v_star(i,j-1))*rh))
            end do
        end do
        !-----------------------------------------------------------------!
        ! GS iterations                                                   !
        !-----------------------------------------------------------------!
        do while (conv_gs > conv_gs_limit .and. n_gs < 100000)
            n_gs    = n_gs+1                ! GS iteration counter
            n_gs_t  = n_gs_t+1              ! GS iteration counter (total)
            do j = 2, N+1
                do i = 2,M+1
                    phi(i,j) = 0.25_WP*(phi(i-1,j)+phi(i+1,j)+phi(i,j-1) &
                        + phi(i,j+1))-gs_RHS(i,j)
                end do
            end do
            !-------------------------------------------------------------!
            ! Update Lagrange multiplier ("pressure") BCs                 !
            !-------------------------------------------------------------!
            phi(:,1)    = phi(:,2)             ! phi bottom wall
            phi(:,N+2)  = phi(:,N+1)           ! phi top wall
            phi(1,:)    = phi(2,:)             ! phi left wall
            phi(M+2,:)  = phi(M+1,:)           ! phi right wall
            !-------------------------------------------------------------!
            ! Check convergence of Lagrange multiplier                    !
            !-------------------------------------------------------------!
            do j = 2, N+1
                do i = 2, M+1
                    res_gs(i,j) = (phi(i+1,j)-2.0_WP*phi(i,j)+phi(i-1,j))*rh2 &
                                    +(phi(i,j+1)-2.0_WP*phi(i,j)+phi(i,j-1))*rh2 &
                                    -gs_RHS2(i,j)   ! residual for entire mesh
                end do
            end do
            conv_gs = maxval(abs(res_gs))           ! check convergence
            if (n_gs > 100000) then
                print *, 'GS did not converge after 100k iterations'
            end if 
        end do
        !=================================================================!
        ! Step III of Fractional Step Method                              !
        !=================================================================!
        !-----------------------------------------------------------------!
        ! Calculating u-velocity solution at n+1                          !
        !-----------------------------------------------------------------!
        do j = 2, N+1
            do i = 2, M
                u(i,j) = u_star(i,j)-(dt*rh)*(phi(i+1,j)-phi(i,j))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Calculating v-velocity solution at n+1                          !
        !-----------------------------------------------------------------!
        do j = 2,N
            do i = 2,M+1
                v(i,j) = v_star(i,j)-(dt*rh)*(phi(i,j+1)-phi(i,j))
            end do
        end do
        !-----------------------------------------------------------------!
        ! Update velocities BCs                                           !
        !-----------------------------------------------------------------!
        u(:,N+2)    = 2.0_WP*Utx-u(:,N+1)           ! u top wall    
        u(:,1)      = 2.0_WP*Ubx-u(:,2)             ! u bottom wall 
        u(1,:)      = Ulx                           ! u left wall   
        u(M+1,:)    = Urx                           ! u right wall  
        v(:,N+1)    = Uty                           ! v top wall    
        v(:,1)      = Uby                           ! v bottom wall 
        v(1,:)      = 2.0_WP*Uly-v(2,:)             ! v left wall   
        v(M+2,:)    = 2.0_WP*Ury-v(M+1,:)           ! v right wall  
        !-----------------------------------------------------------------!
        ! Checking the velocity convergence                               !
        !-----------------------------------------------------------------!
        res_u   = (u-r_u)/dt                    ! rate change of u-velocity
        res_v   = (v-r_v)/dt                    ! rate change of v-velocity
        conv_u  = maxval(abs(res_u))            ! u-vel. convergence
        conv_v  = maxval(abs(res_v))            ! v-vel. convergence
        conv    = max(conv_u, conv_v)           ! maximum vel. rate of change
        !-----------------------------------------------------------------!
        ! Storing convergence history                                     !
        !-----------------------------------------------------------------!
        conv_hist_u(n) = conv_u                ! u-velocity convergence
        conv_hist_v(n) = conv_v                ! v-velocity convergence
        !-----------------------------------------------------------------!
        ! Time step print statement                                       !
        !-----------------------------------------------------------------!
        print "(4X, A, ES10.3, A, ES10.3)", &
                    'Convergence of u = ', conv_u,&
                    ' Convergence of v = ', conv_v
        print "(4X,A, ES10.3, /, 4X, A, I8, /, 4X,  A, ES10.3, /, 4X, A, ES10.3)",&
                'Simulation time ---> ', t, &
                'Iteration --->', n_count, 'GS Iterations -->', real(n_gs), &
                'Total GS Iterations --->', real(n_gs_t)
    end do
    !---------------------------------------------------------------------!
    ! Calculating vorticity                                               !
    !---------------------------------------------------------------------!
    do j = 1, N+1
        do i = 1, M+1
            omega(i,j) = rh*(v(i+1,j)-v(i,j)-u(i,j+1)+u(i,j));
        end do
    end do
    !---------------------------------------------------------------------!
    ! writing vorticity                                                   !
    !---------------------------------------------------------------------!
    open(unit=1, file='../data/vorticity.dat')
    do i = 1, M+1
        write(1,*) ( omega(i,j), j=1,N+1 )
    end do
    close(unit=1)
    !---------------------------------------------------------------------!
    ! writing u-velocity                                                  !
    !---------------------------------------------------------------------!
    open(unit=1, file='../data/u-velocity.dat')
    do i = 1, M+1
        write(1,*) ( u(i,j), j=1,N+2 )
    end do
    close(unit=1)
    !---------------------------------------------------------------------!
    ! writing u-velocity                                                  !
    !---------------------------------------------------------------------!
    open(unit=1, file='../data/v-velocity.dat')
    do i = 1, M+2
        write(1,*) ( v(i,j), j=1,N+1 )
    end do
    close(unit=1)

end program ns_solve
