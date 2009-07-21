! 1D Convection Diffusion equations solver in Fortran
! Solves the equation:
!
! du   d       du         d u
! -- = -- D(x) --- + C(x) ---  + F(x) u - S(x)
! dt   dx      dx         d x
! 
! for u, given functions for D, C, F, and S, defined in this program
!
! Implicit methods are used
!
! Jesse V. Johnson (jvj)
! 6-27-09

program convdiff
    use derivatives
    use solvers

    implicit none
    integer, parameter :: N=100 ! Number of nodes

    real :: dt, t, tfinal ! Time domain
    real :: xl, xr, dx    ! Space domain

    real, dimension(N) :: x,u,rhs,unew,uold
    real, dimension(N) :: d,c,s,f,dddx,sub,dia,sup

    integer i

   ! Set up grid 
    ! Space
    xl = 0.
    xr = 1.
    dx = real((xr - xl) / (N-1))

    do i=1,N
        x(i) = xl + real(dx * (i-1)) * (xr - xl)
    enddo

    ! Time
    t = 0.
    tfinal = 0.5

    ! Initial value
    ! For now, start with a square wave.
    uold = 0.
    uold(N/2-int(N*.05):N/2+int(N*.05)) = 1.

    ! Thankfully, these are not time dependent, so we can simiply fill the
    ! arrays
    d = dfun(x)
    c = cfun(x)
    f = ffun(x)
    s = sfun(x)

    ! Need the derivative of d for calculations
    dddx = dudx(d,dx)

    ! Compute the time step based on CFL condition
    ! dt =min(dx**2./max(d),dx/max(c))
    dt =.005

    ! Zero out the tridiagonal matrix
    sub = 0.
    dia = 0.
    sup = 0.

    ! Loop to place the coefficients into the matrix.
    ! Avoid boundaries this pass.
     do i = 2,N-1
            ! Left of point in question, or subdiagonal
            sub(i) = dt / dx * (dddx(i) / 2. + d(i) / dx + c(i) / 2.)
            ! Point in question, or diagonal
            dia(i) = dt  * (-2. * d(i) / dx**2. + f(i) - 1. / dt)
            ! Right of point in questions, superdiagonal
            sup(i) = dt / dx * (dddx(i) / 2. + d(i) / dx - c(i) / 2.)
            ! Right hand side, or b in A x = b
            rhs(i) = s(i) * dt 
     end do 

    ! Set Dirchlet boundaries
     dia(1) = 1.
     dia(N) = 1.
     rhs(1) = 0.
     rhs(N) = 0.

    time_loop: do while (t <= tfinal)
        call tridiag(sub,dia,sup,unew,rhs - uold)
        write (*,*) unew
        uold = unew
        t = t + dt
    end do time_loop

   contains

    elemental function dfun(x)
     implicit none
     real, intent(in) :: x
     real dfun
     dfun = 0.0
    end function dfun

    elemental function cfun(x)
     implicit none
     real, intent(in) :: x
     real cfun
     cfun = -1.
    end function cfun

    elemental function ffun(x)
     implicit none
     real, intent(in) :: x
     real ffun
     ffun = 0.
    end function ffun

    elemental function sfun(x)
     implicit none
     real, intent(in) :: x
     real sfun
     sfun = 0.
    end function sfun
       
end program convdiff

