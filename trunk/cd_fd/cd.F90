!> 1D Convection Diffusion equations solver in Fortran
!!
!! Solves the equation:
!!
!!\f[
!!\frac{du}{dt}=\frac{d}{dx}D(x)\frac{du}{dx} + C(x)\frac{du}{dx}+F(x)u-S(x)
!!\f]
!! for \f$u\f$, given functions for \f$D\f$, \f$C\f$, \f$F\f$, and \f$S\f$, defined in this program
!!
!! Explicit methods are used
!!
!! \author Jesse V. Johnson (jvj)
!! \date 6-9-09

program explicit_convdiff
  use derivatives
  implicit none

  integer, parameter :: N=50 ! Number of nodes
  real :: dt, t, tfinal ! Time domain
  real :: xl, xr, dx   ! Space
  real, dimension(:), allocatable :: x,u,dudx_store,convection,diffusion
  real, dimension(:), allocatable :: d,c,s,f
  integer i

  allocate(u(N))
  allocate(x(N))
  allocate(d(N))
  allocate(c(N))
  allocate(s(N))
  allocate(f(N))
  allocate(dudx_store(N))
  allocate(convection(N))
  allocate(diffusion(N))

  ! Set up grid 
  ! Space
  xl = 0.
  xr = 1.
  dx = real((xr - xl) / (N-1))
  ! Time
  t = 0.
  tfinal = 0.05
  dt = .0001

  do i=1,N
     x(i) = xl + real(dx * (i-1)) * (xr - xl)
  enddo

  ! Initial value
  u = 0.

  ! For now, start with a square wave.
  u(N/2-3:N/2+3) = 1.

  ! Thankfully, these are not time dependent.
  d = dfun(x)
  c = cfun(x)
  f = ffun(x)
  s = sfun(x)

  time_loop: do while (t <= tfinal)

     dudx_store = dudx(u,dx)

     convection = c * dudx_store 
     diffusion  = dudx(d,dx) * dudx_store + d * d2udx2(u,dx)

     u = u + dt * (convection + diffusion + f * u + s)

     t = t + dt
     u(1)=0.
     u(N)=0.

     write (*,*) u
  end do time_loop

  deallocate(u)
  deallocate(x)
  deallocate(d)
  deallocate(c)
  deallocate(f)
  deallocate(s)

contains

  elemental function dfun(x)
    implicit none
    real, intent(in) :: x
    real dfun
    dfun = 0.
  end function dfun

  elemental function cfun(x)
    implicit none
    real, intent(in) :: x
    real cfun
    cfun = 100.
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

end program explicit_convdiff

