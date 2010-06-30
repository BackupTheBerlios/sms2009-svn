! 1D Diffusion equations solver in Fortran
! Solves the equation:
!
! du     d^2u  
! -- = D ---- 
! dt     dx^2  
! 
! for u, D
!
! Explicit methods are used
!
! Jesse V. Johnson (jvj)
! 6-9-09
! adapted by Ian Rutt
! 7-23-09

program diffusion

  use derivatives

  implicit none

  integer, parameter :: N=50 ! Number of nodes
  real :: dt, t, tfinal ! Time domain
  real :: xl, xr, dx   ! Space
  real, dimension(:), allocatable :: x,u
  real :: d  ! Diffusivity
  real :: s  ! Source
  integer i

  allocate(u(N))
  allocate(x(N))

  ! Set up grid 
  ! Space
  xl = 0.
  xr = 1.
  dx = real((xr - xl) / (N-1))
  ! Time
  t = 0.
  tfinal = 0.5
  dt = .00015

  do i=1,N
     x(i) = xl + real(dx * (i-1)) * (xr - xl)
  enddo

  ! Initial value
  u = 0.

  ! For now, start with a square wave.
  u(N/2-3:N/2+3) = 1.

  ! Thankfully, these are not time dependent.
  d = 1.0
  s = 0.0005
  
  time_loop: do while (t <= tfinal)

     u = u + dt * d * d2udx2(u,dx) + s

     t = t + dt
     u(1)=0.
     u(N)=0.

     write (*,*) u
  end do time_loop

  deallocate(u)
  deallocate(x)

end program diffusion

