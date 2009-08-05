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

  ! local variables
  real, parameter    :: pi = 3.1415926536

  integer, parameter :: N=51                    ! Number of nodes
  real,    parameter :: dt = 0.000125           ! length time step 
  integer, parameter :: nt = 8208               ! number of time steps
  integer :: t                                  ! current time step
  real :: xl                                    ! start of domain
  real :: xr                                    ! end of domain 
  real :: dx                                    ! node spacing

  real, dimension(:), allocatable :: x          ! node positions
  real, dimension(:), allocatable :: u          ! what we're solving for
  real, dimension(:), allocatable :: dudx_store ! space to store du/dx
  real, dimension(:), allocatable :: convection ! convection term
  real, dimension(:), allocatable :: diffusion  ! diffusion term
  real, dimension(:), allocatable :: ua         ! analytical sol'n

  real, dimension(:), allocatable :: d          ! diffusivity coeff
  real, dimension(:), allocatable :: c          ! convection velocity
  real, dimension(:), allocatable :: s          ! source term

  integer :: ii                                 ! a counter
  integer :: errstat                            ! for error checking
  
  ! let's allocate some memory
  allocate(u(N),stat=errstat)
  call checkerr(errstat,"failed to allocate u")

  allocate(ua(N),stat=errstat)
  call checkerr(errstat,"failed to allocate ua")

  allocate(x(N),stat=errstat)
  call checkerr(errstat,"failed to allocate x")

  allocate(d(N),stat=errstat)
  call checkerr(errstat,"failed to allocate d")

  allocate(c(N),stat=errstat)
  call checkerr(errstat,"failed to allocate c")

  allocate(s(N),stat=errstat)
  call checkerr(errstat,"failed to allocate s")

  allocate(dudx_store(N),stat=errstat)
  call checkerr(errstat,"failed to allocate dudx_store")
 
  allocate(convection(N),stat=errstat)
  call checkerr(errstat,"failed to allocate convection")
  
  allocate(diffusion(N),stat=errstat)
  call checkerr(errstat,"failed to allocate diffusion")

  ! Set up grid 
  ! Space
  xl = 0.0
  xr = 1.0
  dx = real((xr - xl) / (N-1))

  ! Initial value
  do ii=1,N
     x(ii)  =  xl + real(dx * (ii-1)) * (xr - xl)
     u(ii)  =  sin(pi * x(ii))
     ua(ii) =  sin(pi * x(ii))
  enddo
  
  ! Set coefficients via function calls
  ! Thankfully, these are not time dependent.
  d = dfun(x)
  c = cfun(x)
  s = sfun(x)

  time_loop: do t=1,nt 

     dudx_store = dudx(u,dx)

     convection = c * dudx_store 
     diffusion  = dudx(d,dx) * dudx_store + d * d2udx2(u,dx)

     u = u + dt * (convection + diffusion + s)
     ua = sin(pi*x) * exp(-pi**2*(t*dt))

     ! boundary conditions enforced
     u(1)=0.0
     u(N)=0.0

     ! Write out normalised error w.r.t. analytical sol'n
     write (*,*) 100.d0 * abs(u(2:N-1) - ua(2:N-1)) / ua(2-N-1)

  end do time_loop

  ! check that arrays are allocated before attempting to clean up
  if(allocated(u))  deallocate(u)
  if(allocated(ua)) deallocate(ua)
  if(allocated(x))  deallocate(x)
  if(allocated(d))  deallocate(d)
  if(allocated(c))  deallocate(c)
  if(allocated(s))  deallocate(s)

contains

  elemental function dfun(x)
    implicit none
    real, intent(in) :: x
    real dfun
    dfun = 1.
  end function dfun

  elemental function cfun(x)
    implicit none
    real, intent(in) :: x
    real cfun
    cfun = 0.
  end function cfun

  elemental function sfun(x)
    implicit none
    real, intent(in) :: x
    real sfun
    sfun = 0.
  end function sfun

  subroutine checkerr(errstat,msg)
    implicit none
    integer,      intent(in) :: errstat
    character(*), intent(in) :: msg 
    if (errstat /= 0) then
       write(*,*) "ERROR:", msg
       stop
    end if
  end subroutine checkerr

end program explicit_convdiff

