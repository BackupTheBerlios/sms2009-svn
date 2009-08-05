!> 1D Convection Diffusion equations solver in Fortran
!!
!! Solves the equation:
!!
!!\f[
!!\frac{du}{dt}=\frac{d}{dx}D(x)\frac{du}{dx} + C(x)\frac{du}{dx}+F(x)u-S(x)
!!\f]
!! for \f$u\f$, given functions for \f$D\f$, \f$C\f$, \f$F\f$, and \f$S\f$, defined in this program
!!
!! Implicit methods are used
!!
!! \author Jesse V. Johnson (jvj)
!! \date 6-27-09

program implicit_convdiff

  use derivatives
  use solvers

  implicit none

  ! local variables

  integer, parameter :: N=100                   ! Number of nodes
  real,    parameter :: dt = 0.005              ! length time step 
  integer, parameter :: nt = 100                ! number of time steps
  integer :: t                                  ! current time step
  real :: xl                                    ! start of domain
  real :: xr                                    ! end of domain 
  real :: dx                                    ! node spacing

  real, dimension(:), allocatable :: x          ! node positions
  real, dimension(:), allocatable :: u          ! what we're solving for
  real, dimension(:), allocatable :: rhs        ! Matrix equation right-hand side
  real, dimension(:), allocatable :: unew       ! New values of u
  real, dimension(:), allocatable :: uold       ! Old values of u
  real, dimension(:), allocatable :: dddx       ! Spatial derivative of diffusivity
  real, dimension(:), allocatable :: sub        ! Subdiagonal terms
  real, dimension(:), allocatable :: dia        ! Diagonal terms
  real, dimension(:), allocatable :: sup        ! Superdiagonal terms

  real, dimension(:), allocatable :: d          ! diffusivity coeff
  real, dimension(:), allocatable :: c          ! convection velocity
  real, dimension(:), allocatable :: s          ! source term

  integer :: ii                                 ! a counter
  integer :: errstat                            ! for error checking
  
  ! let's allocate some memory
  allocate(x(N),stat=errstat)
  call checkerr(errstat,"failed to allocate x")

  allocate(u(N),stat=errstat)
  call checkerr(errstat,"failed to allocate u")

  allocate(rhs(N),stat=errstat)
  call checkerr(errstat,"failed to allocate rhs")

  allocate(unew(N),stat=errstat)
  call checkerr(errstat,"failed to allocate unew")

  allocate(uold(N),stat=errstat)
  call checkerr(errstat,"failed to allocate uold")

  allocate(dddx(N),stat=errstat)
  call checkerr(errstat,"failed to allocate dddx")

  allocate(sub(N),stat=errstat)
  call checkerr(errstat,"failed to allocate sub")

  allocate(dia(N),stat=errstat)
  call checkerr(errstat,"failed to allocate dia")

  allocate(sup(N),stat=errstat)
  call checkerr(errstat,"failed to allocate sup")

  allocate(d(N),stat=errstat)
  call checkerr(errstat,"failed to allocate d")

  allocate(c(N),stat=errstat)
  call checkerr(errstat,"failed to allocate c")

  allocate(s(N),stat=errstat)
  call checkerr(errstat,"failed to allocate s")

  ! Set up grid 
  ! Space
  xl = 0.0
  xr = 1.0
  dx = real((xr - xl) / (N-1))

  ! Initial value
  do ii=1,N
     x(ii)  =  xl + real(dx * (ii-1)) * (xr - xl)
  enddo
  
  ! Initial value
  ! For now, start with a square wave.
  uold = 0.
  uold(N/2-int(N*.05):N/2+int(N*.05)) = 1.

  ! Set coefficients via function calls
  ! Thankfully, these are not time dependent.
  d = dfun(x)
  c = cfun(x)
  s = sfun(x)

  ! Need the derivative of d for calculations
  dddx = dudx(d,dx)

  ! Zero out the tridiagonal matrix
  sub = 0.
  dia = 0.
  sup = 0.

  ! Loop to place the coefficients into the matrix.
  ! Avoid boundaries this pass.
  do ii = 2,N-1
     ! Left of point in question, or subdiagonal
     sub(ii) = dt / dx * (dddx(ii) / 2. + d(ii) / dx + c(ii) / 2.)
     ! Point in question, or diagonal
     dia(ii) = dt  * (-2. * d(ii) / dx**2. - 1. / dt)
     ! Right of point in questions, superdiagonal
     sup(ii) = dt / dx * (dddx(ii) / 2. + d(ii) / dx - c(ii) / 2.)
     ! Right hand side, or b in A x = b
     rhs(ii) = s(ii) * dt 
  end do

  ! Set Dirchlet boundaries
  dia(1) = 1.
  dia(N) = 1.
  rhs(1) = 0.
  rhs(N) = 0.

  time_loop: do t =1,nt
     call tridiag(sub,dia,sup,unew,rhs - uold)
     write (*,*) unew
     uold = unew
  end do time_loop

  ! check that arrays are allocated before attempting to clean up
  if(allocated(x))    deallocate(x)
  if(allocated(u))    deallocate(u)
  if(allocated(rhs))  deallocate(rhs)
  if(allocated(unew)) deallocate(unew)
  if(allocated(uold)) deallocate(uold)
  if(allocated(dddx)) deallocate(dddx)
  if(allocated(sub))  deallocate(sub)
  if(allocated(dia))  deallocate(dia)
  if(allocated(sup))  deallocate(sup)

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
    cfun = -1.0
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

end program implicit_convdiff

