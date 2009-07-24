!> compute derivatives
!! \author Jesse V. Johnson (jvj)

module derivatives

contains
  
  !> compute first derivative using central differences
  function dudx(u,dx)
    implicit none
    real, dimension(:), intent(in) :: u !< the quantity for which the derivative should be computed
    real, intent(in) :: dx              !< the grid spacing
    real, dimension(size(u)) :: dudx
    integer :: N 
    N = size(u)

    dudx(2:N-1) = (u(3:) -  u(1:N-2)) / (2 * dx)
    dudx(1) =  (- 3. * u(1) + 4. * u(2)     - u(3))     / (2. * dx)
    dudx(N) =  (  3. * u(N) - 4. * u(N - 1) + u(N - 2)) / (2. * dx)

  end function dudx

  !> compute the second derivative
  function d2udx2(u,dx)
    implicit none
    real, intent(in), dimension(:) :: u !< the quantity for which the derivative should be computed
    real, intent(in) :: dx              !< the grid spacing
    real, dimension(size(u)) :: d2udx2
    integer :: N
    N = size(u)

    d2udx2(2:N-1) = (u(3:) - 2 * u(2:N-1) + u(1:N-2)) / dx**2
    d2udx2(1) = ( 2. * u(1) - 5. * u(2) + 4. * u(3) - u(4)) / dx**2
    d2udx2(N) = ( 2. * u(N) - 5. * u(N-1) + 4. * u(N-2) - u(N-3)) / dx**2
  end function d2udx2

end module derivatives
