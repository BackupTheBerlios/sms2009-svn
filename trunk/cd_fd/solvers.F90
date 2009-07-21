module solvers
    contains
      subroutine tridiag(a,b,c,x,y)

        !*FD Tridiagonal solver. All input/output arrays should have the 
        !*FD same number of elements.

        real,dimension(:),intent(in) ::  a !*FD Lower diagonal; a(1) is ignored.
        real,dimension(:),intent(in) ::  b !*FD Centre diagonal
        real,dimension(:),intent(in) ::  c !*FD Upper diagonal; c(n) is ignored.
        real,dimension(:),intent(out) :: x !*FD Unknown vector
        real,dimension(:),intent(in) ::  y !*FD Right-hand side

        real,dimension(size(a)) :: aa
        real,dimension(size(a)) :: bb

        integer :: n,i

        n=size(a)

        aa(1) = c(1)/b(1)
        bb(1) = y(1)/b(1)

        do i=2,n
           aa(i) = c(i)/(b(i)-a(i)*aa(i-1))
           bb(i) = (y(i)-a(i)*bb(i-1))/(b(i)-a(i)*aa(i-1))
        end do
        
        x(n) = bb(n)

        do i=n-1,1,-1
           x(i) = bb(i)-aa(i)*x(i+1)
        end do

  end subroutine tridiag

end module solvers
