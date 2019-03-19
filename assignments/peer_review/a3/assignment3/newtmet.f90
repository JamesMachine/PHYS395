! compile with: gfortran -o3 -fdefault-real-8 newtmet.f90 -o newtmet
program newtmet
implicit none

integer :: i, j, n = 8
integer, parameter :: m = 20
real :: x(m), min = -3.0, max = 3.0

forall (i=1:m) x(i) = min + (max - min)*(i-1)/(m-1)

! scans initial values in order to find all 3 roots,
! for each initial value takes n Gauss-Newton steps
do i = 1,m
	write(*,*) "Initial value:", x(i)
	do j = 1,n
		x(i) = x(i) - f(x(i))/df(x(i))
		write(*,*) x(i), f(x(i))
	end do
end do



contains

pure function f(x)
	intent(in) x
	real :: f, x
	f = x**3 - x + 1.0/4.0
end function

pure function df(x)
	intent(in) x
	real :: df, x
	df = 3.0*x**2 - 1.0
end function

end program