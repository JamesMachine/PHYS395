! compile with: gfortran -o3 -fdefault-real-8 goldenmin.f90 -o goldenmin
program goldenmin
implicit none

integer :: i, j, max_iter = 1000
integer, parameter :: m = 20
real :: a, b, c, d
real, parameter :: frac = 0.61803, eps = 10e-10
real :: x(m), min = -3.0, max = 3.0



forall (j=1:m) x(j) = min + (max - min)*(j-1)/(m-1)

a = -2.0; b = -1.5; c = 0.0

write(*,*) "Initial bracketing values:", a, b, c
write(*,*) "Will shift bracket by m =", m, "iterations between", min, ",", max
write(*,*) ""
write(*,*) ""

! loop over m different sets of initial brackets
do i = 1,m

	a = -2.0 + x(i); b = -1.5 + x(i); c = 0.0  + x(i)

	!write(*,*) "Initial values:", a, b, c

	if (.not.(f(b) < f(a) .and. f(b) < f(c))) then
		!write(*,*) "ERROR: minimum is not bracketed"
		!write(*,*) ""
		cycle
	end if

	! interval bracketing loop to find minimum
	j = 1
	do
		if (abs(a-b) > abs(b-c)) then 				! [ab] interval
			d = a + (b - a)*frac
			if (f(d) < f(a) .and. f(d) < f(b)) then
				c = b
				b = d
			else if (f(b) < f(d) .and. f(b) < f(c)) then
				a = d
			end if
		else 										! [bc] interval
			d = b + (c - b)*frac
			if (f(b) < f(a) .and. f(b) < f(d)) then
				c = d
			else if (f(d) < f(b) .and. f(d) < f(c)) then
				a = b
				b = d
			end if
		end if

		if (abs(df(d)) < eps) then
			write(*,*) "Successfully converged to root:"
			write(*,*) "iteration #,     x,     f(x),     df(x)/dx"
			write(*,*) j, d, f(d), df(d)
			exit
		else
			j = j + 1
			if (j == max_iter) then
				write(*,*) "Reached maximum number of iterations, &
				restarting with next set of brackets"
				write(*,*) "iteration #,     x,     f(x),     df(x)/dx"
				write(*,*) j, d, f(d), df(d)
				write(*,*) ""
				j = 1
				exit
			end if
		end if
	end do

end do



contains

pure function f(x)
	intent(in) x
	real :: f, x
	f = (x**2 - 1)**2 + x
end function

pure function df(x)
	intent(in) x
	real :: df, x
	df = 2.0*(x**2 - 1)*2.0*x + 1.0
end function

end program