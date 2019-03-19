! compile with: gfortran -o3 -fdefault-real-8 levmar.f90 -o levmar
program levmar
implicit none

real, parameter :: pi = 4.0*atan(1.0)
real, allocatable :: A(:,:), Areg(:,:), B(:), c(:)
real, allocatable :: xx(:), yy(:), sigma(:)
real :: csq, dcsq, x, y, lambda, dt
integer :: i, j
integer :: n, m, status = 0, arg_status = 0, iterations = 100
character(len=10) :: arg
character(len=20) :: data_filename = "a2.dat"



! n = n + 2 because the 0th basis function and the constant must be included
call get_command_argument(1,arg)
read(arg,*,iostat = arg_status) n
n = n + 2

allocate(A(n,n),Areg(n,n),B(n),c(n))

! read data file to determine its size
open(unit=13,file=data_filename,status="old",action="read")
	m = 0
	do while (status == 0)
		read(13,*,iostat = status) x, y
		if (status < 0) exit
		m = m + 1
	end do
close(13)

allocate(xx(m),yy(m),sigma(m))

! read data again from stdin to accumulate relevant vectors
status = 0
i = 1
do while (status == 0)
	read(*,*,iostat = status) x, y
	xx(i) = x; yy(i) = y
	if (status < 0) exit

	i = i + 1
end do



! begin "Marquardt recipe"

! starting values (sigma_i^2 = 1 for this problem,
! c is function parameter vector)
sigma = 1.0
c = [1.5,0.6,0.2,0.2,-6.8]
lambda = 0.001

! evaluate initial chi-squared
csq = 0.0
do i = 1,m
	csq = csq + chisq(xx(i),yy(i),sigma(i),c,n)
end do

call custom_print(0,2,"Initial reduced chi-squared:",csq/(m-n))

do i = 1,iterations

	! accumulate total Hessian matrix (alpha_kl in textbook)
	A = 0.0
	do j = 1,m
		A = A + ddf(xx(j),c,n)/sigma(j)**2
	end do
	Areg = regular(A,lambda,n)

	! accumulate chisq derivative vector (beta_k in textbook)
	B = 0.0
	do j = 1,m
		B = B - sqrt(chisq(xx(j),yy(j),sigma(j),c,n))*df(xx(j),c,n)/sigma(j)
	end do

	! invert regularlized Hessian matrix to solve for delta_c
	call lsolve(n,Areg,B)

	! evaluate new chi-squared
	dcsq = 0.0
	do j = 1,m
		dcsq = dcsq + chisq(xx(j),yy(j),sigma(j),c+B,n)
	end do

	! update lambda and c based on size of new chi-sq compared to
	! initial chi-sq, reject step if new chi-sq is bigger
	if (dcsq >= csq) then
		lambda = lambda*10
		call custom_print(0,0,"same reduced chisq =",csq/(m-n))
		cycle
	else
		lambda = lambda/10
		c = c + B
		csq = dcsq
		call custom_print(0,0,"next reduced chisq =",csq/(m-n))
		cycle
	end if

end do

! output data/parameters
call custom_print(2,0,"Final fit parameters:")
do i = 1,n
	call custom_print(0,0,"",c(i))
end do
call custom_print(2,0,"x,     y,     bestfit(x)")
do i = 1,m
	write(*,*) xx(i), yy(i), f(xx(i),c,n)
end do

deallocate(A,Areg,B,c,xx,yy,sigma)



contains

! computes fit function
pure function f(x,ca,n)
	intent(in) x, ca, n
	integer :: i, n
	real :: f, x, ca(n), basis

	basis = 0
	do i = 1,n-1
		basis = basis + ca(i)*cos(2*pi*(i-1)*x)
	end do
	f = exp(basis) + ca(n)
end function

! computes 1st derivative vector for single point
pure function df(x,ca,n)
	intent(in) x, ca, n
	integer :: i, n
	real :: df(n), x, ca(n)

	do i = 1,n-1
		df(i) = (f(x,ca,n) - ca(n))*cos(2*pi*(i-1)*x)
	end do
	df(n) = 1.0
end function

! computes 2nd derivative matrix for single point (minus 2nd term)
pure function ddf(x,ca,n)
	intent(in) x, ca, n
	integer :: i, j, n
	real :: ddf(n,n), dfdx(n), x, ca(n)

	dfdx = df(x,ca,n)
	do i = 1,n
		do j = 1,n
			ddf(i,j) = dfdx(i)*dfdx(j)
		end do
	end do
end function

! evaluates chisq for single point
pure function chisq(x,y,sigma,ca,n)
	intent(in) x, y, sigma, ca, n
	integer :: i, n
	real :: chisq, x, y, sigma, ca(n)

	chisq = ((y - f(x,ca,n))/sigma)**2
end function

! adds regularizing lambda parameter to Hessian
pure function regular(A,L,n)
	intent(in) A, L, n
	integer :: i, n
	real :: regular(n,n), A(n,n), L

	regular = A
	do i = 1,n
		regular(i,i) = regular(i,i)*(1.0 + L)
	end do
end function

! from Dr.Frolov's "linear.f90" code
! solve A.x = B using LAPACK xGESV
! A gets destroyed, answer is returned in B
subroutine lsolve(n, A, B)
	integer n, pivot(n), status; real A(n,n), B(n)

	! initialize status to all clear
	status = 0

	! call appropriate version of xGESV
	select case (kind(A))
		case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
		case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
		case default; call abort
	end select

	! abort at first sign of trouble
	if (status /= 0) stop "singular matrix in lsolve()"
end subroutine

! prints m lines and custom message + optional numerical value
subroutine custom_print(m,k,msg,val)
implicit none

	integer :: i, m, k
	character(len=*) :: msg
	real, optional :: val

	if (m >= 1) then
		do i = 1,m
			write(*,*) ""
		end do
	end if
	if (present(val)) then
		write(*,*) "#", msg, val
	else
		write(*,*) "#", msg
	end if
	if (k >= 1) then
		do i = 1,k
			write(*,*) ""
		end do
	end if
end subroutine

end program
