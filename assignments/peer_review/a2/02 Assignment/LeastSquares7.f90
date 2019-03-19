! general linear least squares fit
! compile with:
! gfortran -O3 -fdefault-real-8 LeastSquares7.f90 -llapack && ./a.out < Assignment\ #2.dat
! plot with:
! p'Assignment #2.dat' w d, 'fort.1' w l

program leastsq
implicit none

! number of basis functions, which run from 0 to n-1
integer, parameter :: n = 8
real, parameter :: pi = 4.D0*DATAN(1.D0)

real x, y, A(n,n), B(n), U(n)

integer i, j
integer :: status = 0

! initialize accumulators
A = 0.0
U = 0.0

! Open the data file that we want to play with
open (unit=2, file="Assignment #2.dat")
! read data from standard input
do while (status == 0)
	read (2,*,iostat=status) x, y
	if (status < 0) exit
	! evaluate basis functions at point x
	call evalb(x, B, n)
	
	! accumulate least square problem matrices
	forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j)
	U = U + y*B
end do
close (2)

write (*,*) "For n = ", n-1

! solve for best fit parameters
!call lsolve(n, A, U)
call svdsolve(n, A, U, 1.0e-6)

! output best fit parameters
write (*,*) "The best fit coefficients are:", U

call dump (U, n)

call chisq (U, n)

contains

! basis functions we are fitting
subroutine evalb(x, B, n)
	real x, B(n);	integer i, n

	! Evalute the cosine functions
	forall (i=1:n) B(i) = cos(2*pi*(i-1)*x)
end subroutine

! solve A.x = B using LAPACK xGESVD
! A gets destroyed, answer is returned in B
subroutine svdsolve(n, A, B, epsilon)
	integer n, status; real A(n,n), B(n), epsilon
	
	! SVD matrices
	real U(n,n), V(n,n), S(n), W(6*n)
	
	! initialize status to all clear
	status = 0
	
	! call appropriate version of xGESV
	select case (kind(A))
		case(4); call sgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
		case(8); call dgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
		case default; call abort
	end select
	
	! abort at first sign of trouble
	if (status /= 0) stop "singular matrix in svdsolve()" 

	write (*,*) "The condition number is:", S(1)/S(n)
	
	! compute the solution using pseudo-inverse
	B = matmul(transpose(U),B)
	where (S > epsilon*S(1)); B = B/S; elsewhere; B = 0.0; end where
	B = matmul(transpose(V),B)
end subroutine

! Write results of fit
subroutine dump (U, n)
	integer n, k, l;	integer :: g = 10000;	real U(n), y
	
	do k = 1,g
		x = (k-1.0)*1.0/(g-1.0)
		y = f (x, U, n)
		write (1,*) x, y
	end do
	write (1,*) "";	write (1,*) ""
end subroutine

subroutine chisq (U, n)
	integer n, k, l, counter;	real U(n), x, y, chi, fox, q
	
	! Must reset status from when it failed in the main body
	status = 0;	chi = 0.0;	counter = 0
	! Since we use the same file as in the body, we need to
	! start reading from its beginning again
	rewind (2)
	open (unit = 2, file = "Assignment #2.dat")
	do while (status == 0)
		read (*,*,iostat=status) x, y
		if (status < 0) exit
		counter = counter + 1
		! evaluate chi squared at point x
		fox = f (x, U, n)
		chi = chi + (y-fox)**2.0
	end do
	close (2)
	write (*,*) "Chi squared is:", chi
	write (*,*) "Reduced chi squared is:", chi/(counter-(n+1))
end subroutine

! Evaluate the approximation at a certain x value
real function f(x, U, n)
	integer i, n;	real x, U(n)
	f = 0
	do i = 1,n
		f = f + cos(2.0*pi*(i-1.0)*x)*U(i)
	end do
	return
end function


end program
