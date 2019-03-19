! general linear least squares fit
! compile with:
! gfortran -O3 -fdefault-real-8 LapackLeastSquares.f90 -llapack && ./a.out < Assignment\ #2.dat
! plot with:
! p'Assignment #2.dat' w d, 'fort.3' w l

program LapackLeastSquares
implicit none
real, parameter :: pi = 4.D0*DATAN(1.D0)
integer, parameter :: n = 4

real, allocatable :: A(:,:), B(:)
real x, y

integer i, j, k, counter
integer :: status = 0

counter = lines (7, "Assignment #2.dat")

allocate (A(counter, n), B(counter))


! Open the data file that we want to play with
open (unit=2, file="Assignment #2.dat")
! read data from standard input
i = 1
do while (status == 0)
	read (2,*,iostat=status) x, y
	if (status < 0) exit
	B(i) = y
	do j = 1, n
		A(i, j) = basis(x, j)
	end do
	i = i + 1
end do
close (2)

call solve (A, B, counter, n)

call dump (A, B(1:n), counter, n)

contains

! Find the number of lines in the file
function lines (window, filename)
	integer lines, window, status;	character(len=100) :: filename
	
	! Initiate counters
	lines = 0;	status = 0
	
	! Open the data file that we want to play with and count the lines of
	open (unit=window, file=filename)
	do while (status == 0)
		read (window,*,iostat=status)
		if (status < 0) exit
		lines = lines + 1
	end do
	close (window)
end function

! Evaluate the approximation at a certain x value
real function basis(x, k)
	integer k;	real x
	basis = cos(2.0*pi*(k-1.0)*x)
	return
end function


subroutine solve (A, B, M, N)
	integer M, N
	real A(M, N), B(M), C(M, N)
	! Out
	real S(N), W(6*M)
	integer R, status
	
	C = A
	
	call dgelss(M, N, 1, C, M, B, M, S, 1.0e-6, R, W, 6*M, status)
end subroutine

subroutine dump  (A, B, counter, n)
	real A(counter, n), B(n), C(counter)
	integer counter, n, i
	status = 0; i = 1
	
	C = matmul(A,B)
	
	rewind (2)
	open (unit = 2, file = "Assignment #2.dat")
	do while (status == 0)
		read (*,*,iostat=status) x
		if (status < 0) exit
		write (1,*) x, C(i)
		i = i + 1
	end do
	close (2)
	

end subroutine

end program
