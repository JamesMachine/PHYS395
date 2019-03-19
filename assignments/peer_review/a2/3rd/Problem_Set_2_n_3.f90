! This program will fit to a set of data points subjected
! to potentially noisy data points {x-i, y-i}

! Compile	: gfortran -O3 -fdefault-real-8 Problem_Set_2_n_3.f90 -llapack && ./a.out < "Assignment #2.dat"
! Run		: ./a.out

program Data_Fitting_3
implicit none

! Initialize parameters

! Allocate memory for size --> n+1
integer, parameter :: allocation = 1024, n = 3
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

! Size of array
integer :: input_size = 0, array_counter = 1, modcounter = 1, i = 1
integer j, k, status

! For initializing and allocation of array
real s(n+1,n+1), c(n+1), chi_sq, R_chi_sq
real,dimension(:),allocatable :: x, y, B, evaluated_func
real,dimension(:,:),allocatable :: A
allocate(x(allocation), y(allocation), B(allocation), evaluated_func(allocation), A(allocation, n+1))

! Read data from input file and count number of data points
do
	read (*,*,iostat=status) x(i), y(i)
	if (status < 0) exit

  	input_size = input_size + 1; i = i + 1

	! Keeping track of size allocation
  	if (modulo(input_size,allocation) == 0) then
		modcounter = modcounter + 1
    		call resize_array
  	end if
end do

! Finalize the array to numbers
call finalize_array

! Generate array of 3 basis functions for each row of A
do i=1,input_size; A(i,:) = basis(n, x(i)); end do

! Solve using lapack
call lsrsolve(input_size,n,A)

! Gather the coefficients from B to c
do i=1,n+1; c(i) = B(i); end do

! Print the basis coefficients
write (*,*)
write (*,*)	'***** Basis functions output *****'
write (*,*)

do i=1,n+1; write (*,*) 'Coefficient', i, '=', c(i); end do

! Calculate the goodness of fit here
call analyze(input_size, n, x, y, c)


contains


!  Evaluate basis function B_i(x) at point x
pure function basis(n, x)

	integer k, n; real x, basis(n+1); intent(in) n, x

	forall (k=1:n+1) basis(k) = cos(2*pi*(k-1)*x)

end function


! This will solve for the x that minimizes chi^2
subroutine lsrsolve(input_size, n, A)

	integer input_size, n, rank, info, lwork; real A(input_size,n+1), y(input_size), s(n+1,n+1), work
  	lwork = 3*input_size

  	call dgelss(input_size, n+1, 1, A, input_size, B, input_size, s, -1, rank, work, lwork, info)
  	if (info /= 0) then; call abort(); end if

end subroutine


subroutine resize_array
  	real,dimension(:),allocatable :: x_temp,y_temp

  	allocate(x_temp(modcounter*allocation), y_temp(modcounter*allocation))
  	x_temp( 1:size(x)) = x; y_temp(1:size(y))=y
  	deallocate(x, y)
  	allocate(x(size(x_temp)), y(size(y_temp)))
  	x = x_temp; y = y_temp

end subroutine resize_array


! Will finalize the array and resize to number of data points
subroutine finalize_array
	real,dimension(:),allocatable :: x_temp, y_temp

  	allocate(x_temp(input_size), y_temp(input_size))
  	x_temp(1:input_size) = x; y_temp(1:input_size) = y
  	deallocate(x,y)
  	allocate(x(size(x_temp)), y(size(y_temp)))

  	x = x_temp; y = y_temp

	! Resize A and B
	deallocate(A, B, evaluated_func)
	allocate(A(input_size, n+1), B(input_size), evaluated_func(input_size))

  	B = y

end subroutine finalize_array


! Here we will obtain value for approximated x(i) and the corresponding chi
subroutine analyze(input_size, n, x, y, c)! , evaluated_func)
  	integer input_size, n, i; real x(input_size), y(input_size), evaluated_func(input_size), c(n+1)
  	chi_sq = 0.0
  	do i = 1, input_size

    	! Evaluate the function at x(i)
    	evaluated_func(i) = sum(c*basis(n,x(i)))

    	! Obtain chi squared value by adding recursively
    	chi_sq = chi_sq + (y(i) - evaluated_func(i))**2

    	! Write to file 'fort.3'
    	write (3,*) x(i), y(i), evaluated_func(i)

  	end do

	! calculate reduced chi-squared
  	R_chi_sq = chi_sq/(input_size-1)

  	write (*,*) '***** Chi Squared *****'
  	write (*,*)
  	write (*,*) 'Value chi-squared = ' , chi_sq
	write (*,*)
  	write (*,*) '***** Reduced Chi Squared *****'
  	write (*,*)
  	write (*,*) 'Number of points used = ', input_size
  	write (*,*)
 	write (*,*) 'Value reduced-chi-squared = ', R_chi_sq
  	write (*,*)
	write (*,*) '"Good Fit" dictates that value of chi-squared must be on the order of 1'
	write (*,*)

end subroutine

end program
