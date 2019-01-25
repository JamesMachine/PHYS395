!use the following command to compile:
!gfortran -O3 -fdefault-real-8 Fibonacci.f90 -o Fibonacci  && ./Fibonacci > OUTPUT
!outputs to text file "OUTPUT"

program fibonacci
	!calculates and outputs the first 50 fib numbers using an array
	!starting from 0,1... (so first number defined as 0)
	integer :: max_n = 50
	integer*8, dimension(50) :: fib_numbers !needs *8 or last couple numbers overflow
	integer :: n = 3 !index starts from 3 since 1 and 2 are hardcoded

	!hardcode for first 2 fib numbers
	fib_numbers(1) = 0
	write (*,*) fib_numbers(1)
	fib_numbers(2) = 1
	write (*,*) fib_numbers(2)

	!loop through the array and calculate each fib number
	do
		fib_numbers(n) = fib_numbers(n-1) + fib_numbers(n-2)
		write (*,*) fib_numbers(n)
		if (n == max_n) exit
		n = n+1
	end do
  
end program
