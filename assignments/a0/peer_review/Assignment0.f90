! To compile script: gfortran -O3 -fdefault-integer-8 Assignment0.f90 -o fib
! To run with output in terminal: ./fib
! To run with output in file: ./fib > Fib50

program fibonacci
implicit none

! number of elements in sequence (n)
integer, parameter :: n=50

integer i ! loop variable

! empty array which will be filled with Fibonacci numbers with n elements
integer x(0:n) ! array of integers with indices from 0 to 50


! initialize first two numbers of the sequence
x(0)=0
x(1)=1
! print first two numbers of the sequence
write (*,*) x(0)
write (*,*) x(1)

! loop over
do i = 2,n
    x(i)=x(i-1)+x(i-2)
    write (*,*) x(i)
end do

end
