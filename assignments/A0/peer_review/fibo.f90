program fibo
implicit none
! NOTE: Iteration 47 and above are erroneous, becoming negative and incorrect values.



! Controls how many points are calculated. Must be at least 2 iterations.
integer, parameter :: n = 50

! Creates blank matrix to input values.
integer y(0:n)



! Creates variables to be used.
integer :: a, b, i

a = 0

b = 1

y(0) = a
y(1) = b


! Prints values in same format as the rest.
write (*,*) y(0)
write (*,*) y(1)


! Creates and prints value of each fibonacci iteration.
do i = 2,n

y = a + b


a = b

b = y(i)

  write (*,*) y(i)

end do


end


! Creator: Darryl Arab
!STD.Num:  301289807
! Assignment # 0
