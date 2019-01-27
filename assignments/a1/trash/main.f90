!Assignment 1

!Running Instruction
!gfortran -fdefault-real-8 main.f90 -o main && ./main

program main

  use tools
  implicit none



  !Question 1
  integer, parameter :: n = 3
  integer ::i


  integer ::terms
  real ::x(100)

  real A(n,n), B(1,n)

  A(1,:) = [0.0, 1.0, 3.0]
  A(2,:) = [3.0, 300.0, 9.0]
  A(3,:) = [4.0, 4.0, 1.0]

  B(1,1) = 5.0
  B(1,2) = 600.0
  B(1,3) = 8.0


  call printMatrix(n,n,A)
  call printMatrix(n,1,B)
  call gaussj(n,A,B)
  call printMatrix(n,n,A)
  call printMatrix(n,1,B)



end program
