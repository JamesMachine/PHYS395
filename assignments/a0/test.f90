! Running Instruction following:
! gfortran -o3 -fdefault-real-8 test.f90 -o test && ./test

!program start as following
program fibonacci
implicit none

  integer::n

  write(*,*) 'program start!'
  do n = 1,10
    write(*,*) f(n)
  end do


contains

  integer function f(n)
    implicit none
    integer,intent(in)::n
      f = n+1
  end function f

end program fibonacci

