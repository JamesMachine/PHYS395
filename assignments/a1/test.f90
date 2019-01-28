!Running Instruction
!gfortran -fdefault-real-8 test.f90 -o test && ./test

program test

  implicit none
  integer:: i,n

  real,dimension(:),allocatable::x,y




  n=10
  allocate(x(n),y(n))

  x=populateArray(-1.0,1.0,n)
  do i=1,n
    print*, x(i), f(x(i))
  end do

  deallocate(x,y)


contains


  elemental function f(x)
    real, intent(in) :: x
    real ::f

    f=(1.0+10.0*x**2.0)**(-1.0)
  end function


  pure function populateArray(a,b,n)
    real,intent(in) :: a,b
    integer, intent(in) ::n
    integer ::i
    real ::populateArray(n)

    populateArray(1:n) = (/(a+(b-a)*i/real((n-1)),i=0,n-1)/)
  end function


end program
