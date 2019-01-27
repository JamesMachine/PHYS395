! Running Instruction following:
! gfortran -o3 -fdefault-real-8 rec_fib.f90 -o rec_fib && ./rec_fib

!program start as following
program rec_fib
  implicit none
  integer::n



  write(*,*) 'program start!!----------------'
  do n=0,20
    write(*,*) 'f(',n,')=',f(n)
  end do


contains

  recursive integer function f(n) result(res)
    implicit none
    integer, intent(in)::n
    if(n<=1) then
      res = n
    else
      res = f(n-1) + f(n-2)
    end if

  end function f

end program rec_fib

