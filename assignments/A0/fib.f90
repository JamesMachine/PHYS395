!Assignment 0
!Coder: Takehiro Tanaka

! Running Instruction following:
! gfortran -o3 -fdefault-real-8 fib.f90 -o fib && ./fib

program fib
  implicit none

  integer(kind = 8)::n,fibnum,temp1,temp2 ! declare  8 byte (64-bit) variables
  ! fibnum: f(n)
  ! temp1: f(n-2)
  ! temp2: f(n-1)

  integer::maxN =  50 ! choose maximum n to be calculate for fibnum
  fibnum = 0
  do n = 0,maxN
    if (n<=1) then
      temp1 = 0 !set temp1 always 0 for n<=1
      temp2 = 1 !set temp2 always 1 for n<=1
      fibnum = n
    else
      fibnum = temp1 + temp2
      temp1 = temp2
      temp2 = fibnum
    end if

    write(*,*) n, fibnum ! print the result: n and fibnum
  end do

end program fib

