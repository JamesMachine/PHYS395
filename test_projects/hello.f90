program hello
implicit none

! number of points to print
integer :: i=& !go to next line
10
real a,b(3)

i = 10
i=i+3.0
a = i/3.0
b = [1.0,2.0,3.0] + 10.0

do i=1,3
  if (mod(int(b(i)),2)==1) then
    write(*,*) "'''"
  else
    write(*,*) "???"
  end if
  write(*,'(g32.3)') b(i)
end do

!forall write b(i)


!write (*,*) "hello, world", i/3, a
!write (*,*) sum(b)

!number of points to print
integer, parameter :: n =100
real x(0:n), y(0:n)

!x = [0:n]/real(n)
!y = sin(10.0*x)

forall (i=0:n) x(i)=i/real(n)

do i=0,n
  write(*,*) x(i), y(i)
end do




end
