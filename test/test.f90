! Assignment3 - Q2, Q3

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 test.f90 -llapack -o test && ./test

program a3code; implicit none

  !constants
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
	real, parameter :: eps = 1.0e-18 ! small number epsilon
  integer,parameter ::numpts = 1000

  !declarations
  real :: a(2),x,y
  !real,dimension(:) :: xapprox,yapprox


  integer ::i,status,lines
  character(len=8000) :: arg   ! buffer to store arguments in

  !read number of lines
  !call get_command_argument(1, arg)
  !read(arg,*,iostat=status) lines
  !write(*,*) lines

  open(1, file='data.dat',status = 'old')




  close(1)




contains

!it performs non-element wise operation --> gives non-scalar, so cannot be elemental, must be pure
  pure function populateArray(a,b,numpts)
    real,intent(in) :: a,b
    integer, intent(in) ::numpts
    integer ::i
    real ::populateArray(numpts)

    populateArray(1:numpts) = (/(a+(b-a)*i/(numpts-1),i=0,numpts-1)/)
  end function

  ! basis functions we are fitting
  subroutine evalb(x,B)
    real x, B(n)
    integer::alpha

    do alpha = 1,n
      B(alpha) =cos(2.0*pi*(real(alpha-1.0)*x))
    end do
  end subroutine


end program
