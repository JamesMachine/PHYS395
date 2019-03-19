! argv.f90  -  command line arguments in Fortran 2003
! compile with: gfortran -O3 test.f90 -o test && ./test


program test; implicit none

  !constants
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
	real, parameter :: eps = 1.0e-18 ! small number epsilon
  integer,parameter ::numpts = 1000, n=3, m=9130
  !numpts: number of points for the approximation to draw
  !n: number of terms + number of const added (in our case, 3 terms + 1 const = 4)
  !m: number of sample points
  !r: residuals vector

  !declarations
  real ::x(m),y(m),c(n+1),r(m),I(n,n)
  integer::i,status, lamda,decreaseFactor,increaseFactor,j
  !real,dimension(:) :: xapprox,yapprox


  integer ::lines
  !character(len=8000) :: arg   ! buffer to store arguments in

  !read number of lines
  !call get_command_argument(1, arg)
  !read(arg,*,iostat=status) lines
  !write(*,*) lines

  !read data points
  open(1, file='data.dat',status = 'old')
  do i=1,m
    read (1,*,iostat=status) x(i), y(i)
    if (status < 0) exit
  end do
  close(1)

  !initial guess of coefficients vector c, lamda
  c(:) = [1.0, 1.0, 1.0, 1.0, 1.0]

  lamda = 1.0
  decreaseFactor = 1.5
  increaseFactor = 1.5

  I=indentity(n)
  do i=1,n
   do j=1,n
    write(*,*) I(i,j)
   end do
  end do








contains

!it performs non-element wise operation --> gives non-scalar, so cannot be elemental, must be pure
  pure function populateArray(a,b,numpts)
    real,intent(in) :: a,b
    integer, intent(in) ::numpts
    integer ::i
    real ::populateArray(numpts)

    populateArray(1:numpts) = (/(a+(b-a)*i/(numpts-1),i=0,numpts-1)/)
  end function

  function identity(n); intent(in)::n
    integer::n,i,j
    real::identity(n,n)
    do i=1,n
      indentity(i,i) = 1.0
    end do
  end function



end program
