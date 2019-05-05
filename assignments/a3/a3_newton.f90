! Assignment3 - Q1

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 a3_newton.f90 -llapack -o a3_newton && ./a3_newton


program a3_newton; implicit none

  !constants
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
	real, parameter :: eps = 1.0e-18 ! small number epsilon
	real, parameter :: phi = (sqrt(5.0) + 1.0) / 2.0 !Golden Ratio: 1.618...

  integer,parameter ::numpts = 100


  !declarations
  real ::x,xexpect(numpts),yexpect(numpts)
  integer ::i

  !write f(x) into output1
  open(1,file='output1')
  xexpect = populateArray(-2.0,2.0,numpts)
  do i=1,numpts
    yexpect(i) = f(xexpect(i))
    write(1,*) xexpect(i),yexpect(i)
  end do
  close(1)

  !write roots into output1_roots
  open(1,file='output1_roots')
  x= newton(1.0)
  write(1,*) x, f(x)
  x= newton(0.0)
  write(1,*) x, f(x)
  x= newton(-1.0)
  write(1,*) x, f(x)
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

  ! function to find a root of...
  elemental function f(x); intent(in) x
    real f, x
    f = x ** 3 - x + 1.0 / 4.0
  end function

  ! derivative of a function to find a root of
  elemental function df(x); intent(in) x
    real df, x
    df = 3 * x ** 2 - 1
  end function

  ! newton's method function
  elemental function newton(x0)
    real,intent(in) ::x0 ! intent(in) of x0 to make x0 cannot be changed
    real::x !making copy of x but not global variable x
    real::newton ! you need this for output

    ! initial guess
    x=x0
    do while (abs(f(x)) > eps)
      x = x - f(x) / df(x)
    end do
    newton = x ! result to be returned
  end function




end program
