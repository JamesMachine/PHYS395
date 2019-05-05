! q2

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 q2.f90 -llapack -o q2 && ./q2


program q2; implicit none

  !constants
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
	real, parameter :: eps = 1.0e-18 ! small number epsilon
	real, parameter :: phi = (sqrt(5.0) + 1.0) / 2.0 !Golden Ratio: 1.618...

  integer,parameter ::numpts = 1000


  !declarations
  real ::x,xexpect(numpts),yexpect(numpts)
  integer ::i


  write(*,*) "q2 started..."

  !write P(x)
  open(1,file='./output/output_q2/function.dat')
  xexpect = populateArray(-20.0,20.0,numpts)
  do i=1,numpts
    yexpect(i) = f(xexpect(i))
    write(1,*) xexpect(i),yexpect(i)
  end do
  close(1)

  write(*,*) "extrema:"
  write(*,*) "(x,f(x))"
  open(1,file='./output/output_q2/extrema.dat')
    x = bracket_min(-4.0, -2.0)
    write(1,*) x, f(x)
    write(*,*) x, f(x)
    x = bracket_max(-2.0, 0.0)
    write(1,*) x, f(x)
    write(*,*) x, f(x)
    x = bracket_min(0.0, 2.0)
    write(1,*) x, f(x)
    write(*,*) x, f(x)
  close(1)

  write(*,*) ""
  write(*,*) ""

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
    f = x**4 + 3*x**3 - 4*x**2 - 3*x + 4
  end function


  ! bracketed minimum search function
  elemental function bracket_min(a0,b0); intent(in)::a0,b0 ! you can also do intent(in) here
    real::a0,b0,a,b,d,xleft,xright
    real::bracket_min

    a = a0 !store into a and b, since you cannot modify a0 and b0
    b = b0

    d = (b - a) / phi
    xleft = b - d
    xright = a + d
    do while (abs(xright - xleft) > eps)
      if (f(xleft) < f(xright)) then
        b = xright ! b becomes xright
      else
        a = xleft  ! a becomes xleft
      end if

      ! after narrowing [a,b], update xleft and xright
      d = (b - a) / phi
      xleft = b - d
      xright = a + d
    end do

    bracket_min = (a + b) / 2
  end function

  ! bracketed maximum search function
  elemental function bracket_max(a0,b0); intent(in)::a0,b0 ! you can also do intent(in) here
    real::a0,b0,a,b,d,xleft,xright
    real::bracket_max

    a = a0 !store into a and b, since you cannot modify a0 and b0
    b = b0

    d = (b - a) / phi
    xleft =  a + d
    xright = b - d
    do while (abs(xright - xleft) > eps)
      if (f(xleft) < f(xright)) then
        b = xright ! b becomes xright
      else
        a = xleft  ! a becomes xleft
      end if

      ! after narrowing [a,b], update xleft and xright
      d = (b - a) / phi
      xleft =  a + d
      xright = b - d
    end do

    bracket_max = (a + b) / 2
  end function


end program

