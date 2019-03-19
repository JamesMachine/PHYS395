! Assignment3 - Q2, Q3

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 a3_bracket.f90 -llapack -o a3_bracket && ./a3_bracket


program a3code; implicit none

  !constants
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
	real, parameter :: eps = 1.0e-18 ! small number epsilon
	real, parameter :: phi = (sqrt(5.0) + 1.0) / 2.0 !Golden Ratio: 1.618...

  integer,parameter ::numpts = 100


  !declarations
  real ::x,xexpect(numpts),yexpect(numpts)
  integer ::i


  !Q2, Q3
  !write f(x) into output2
  open(1,file='output2')
  xexpect = populateArray(-2.0,2.0,numpts)
  do i=1,numpts
    yexpect(i) = f(xexpect(i))
    write(1,*) xexpect(i),yexpect(i)
  end do
  close(1)


  open(1,file='output2_minima')
    x = bracket(-1.5,-0.5)
    write(1,*) x, f(x)
    x = bracket(0.5 , 1.5)
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
    f = (x**2 - 1) ** 2 + x
  end function


  ! bracketed minimum search function
  elemental function bracket(a0,b0); intent(in)::a0,b0 ! you can also do intent(in) here
    real::a0,b0,a,b,d,xleft,xright
    real::bracket

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

    bracket = (a + b) / 2
  end function


end program
