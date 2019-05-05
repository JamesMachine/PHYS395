! Final Exam - root finding q1

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 q1.f90 -llapack -o q1 && ./q1


program q1; implicit none

  !constants
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
	real, parameter :: eps = 1.0e-10 ! small number epsilon
  real, parameter :: bisect_interval = 0.5
  integer,parameter ::numpts = 1000


  !declarations
  real :: x, xexpect(numpts), yexpect(numpts)
  real :: a, b

  integer ::i

  !output files
  open(1,file='./output/output_q1/function.dat')
  open(2,file='./output/output_q1/roots.dat')

  write(*,*) "q1 started..."

  !Step1: Write f(x) into function.dat
  xexpect = populateArray(-20.0,20.0,numpts)
  do i=1,numpts
    yexpect(i) = f(xexpect(i))
    write(1,*) xexpect(i),yexpect(i)
  end do

  !Step2: interations with 0.5 interval to find roots
  a = -20.0
  b = a + bisect_interval

  write(*,*) "roots:"
  write(*,*) "(x,f(x))"
  do i=1,100
    x = bisect(a, b)
    ! If found bracketed value, write into file roots.dat
    if (x/=-9999.0) then
      write(*,*) x, f(x)
      write(2,*) x, f(x)
    end if

    a = a+bisect_interval
    b = b+bisect_interval
  end do

  close(1)
  close(2)

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
    f = cos(x) - x / 5.0
  end function



  ! derivative of a function to find a root of
  !elemental function df(x); intent(in) x
  !  real df, x
  !  df = sin(x) - 1.0 / 5.0
  !end function




  function bisect(a0,b0); intent(in) ::a0,b0
    real ::a0, b0, a, b, c, fa, fb, fc
    integer ::i
    real ::bisect


    ! initial interval
    a = a0; fa = f(a)
    b = b0; fb = f(b)
    if (fa*fb > 0.0) then
      !write(*,*) "The root is not bracketed, bailing out..."
      bisect = -9999.0 ! -9999.0 is a signal to indicate that root is not found
      return
    end if

    fc = 1 ! initialize fc to 1, because fortran contains small value in fc before initializing it
    ! bisect the interval
    do while(abs(fc) > eps)
      c = (a+b)/2.0; fc = f(c); if (fc == 0.0) exit

      ! Ridder's variation on the basic method
      c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(c)

      if (fa*fc < 0.0) then; b = c; fb = fc; end if
      if (fc*fb < 0.0) then; a = c; fa = fc; end if
    end do

    bisect = c
  end function



end program
