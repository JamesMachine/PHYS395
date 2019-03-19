! Assignment3 - Q4,Q5 - Gradient Descent method

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 a3_nonlinear_fit_gradDec.f90 -llapack -o a3_nonlinear_fit_gradDec && ./a3_nonlinear_fit_gradDec

program a3_nonlinear_fit_gradDec; implicit none

  !constants
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
	real, parameter :: eps = 1.0e-18 ! small number epsilon
  integer,parameter ::numpts = 1000, n=5, m=9130
  real,parameter :: dt = 0.1
  !numpts: number of points for the approximation to draw
  !n: number of coefficients (# of terms + const coeffi)
  !m: number of sample points
  !r: residuals vector

  !declarations
  integer::i,status, iterations
  real ::x(m),y(m),c(n),cnew(n),r(m),rnew(m), J(m,n), g(n,n), E(n,n), gradChi(n),delta(n)
  real ::chi,chinew
  !real,dimension(:) :: xapprox,yapprox
  !n_modified: number of terms + number of const added (in our case, 3 terms + 1 const = 4)


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

  !initial guess of coefficients vector c, lambda
  c(:) = [1.0, 1.0, 1.0, 1.0, 1.0]

  !initialize residual vector
  do i=1,m
    r(i) = y(i) - f(x(i),c)
  end do

  !print initial chi
  chi = 0.5 * dot_product(r,r)
  write(*,*) "initial chi",chi

  !initialize jacobian matrix
  J = jacobian(m,n,x,c)

  do iterations=1,100

    !evaluate the gradient of chi
    gradChi = matmul(transpose(J),r)

    !evaluate the chi
    chi = 0.5 * dot_product(r,r)

    !solve for delta
    delta = gradChi*dt

    ! calculate new coefficients
    cnew = c - delta

    ! calculate new residuals
    do i=1,m
      rnew(i) = y(i) - f(x(i),cnew)
    end do

    ! calculate new chi square
    chinew = 0.5 * dot_product(rnew,rnew)

    c = cnew
    r = rnew
    write(*,*) chi, chinew

  end do

  write(*,*) "final chi", chi
  write(*,*) "cofficients", c

contains

!it performs non-element wise operation --> gives non-scalar, so cannot be elemental, must be pure
  pure function populateArray(a,b,numpts)
    real,intent(in) :: a,b
    integer, intent(in) ::numpts
    integer ::i
    real ::populateArray(numpts)

    populateArray(1:numpts) = (/(a+(b-a)*i/(numpts-1),i=0,numpts-1)/)
  end function

  function basis(xi,alpha); intent(in)::xi,alpha
    real::xi
    integer::alpha
    real::basis
    basis = cos(2*pi*(alpha-1)*xi) ! index is adjust to (alpha-1) because Fortran starts index at 1
  end function

  function f(xi,c); intent(in)::xi,c
    real::xi,c(n),f, sigma
    integer ::alpha

    sigma = 0
    do alpha = 1,n-1 ! alpha is 1 to n-1 because Fortran starts index from 1 not 0 (note basis is adjusted)
      sigma = sigma + c(alpha) * basis(xi,alpha)
    end do
    f = exp(sigma) + c(n) !add coefficient at the last index, which is const
  end function

  !partial derivative at alpha index
  function pdf(xi,c,pd_at_alpha); intent(in)::xi,c,pd_at_alpha
    real::xi,c(n),pdf, sigma
    integer ::alpha,pd_at_alpha

    sigma = 0
    if (pd_at_alpha /= n) then
      do alpha = 1,n-1
        sigma = sigma + c(alpha) * basis(xi,alpha)
      end do
      pdf = basis(xi,pd_at_alpha) * exp(sigma)
    else
      pdf = 1.0
    end if

  end function


  function jacobian(m,n,x,c); intent(in)::m,n,x,c
    integer::m,n,i,alpha
    real::x(m),c(n)
    real::jacobian(m,n)

    do alpha=1,n
      do i=1,m
        jacobian(i,alpha) = -pdf(x(i),c,alpha)
      end do
    end do
  end function

end program
