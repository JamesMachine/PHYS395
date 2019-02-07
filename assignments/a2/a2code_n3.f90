!Assignment 2

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 a2code_n3.f90 -llapack -o a2code_n3 && ./a2code_n3


program a2code
  implicit none

  ! number of basis functions & constants
  integer, parameter :: n = 4 !(n=0, 1, 2, 3)
  integer,parameter :: numpts = 1000
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0


  real x, y, yexpect, A(n,n), B(n), U(n) !this is professor way to use do while loop
  real::xapprox(numpts),yapprox(numpts)
  real::condnum, chisq,counter,degOfFreedom


  integer i, j, k
  integer :: status

  ! initialize accumulators
  A = 0.0
  U = 0.0

  ! read data from standard input
  open(1,file='Assignment #2.dat',status='old')
  do while (status == 0) ! loop through each point x and y (or say xi)
    read (1,*,iostat=status) x, y
    if (status < 0) exit

    ! evaluate basis functions at point x
    call evalb(x, B)

    ! accumulate least square problem matrices
    forall (k=1:n,j=1:n) A(k,j) = A(k,j) + B(k)*B(j) !every time add Bk*Bj in the matrix A
    !U = U + y*B
    forall(j=1:n) U(j)=U(j)+y*B(j) !everytime, add y*Bj in the vector U(rho)
  end do
  close(1)

  ! output accumulated matrices
  !write (*,*) A
  !write (*,*) U !rho


  !solve for best fit parameters --> get corefficients U
  !call lsolve(n, A, U)
  call svdsolve(n, A, U,condnum, 1.0e-6)

  ! output best fit parameters
  !write (*,*) U !coefficient vec


  xapprox = populateArray(0.0,1.0,numpts)

  open(1,file='output1.dat',status='unknown')
  call fapprox_array(U,xapprox,yapprox)
  do i=1,numpts
    write(1,*) xapprox(i),yapprox(i)
  end do
  close(1)

  write(*,*) 'Condition number is ',condnum
  write(*,*) 'Best fit parameter is:',U

  !Chi-squared test
  chisq = 0
  counter = 0
  status = 0
  open(1, file='Assignment #2.dat',status = 'old')
  do while(status == 0)
    read(1,*,iostat=status) x,y
    if (status < 0) exit

    call fapprox(U,x,yexpect)
    chisq = chisq + abs((y - yexpect)**2)
    counter = counter + 1
    !write(*,*) x,y,yexpect
  end do
  close(1)

  degOfFreedom = counter - n
  write(*,*) 'Degree of Freedom is: ', degOfFreedom
  write(*,*) 'Chi Squared is: ', chisq
  write(*,*) 'Error % is:', abs((chisq-degOfFreedom)/degOfFreedom)*100




contains

  ! basis functions we are fitting
  subroutine evalb(x,B)
    real x, B(n)
    integer::alpha

    do alpha = 1,n
      B(alpha) =cos(2.0*pi*(real(alpha-1.0)*x))
    end do
  end subroutine


  ! solve A.x = B using LAPACK xGESV
  ! A gets destroyed, answer is returned in B
  subroutine lsolve(n, A, B)
    integer n, pivot(n), status; real A(n,n), B(n)

    ! initialize status to all clear
    status = 0

    ! call appropriate version of xGESV
    select case (kind(A))
      case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
      case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
      case default; call abort
    end select

    ! abort at first sign of trouble
    if (status /= 0) stop "singular matrix in lsolve()"
  end subroutine



  ! solve A.x = B using LAPACK xGESVD
  ! A gets destroyed, answer is returned in B
  subroutine svdsolve(n, A, B,condnum, epsilon)
    integer n, status; real A(n,n), B(n), epsilon

    ! SVD matrices
    real U(n,n), V(n,n), S(n), W(6*n)
    real ::condnum

    ! initialize status to all clear
    status = 0

    ! call appropriate version of xGESV
    select case (kind(A))
      case(4); call sgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
      case(8); call dgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
      case default; call abort
    end select

    ! abort at first sign of trouble
    if (status /= 0) stop "singular matrix in svdsolve()"

    ! compute the solution using pseudo-inverse
    ! B is rho
    B = matmul(transpose(U),B)
    !B is U^-1 * rho
    where (S > epsilon*S(1)); B = B/S; elsewhere; B = 0.0; end where
    !B is S^-1 * U^-1 * rho
    B = matmul(transpose(V),B)
    !B = V * S^-1 * U^-1 * rho

    !now computer the condition number
    condnum = maxval(S)/minval(S)

  end subroutine

  pure function populateArray(a,b,numpts)
    real,intent(in) :: a,b
    integer, intent(in) ::numpts
    integer ::i
    real ::populateArray(numpts)

    populateArray(1:numpts) = (/(a+(b-a)*i/(numpts-1),i=0,numpts-1)/)
  end function

  !approximation using array
  subroutine fapprox_array(c,xapprox,yapprox)
    real::c(n),xapprox(numpts),yapprox(numpts)
    integer::i,alpha
    !forall(i=1:numpts,j=1:n) yapprox(i) = yapprox(i) + c(j)*cos((j-1)*acos(xapprox(i))) !this is syntax error

    do i=1,numpts
      yapprox(i)=0
      do alpha=1,n
        yapprox(i) = yapprox(i) + c(alpha)*cos(2*pi*(alpha-1)*xapprox(i))
      end do
    end do
  end subroutine

  !approximation at one value
  subroutine fapprox(c,xapprox,yapprox)
    real::c(n),xapprox,yapprox
    integer::alpha

    yapprox=0
    do alpha=1,n
      yapprox = yapprox + c(alpha)*cos(2*pi*(alpha-1)*xapprox)
    end do
  end subroutine





end program
