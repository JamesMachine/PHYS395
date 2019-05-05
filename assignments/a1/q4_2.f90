!Assignment 1 - q4 - 100 terms

!Running Instruction for this specific file
!gfortran -fdefault-real-8 q4_2.f90 -o q4_2 && ./q4_2


program a1code
implicit none

!declaration

integer:: i,n, k
!n: number of samples
!k: number of output points

real,dimension(:),allocatable::x,y,xapprox,yapprox,xexpect,yexpect,error


!choose number of samples
n=100

!choose output points as approximation
k=1000

allocate(x(n),y(n),xapprox(k),yapprox(k),xexpect(k),yexpect(k),error(k))

x=populateArray(-1.0,1.0,n)
xapprox=populateArray(-1.0,1.0,k)
xexpect=populateArray(-1.0,1.0,k)

write(*,*)
write(*,*) 'Approximating 100 sample points using 100 terms (q4) -------------------'
!plotting samples
open(1, file = 'output4_1_2.csv', status='unknown')
do i=1,n
  y(i)=f(x(i))
  write(1,*)  x(i), y(i)
end do
close(1)

! plotting expected curve
open(1, file = 'output4_2_2.csv', status='unknown')
do i=1,k
  yexpect(i) = f(xexpect(i))
  write(1,*) xexpect(i), yexpect(i)
end do
close(1)

! plotting approximated points
open(1, file = 'output4_3_2.csv', status='unknown')
yapprox = fapprox(x,y,xapprox,n,k)
do i=1,k
  write(1,*) xapprox(i), yapprox(i)
end do
close(1)

!plotting error taking abolute value of expected values and approximated values
open(1, file = 'output4_4_2.csv', status='unknown')
do i=1,k
  error(i)=abs(yexpect(i)-yapprox(i))
  write(1,*) xexpect(i), error(i)
end do
close(1)
write(*,*) " 1. max error for 100 terms:" ,maxval(error)
write(*,*) " 2. x value at the max error value for 100 terms:" ,xexpect(maxloc(error))

! now doing derivative -----------------------
! plotting expected curve
open(1, file = 'output4_5_2.csv', status='unknown')
do i=1,k
  yexpect(i) = df(xexpect(i))
  write(1,*) xexpect(i), yexpect(i)
end do
close(1)

! plotting approximated derivative points
open(1, file = 'output4_6_2.csv', status='unknown')
yapprox = dfapprox(x,y,xapprox,n,k)
do i=1,k
  write(1,*) xapprox(i), yapprox(i)
end do
close(1)

!plotting error taking abolute value of expected values and approximated values
open(1, file = 'output4_7_2.csv', status='unknown')
do i=1,k
  error(i)=abs(yexpect(i)-yapprox(i))
  write(1,*) xexpect(i), error(i)
end do
close(1)
write(*,*) " 3. max error for 10 terms (derivative):" ,maxval(error)
write(*,*) " 4. x value at the max error value for 10 terms (derivative):" ,xexpect(maxloc(error))

deallocate(x,y,xapprox,yapprox,xexpect,yexpect,error)


contains

  subroutine gaussj(n, A, B)
    integer ::n
    real:: A(n,n), B(n)
    integer ::i,j,pivot
    real ::temp(n,n), tempb, divnum !this is used for swapping
    integer ::loop

    i=1 !initialization of pivot row
    j=1 !initialization of pivot column
    do while(i<=n .and. j<=n)

      !step 1: if the aij==0, needs to swap the rows
      if (A(i,j)==0) then
        pivot=maxloc(abs(A(i,j:)), dim=1)

        if (A(i-1+pivot,j)==0) then ! special case: all row entries at j-th column is 0
          j=j+1 !go to the next column
          cycle ! go to the next iteration
        else
          !if all other rows entries are not zero, swap the rows
          temp(i,:) = A(i,:)
          A(i,:) = A(i-1+pivot,:)
          A(i-1+pivot,:) = temp(i,:)

          tempb = B(i)
          B(i) = B(i-1+pivot)
          B(i-1+pivot) = tempb

        end if
      end if

      !step2: divide row to obtain 1
      divnum=A(i,j)
      A(i,:) = A(i,:)/divnum
      B(i) = B(i)/divnum

      !step3: eliminate all other entries at j-th column --> 0
      do loop=1,n
        if(loop /= i) then
          B(loop)=B(loop)-A(loop,j)*B(i) !this one has to be first as you are using A(j,loop)
          A(loop,:)=A(loop,:)-A(loop,j)*A(i,:)
        end if
      end do


      i=i+1
      j=j+1
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


  subroutine printMatrix(row, col, M)
    integer :: row, col,i
    real :: M(row,col)

    do i=1,row
      write(*,*) M(i,:)
    end do
  end subroutine



  elemental function f(x)
    real,intent(in) ::x
    real ::f

    f=(1.0+10.0*x**2.0)**(-1.0)
  end function

  elemental function df(x)
    real,intent(in) ::x
    real ::df

    df=-20.0*x/(10.0*x**2.0+1.0)**2.0
  end function


  pure function populateArray(a,b,n)
    real,intent(in) :: a,b
    integer, intent(in) ::n
    integer ::i
    real ::populateArray(n)

    populateArray(1:n) = (/(a+(b-a)*i/real((n-1)),i=0,n-1)/)
  end function


  ! Chebyshev Term
  function ChebyshevT(x, n); intent(in) x, n
    real :: ChebyshevT, x
    integer :: n

    ChebyshevT = cos(n*acos(x))
  end function


  ! derivative of Chebyshev Term
  function  dChebyshevT(x, n); intent(in) x, n
    real :: dChebyshevT, x
    integer :: n

    dChebyshevT = n * sin( n * acos(x)) / sqrt( 1 - x**2 )
  end function


  ! Approximated point using sum of Chebyshev terms
  function ChebyshevTSums(c,xapprox,n,k)
    real,intent(in) ::c(n), xapprox(k)
    integer, intent(in) ::n,k

    integer ::i,j
    real::ChebyshevTSums(k), sums

    do i=1,k ! i-th number of approximating point
      sums=0
      do j=1,n !j-th number of term
        sums = sums + c(j) * ChebyshevT(xapprox(i), j-1) ! make sure j-1 !!!!!
      end do
      ChebyshevTSums(i)=sums
    end do
  end function


  !derivative of ChebyshevTSums
  function dChebyshevTSums(c,xapprox,n,k)
    real,intent(in) ::c(n), xapprox(k)
    integer, intent(in) ::n,k

    integer ::i,j
    real::dChebyshevTSums(k), sums

    do i=1,k ! i-th number of approximating point
      sums=0
      do j=1,n !j-th number of term
        sums = sums + c(j) * dChebyshevT(xapprox(i), j-1) ! make sure j-1 !!!!!
      end do
      dChebyshevTSums(i)=sums
    end do
  end function


  ! approximating value using ChebyshevT given x and y array sample points
  function fapprox(x, y, xapprox, n, k)
    real,intent(in) ::x(n), y(n), xapprox(k) !x, y are sample points
    integer, intent(in) ::n,k
    real ::fapprox(k)
    real ::A(n,n)
    real ::B(n)
    integer ::i,j !
    ! i: i-th sample point
    ! j: j-th term

    !Based on the sample points, compute the coefficients
    do i=1,n
      do j=1,n
        A(i,j) = ChebyshevT(x(i),j-1) ! make sure j-1 !!!!!
      end do
      B(i) = y(i)
    end do

    call gaussj(n,A,B)
    !call lsolve(n,A,B) ! this also works
    !now you got the coefficients

    fapprox(1:k) = ChebyshevTSums(B, xapprox, n, k)
  end function




  ! derivative approximation
  function dfapprox(x,y,xapprox,n,k)
    real,intent(in) ::x(n), y(n), xapprox(k)
    integer, intent(in) ::n,k
    real ::dfapprox(k)
    real ::A(n,n)
    real ::B(n)
    integer ::i,j

    do i=1,n
      do j=1,n
        A(i,j) = ChebyshevT(x(i),j-1) ! make sure j-1 !!!!!
      end do
      B(i) = y(i)
    end do

    call gaussj(n,A,B)
    !call lsolve(n,A,B) ! this also works
    !now you got the coefficients

    dfapprox(1:k) = dChebyshevTSums(B, xapprox, n, k)
  end function



end program

