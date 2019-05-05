!Assignment 1 - q2_2 - 100 terms

!Running Instruction for this specific file
!gfortran -fdefault-real-8 q2_2.f90 -llapack -o q2_2 && ./q2_2

program q2
implicit none

!constants
integer, parameter :: n=100 !choose number of samples
integer, parameter :: k=1000 !choose output points as approximation
!n: number of samples
!k: number of output points

!declaration
integer:: i
real,dimension(:),allocatable::x,y,xapprox,yapprox




allocate(x(n),y(n),xapprox(k),yapprox(k))
x=populateArray(-1.0,1.0,n)
xapprox=populateArray(-1.0,1.0,k)


open(1, file = './output/output_q2/q2_100terms_sample.dat', status='unknown')
do i=1,n
  y(i)=f(x(i))
  write(1,*)  x(i), y(i)
end do
close(1)

open(1, file = './output/output_q2/q2_100terms_approx.dat', status='unknown')
yapprox = fapprox(x,y,xapprox,n,k)
do i=1,k
  write(1,*) xapprox(i), yapprox(i)
end do
close(1)


deallocate(x,y,xapprox,yapprox)

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


  pure function populateArray(a,b,n)
    real,intent(in) :: a,b
    integer, intent(in) ::n
    integer ::i
    real ::populateArray(n)

    populateArray(1:n) = (/(a+(b-a)*i/real((n-1)),i=0,n-1)/)
  end function


  elemental function ChebyshevT(x, n)
    real ChebyshevT, x; integer n
    intent(in) x, n

    ChebyshevT = cos(n*acos(x))
  end function



  function ChebyshevTSums(c,xapprox,n,k)
    real,intent(in) ::c(n), xapprox(k)
    integer, intent(in) ::n,k

    integer ::i,j
    real::ChebyshevTSums(k), sums

    do i=1,k ! i-th number of approximating point
      sums=0
      do j=1,n !j-th number of term
        sums=sums+c(j)*cos((j-1)*acos(xapprox(i))) ! make sure j-1 !!!!!
        ChebyshevTSums(i)=sums
      end do
    end do


  end function


  ! approximating value using ChebyshevT given x and y array sample points
  function fapprox(x,y,xapprox,n,k)
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
        A(i,j)=ChebyshevT(x(i),j-1) ! make sure j-1 !!!!!
      end do
      B(i) = y(i)
    end do

    call gaussj(n,A,B)
    !call lsolve(n,A,B) ! this also works
    !now you got the coefficients

    fapprox(1:k) = ChebyshevTSums(B, xapprox, n, k)
  end function


end program

