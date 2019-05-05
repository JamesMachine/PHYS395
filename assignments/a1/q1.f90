!Assignment 1 - q1 -Solving Gauss-Jordan Elimination

!Running Instruction for this specific file
!gfortran -fdefault-real-8 q1.f90 -llapack -o q1 && ./q1 > output1 && cat output1


program q1
implicit none

!constants
integer,parameter :: n = 3

!declaration

real ::A(n,n), B(n)

integer:: i


!Question 1
write(*,*) 'Example problem1:'

A(1,:) = [0.0, 1.0, 3.0]
A(2,:) = [3.0, 300.0, 9.0]
A(3,:) = [4.0, 4.0, 1.0]

B(:) = [5.0, 600.0, 8.0]

write(*,*) 'A'
call printMatrix(n, n, A)

write(*,*) 'b'
call printMatrix(n,1,B)

! solve either way:
call gaussj(n, A, B)
!call lsolve(n, A, B)

write(*,*) 'After Gaussian Elimination:'

write(*,*) 'A'
call printMatrix(n, n, A)

write(*,*) 'b'
call printMatrix(n, 1, B)


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


end program

