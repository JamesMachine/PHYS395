!Assignment 1 - q1

!Running Instruction for this specific file
!gfortran -fdefault-real-8 q1.f90 -o q1 && ./q1 > output1 && cat output1


program a1code
implicit none

!declaration

real ::A(3,3), B(1,3)

integer:: i,n


!Question 1
write(*,*) 'Example problem1:'

A(1,:) = [0.0, 1.0, 3.0]
A(2,:) = [3.0, 300.0, 9.0]
A(3,:) = [4.0, 4.0, 1.0]

B(1,1) = 5.0
B(1,2) = 600.0
B(1,3) = 8.0

write(*,*) 'A'
call printMatrix(3,3,A)

write(*,*) 'b'
call printMatrix(3,1,B)

call gaussj(3,A,B)

write(*,*) 'After Gaussian Elimination:'

write(*,*) 'A'
call printMatrix(3,3,A)

write(*,*) 'b'
call printMatrix(3,1,B)


contains

  subroutine gaussj(n, A, B)
    integer ::n; real:: A(n,n), B(1,n)
    integer ::i,j,pivot
    real ::temp(n,n), tempb, divnum !this is used for swapping
    integer ::loop

    i=1 !initialization of pivot row
    j=1 !initialization of pivot column
    do while(i<=n .and. j<=n)

      !step 1: if the aij==0, needs to swap the rows
      if (A(j,i)==0) then
        pivot=maxloc(abs(A(j,i:)), dim=1)

        if (A(j,i-1+pivot)==0) then ! special case: all row entries at j-th column is 0
          j=j+1 !go to the next column
          cycle ! go to the next iteration
        else
          !if all other rows entries are not zero, swap the rows
          temp(:,i) = A(:,i)
          A(:,i) = A(:,i-1+pivot)
          A(:,i-1+pivot) = temp(:,i)

          tempb = B(1,i)
          B(1,i) = B(1,i-1+pivot)
          B(1,i-1+pivot) = tempb

        end if
      end if

      !step2: divide row to obtain 1
      divnum=A(j,i)
      A(:,i) = A(:,i)/divnum
      B(1,i) = B(1,i)/divnum

      !step3: eliminate all other entries at j-th column --> 0
      do loop=1,n
        if(loop /= i) then
          B(1,loop)=B(1,loop)-A(j,loop)*B(1,i) !this one has to be first as you are using A(j,loop)
          A(:,loop)=A(:,loop)-A(j,loop)*A(:,i)
        end if
      end do


      i=i+1
      j=j+1
    end do

  end subroutine


  subroutine printMatrix(row, col, M)
    integer :: row, col,i
    real :: M(col,row)

    do i=1,row
      write(*,*) M(:,i)
    end do
  end subroutine


end program

