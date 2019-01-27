!Running Instruction
!gfortran -fdefault-real-8 test.f90 -o test && ./test

program test

real ::A(2,2), temp(2,2)

A(1,1) = 7.0
A(1,2) = 2.0
A(2,1) = 5.0
A(2,2) = 10.0

call printMatrix(2,2,A)
!temp(:,2) = A(:,1)
!A(:,1) = A(:,2)
!A(:,2) = temp(:,2)

A(:,1)=A(:,1)*2
call printMatrix(2,2,A)



contains

  subroutine printMatrix(row, col, M)
    integer :: row, col
    real :: M(col,row)

    do i=1,row
      write(*,*) M(:,i)
    end do
  end subroutine


end program
