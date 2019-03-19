! Assign1.f90  -  Solve linear square matrixes.
! Use text file in Row-formatting, with solutions transposed as a row. (Add another row instead of column)
! Text file must be named "matrix.txt"
! compile with: gfortran -O3 -fdefault-real-8 Assign.f90 -llapack -o Assign1
!
! NOTE: This only works on matrixes of form (n,n+1) due to the solution column being added as a row.


program gaussjordan
implicit none

integer k, j, i, l
! Must manually change matrixes to desired form. Give C(n,n+1) for data.
real A(3,3), B(3), C(3,4), D(3,3)

! This is the dimensions of the linear system matrix
k = 3

l = k+1

! Opens up text file with linear system. Name the file to "matrix.txt"
  open(12, FILE="matrix.txt")


! Places text data into a matrix format.
read(12,*) C

! Creates the solution and linear system matrixes.
  do i=1,k

    B(i) = C(i,l)

      do j = 1,k

        A(i,j) = C(i,j)

      end do

  end do

write(*,*) A

! Undoes the tranpose from fortran reading the data.
A = transpose(A)

write(*,*)

! Displays the linear system.
write (*,*) B

write (*,*)



! Calls the lapack subroutine to solve the matrix.
! Alter the 3 to your dimensions needed.
call lsolve(3, A, B)

! Shows the solutions under their respective columns.
write (*,*) B

close(12)

contains

! solve A.x = B using LAPACK xGESV
! A gets destroyed, answer is returned in B
subroutine lsolve(n, A, B)
	integer n, pivot(n), status; real A(n,n), B(n)

	! initializes status to all clear
	status = 0

	! calls appropriate version of xGESV
	select case (kind(A))
		case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
		case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
		case default; call abort
	end select

	! abort at first sign of trouble
	if (status /= 0) stop "singular matrix in lsolve()"
end subroutine


end

