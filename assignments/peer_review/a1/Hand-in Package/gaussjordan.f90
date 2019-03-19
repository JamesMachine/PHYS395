!module containing code for Gauss-Jordan elimination
!compile with:
!gfortran -c -fdefault-real-8 gaussjordan.f90

module gaussjordan
implicit none
contains

! solve A.x = B
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n; real A(n,n), B(1,n)
	integer i, j, k, pivot
	real tempA(n,n), tempB
 
	! get matrix into upper-triangular form
	i = 1
	j = 1
	do while (i <= n .and. j <= n)
    ! find pivot and swap if 0 on diagonal 
		if (A(j,i) == 0) then
			pivot = maxloc(abs(A(j,i:)), dim=1) + i-1		
			! edge case for entire col is 0
			if (A(j,pivot) == 0) then
				j = j+1
				cycle
			end if
			tempA(:,i) = A(:,i)
			A(:,i) = A(:,pivot)
			A(:,pivot) = tempA(:,i)
			tempB = B(1,i)
			B(1,i) = B(1,pivot)
			B(1,pivot) = tempB
		end if

		! get row reduced form
		B(1,i) = B(1,i)/A(j,i)
		A(:,i) = A(:,i)/A(j,i)
		! "back sub" to echelon form
		do k = 1,n
			if(k /= i) then
				B(1,k) = B(1,k) - A(j,k)*B(1,i)
				A(:,k) = A(:,k) - A(j,k)*A(:,i)
			end if
		end do
		i=i+1
		j=j+1
	end do
end subroutine
end module
