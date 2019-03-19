!10 term zeros grid

program tentermzero
implicit none

integer  n, i, j
real F(10), G(10), x(10), T(10,10), U(10,10), C(10), D(10)

n = 10

call Initialize(n,x,T,U,F,G)

call gaussj(n, T, F)


do i=1,n
	do j = 1,n
		C(i) = C(i) + F(j)*ChebyshevT(x(i),j-1)
		D(i) = D(i) + F(j)*DerChebyshevT(x(i),j-1)
	end do
end do

do i=1,n
	write (*,*) x(i), (1.0/(1.0 + 10.0*x(i)*x(i)))
end do

write (*,*) ""
write (*,*) ""

do i=1,n
	write (*,*) x(i),  C(i)
end do

write (*,*) ""
write (*,*) ""

do i=1,n
	write (*,*) x(i),  D(i)
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
!Initialize all the needed matrices
subroutine Initialize(n,x,T,U,F,G)
	integer n, i, j
	real x(n), T(n,n), F(n), U(n,n), G(n)

	do i = 1,n
		x(i) = cos((i-0.5)*(335.0/113.0)/n)
	end do

	forall (i=1:n) F(i) = (1.0/(1.0 + 10.0*x(i)*x(i)))
	forall (i=1:n) G(i) = (1.0/(1.0 + 10.0*x(i)*x(i)))

	do i=1,n
		do j=1,n
			T(i,j) = ChebyshevT(x(i),j-1)
		end do
	end do

	do i=1,n
		do j=1,n
			U(i,j) = DerChebyshevT(x(i),j)
		end do
	end do

end subroutine


! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n, i, j; real A(n,n), B(n), M(n,n+1)

!Make matrix M which is just [A|B]
	M(:,1:n) = A(:,:)
	M(:,n+1) = B(:)

	!Triangularize the A part of the matrix
	do j = 1,n-1
		do i = j+1,n
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
	end do

	!Diagonalize and normalize the A part of the matrix
	do j = n,2,-1
		do i = j-1,1,-1
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
		M(j,:) = M(j,:)/M(j,j)
	end do


	!The last column is now the solution
	B(:) = M(:,n+1)

end subroutine


! Chebyshev polynomial T_n(x)
elemental function ChebyshevT(x, n)
	real ChebyshevT, x; integer n
	intent(in) x, n

	ChebyshevT = cos(n*acos(x))
end function

!Derivative of Chebyshev polynomial
elemental function DerChebyshevT(x, n)
	real DerChebyshevT, x; integer n
	intent(in) x, n

	DerChebyshevT = n*sin(n*acos(x))/sqrt(1.0-x*x)
end function

end
