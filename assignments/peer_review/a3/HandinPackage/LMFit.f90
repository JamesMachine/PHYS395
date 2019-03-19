!module containing code for bracketed minimum search by golden sections
!compile with:
!gfortran -c -fdefault-real-8 LMFit.f90 -llapack
!Last edit: 3/5 12:24pm

module LMFit
implicit none

! Interface for calling the f&df functions
! Forall block needs these to be pure
interface
	pure function f3Int(c, x)
		real, intent(in) :: c(:), x
		real :: f3Int
	end function f3Int
	
	pure function df3Int(c, x)
		real, intent(in) :: c(:), x
		real :: df3Int(size(c))
	end function df3Int
end interface

contains

!Use Levenberg  Marquardt algorithm to minimize chi2
function levenMarq(x, y, f, df, c0, lambda0, eps, maxLoops) result(c)
	procedure(f3Int) :: f
	procedure(df3Int) :: df
	real, intent(in) :: x(:), y(:), c0(:), lambda0, eps
	integer, intent(in) :: maxLoops
	integer i, loop, numPoints
	real, parameter :: lambdaFactor = 1.1

	real :: sumSq, sumSqNext, lambda, epsNext
	real :: J(size(x), size(c0)), JT(size(c0),size(x))
	real, dimension(size(x)) :: diff, diffNext
	real, dimension(size(c0)) :: c, cNext, JTxDiff
	real, dimension(size(c0), size(c0)) :: deltaM, JTxJ

	!initialize vars
	numPoints = size(x)
	loop = 0
	lambda = lambda0
	c = c0
	epsNext = eps+1 !dummy initial so we dont stop at first run
	forall (i=1:numPoints) diff(i) = y(i) - f(c, x(i)) 	!first sumSq use data
	sumSq = sum(diff**2.0)

	do while (epsNext >= eps)
		!Jacobians
		forall (i=1:size(x)) JT(:, i) = df(c, x(i))
		J = transpose(JT)
		JTxJ = matmul(JT, J)

		!Right and left sides of eq
		JTxDiff = matmul(JT, diff) !right side
		deltaM = JTxJ !deltaM is Matrix next to delta in the eq  
		forall (i=1:size(c0)) deltaM(i,i) = lambda * JTxJ(i,i)
		
		!Solve eq with lss, get next parameters
		cNext = c + lss(deltaM, JTxDiff)
		forall (i=1:size(x)) diffNext(i) = y(i) - f(cNext, x(i))
		sumSqNext = sum(diffNext**2.0)
		epsNext = abs(sumSqNext - sumSq)
		
		!accept next coeffs and down lambda if sumSq is smaller (better fit)
		if (sumSqNext < sumSq) then
			c = cNext
			diff = diffNext
			sumSq = sumSqNext
			lambda = lambda / lambdaFactor
		else
			lambda = lambda * lambdaFactor
		end if
		
		loop = loop + 1
		if (loop >= maxLoops) then !easy to accidentally inf loop here
			write(0,*) "Max runs on LM fitting algorithm reached"			
			exit
		end if
	end do
end function

! Use dgelss to solve A(x) = B
function lss(A, B)
	integer :: numPoints, n, LWork, rank, info
	real, intent(in) :: A(:,:), B(:)
	real rcond
	real, dimension(:,:), allocatable :: ACopy
	real, dimension(:), allocatable :: BCopy, sig, work, lss
	
	numPoints = size(A,1)
	n = size(A,2)
	LWork = 6*n
	allocate(ACopy(numPoints,n), BCopy(numPoints), sig(n), work(LWork), lss(n))
	ACopy = A
	BCopy = B

	call dgelss(numPoints, n, 1, ACopy, numPoints, BCopy, numPoints, sig, rcond, rank, work, LWork, info)	
	if (info /= 0) call abort

	lss = BCopy(1:size(lss))
	deallocate(ACopy, BCopy, sig, work)
end function

! basis is b_a(x)=cos(2pi*a*x)
pure function basis(a, x)
	integer, intent(in) :: a
	real, intent(in) :: x
	real :: basis
	real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510

	basis = cos(2 * pi * a * x)
end function

! non-linear model for fit f3 = e^(sum(c*b)) + const
pure function f3(c, x)
	real, intent(in) :: c(:), x
	real :: f3, b(size(c) - 1)
	integer :: a, n
	n = size(c) - 1
	
	forall (a=1:n) b(a) = basis(a-1 ,x)
	f3 = exp(sum(c(1:n) * b)) + c(n+1)
end function

! df3 = f3' wrt b
pure function df3(c, x)
	real, intent(in) :: c(:), x
	real :: df3(size(c)), b(size(c) -1), factor
	integer :: a, n
	n = size(c) - 1
	
	forall (a=1:n) b(a) = basis(a-1, x)
	df3(1:n) = exp(sum(c(1:n) * b)) * b
	df3(n+1) = 1.0
end function
end module
