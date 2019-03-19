!module containing code for polynomial approximation
!compile with:
!gfortran -c -fdefault-real-8 polyapprox.f90

module polyapprox
use gaussjordan
implicit none

contains

! uniformly spaced grid of n points covering interval [-1,1]
pure function uniform(n)
	integer, intent(in) :: n
	real uniform(n)
	integer i
	forall (i=1:n) uniform(i) = (2*i-n-1.0)/(n-1.0)
end function

! grid of n points given by zeroes of T_n on [-1,1]
pure function chebyGrid(n)
	integer, intent(in) :: n
	real chebyGrid(n)
	integer i
	real, parameter :: pi = 3.1415926535897932384
	forall (i=1:n) chebyGrid(i) = cos((pi/n)*(i-0.5))
end function


! function to be approximated f(x)
elemental function f(x)
	real f, x; intent(in) x
	f = 1.0/(1.0 + 10.0*x*x)
end function

! derivative of f(x)
elemental function g(x)
	real g, x; intent(in) x
	g = (-20.0*x)/((1.0+10.0*x*x)**2.0)
end function

! approximate f(x) or g(x) using chebyshev basis
function fapprox(x, y, xapprox, n, m, deriv)
	real, intent(in) :: x(n), y(n), xapprox(m)
	integer, intent(in) :: n, m, deriv
	integer i, j
	real fapprox(m), A(n,n), B(1,n)

	!populate A using basis
	do i=1,n
		do j=1,n
			A(j,i)=basis(j-1, x(i))
		end do
		B(1,i) = y(i)
	end do

	!solve with gauss jordan
	call gaussj(n,A,B)
	fapprox(1:m) = basisSums(B(1,:), xapprox, n, m, deriv)
end function

!basis with chebyshev polynomials
elemental function basis(n, x)
	integer, intent(in) :: n
	real, intent(in) :: x
	real basis
	basis = cos(n*acos(x))
end function

!summation with basis polynomial, deriv toggles derivative
function basisSums(b, xapprox, n, m, deriv)
	real, intent(in) :: b(n), xapprox(m)
	integer, intent(in) :: n, m, deriv
	integer i,j
	real basisSums(m), sums, numerator, denominator

	do i=1,m
	  sums = 0
		!term number
		do j=1,n
			if (deriv == 0) then
				sums = sums + b(j)*cos((j-1)*acos(xapprox(i)))
			else
				numerator = (j-1)*sin((j-1)*acos(xapprox(i)))
				denominator = sqrt(1.0-xapprox(i)**2.0)
				sums = sums + b(j)*numerator/denominator
			end if
		end do
		basisSums(i)= sums
	end do
end function
end module
