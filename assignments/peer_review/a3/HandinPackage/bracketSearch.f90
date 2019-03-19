!module containing code for bracketed minimum search by golden sections
!compile with:
!gfortran -c -fdefault-real-8 bracketSearch.f90 -llapack
!Last edit: 3/5 12:34pm

module bracketSearch
implicit none

! Interface for calling the f&df functions
interface
    function f2Int(x)
        real, intent(in) :: x
        real :: f2Int
    end function f2Int
end interface

contains

function bracketMin(f, a0, b0, eps)
	procedure(f2Int) :: f
	real a0, b0, eps, bracketMin, a, b, c, d; intent(in) :: a0, b0, eps
	real, parameter :: phi = 0.5 + sqrt(5.0)/2.0
!	write(0,*) "phi=", phi
	
	!starting points
	a = a0
	b = b0
	c = b - (b-a)/phi
	d = a + (b-a)/phi
	
	!golden-section search, from while (eps < abs(c-d))
	do while(abs(c-d) > eps)
		if (f(c) < f(d)) then
			b = d
		else
			a = c
		end if
		
		c = b - (b-a)/phi
		d = a + (b-a)/phi
	end do
	
	bracketMin = (b+a)/2.0
end function

pure function f2(x)
	real f2, x; intent(in) x
	
	f2 = (x**2.0 - 1.0)**2.0 + x
end function
end module
