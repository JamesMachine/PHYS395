!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o newtonMethod.o bracketSearch.o && ./tester > output && cat output

program tester
use newtonMethod
use bracketSearch
use LMFit
implicit none

!Variable declarations
real min1X, min2X, chi2
real, dimension(:), allocatable :: x, y, yFit, fitResults
integer i, iost, numPoints
integer, parameter :: n = 3, maxPoints = 1000000
real xFile(maxPoints), yFile(maxPoints)

!Format labels, x = space, a = characters, nfw.d = n (f)loats (w)idth (d)ecimals
1 format(x, a, 1f20.16)
2 format(/, a)
3 format(x, a, 1f8.5, a, 1f8.5, a)
4 format(a5, i1, a, 1f11.8)
5 format(a9, 1f11.8)
6 format(a9, 1f11.5)

!Problem 1:
write (*,*) "Problem 1: Use Newton's method to solve x^3 -x +1/4 = 0 to double precision (eps = 1.0e-16)"
write (*,1) "Root x1 = ", newton(f1, df1, x0=-1.0, eps=1.0e-16)
write (*,1) "Root x2 = ", newton(f1, df1, x0=0.0, eps=1.0e-16)
write (*,1) "Root x3 = ", newton(f1, df1, x0=1.0, eps=1.0e-16)

!Problem 2&3:
write (*,2) "Problem 2&3: Use bracketed search with phi to minimize (x^2-1)^2 -x"
min1X = bracketMin(f2, a0=-1.5, b0=-0.5, eps=1.0e-6)
min2X = bracketMin(f2, a0=0.5, b0=1.5, eps=1.0e-6)
write (*,3) "Minimum (x1, y1) = (", min1X, ", ", f2(min1X), ")"
write (*,3) "Minimum (x2, y2) = (", min2X, ", ", f2(min2X), ")"

!Problem 4&5:
write (*,2) "Problem 4&5: Use Levenberg-Marquardt method to fit non-linear model"
! read file "data.dat", limit maxPoints
open(unit=1, file="data.dat", status="old", action="read")
do i = 1, maxPoints
	read(1, *, iostat=iost) xFile(i), yFile(i)
	if(iost < 0) exit
end do
close(1)
! check actual data extent
numPoints = i-1
if (numPoints == maxPoints) write(0,*) "Read data extent was truncated, recompile with larger limit"
! allocate data arrays
allocate(x(numPoints), y(numPoints), yFit(numPoints), fitResults(n+2))
x = xFile(1:i)
y = yFile(1:i)

! actual fitting
fitResults = 0.0
fitResults = levenMarq(x, y, f3, df3, c0=fitResults, lambda0=0.7, eps=0.1, maxLoops=100)
forall (i=1:size(x)) yFit(i) = f3(fitResults, x(i))
chi2 = sum((y-yFit)**2.0) !sigma = 1

! output results
do i = 1, n+1
	write(*,4) "c", i-1, " = ", fitResults(i)
end do
write(*,5) "Const = ", fitResults(n+2)
write(*,6) "Chi2 = ", chi2

! write datapoints for plotting
open(unit=2, file="LMFit.csv", status="replace", action="write")
write(2,*) "x, y, yFit"
do i = 1, numPoints
	write(2,*) x(i), ",", y(i), ",", yFit(i)
end do
close(2)

deallocate(x, y, yfit, fitResults)
end program

