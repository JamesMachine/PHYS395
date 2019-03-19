!compile with:
!gfortran -c -fdefault-real-8 tester.f90
!link and run with:
!gfortran tester.o gaussjordan.o polyapprox.o -o tester && ./tester > output && cat output

program tester
use gaussjordan
use polyapprox
implicit none

!All variables
integer i, n, m
real A(3,3), B(3)
real, dimension(:), allocatable :: x, y, xapprox, yapprox, error
!Problem 1:
n = 3
A(:,1) = [0.0, 2.0, 1.0]
A(:,2) = [1.0, -2.0, -3.0]
A(:,3) = [-1.0, 1.0, 2.0]
B = [-8.0, 0.0, 3.0]

write (*,*) "Problem 1: Solve the linear system A|B:"
do i=1,3
	write (*,*) A(:,i), "|", B(i)
end do
call gaussj(n, A, B)
write (*,*) "System solution:"
write (*,*) B

!Problem 2 (Sample size 10):
n = 10
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "freal_10.csv", status = "unknown")
x = uniform(n)
do i=1,n
	y(i) = f(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "fapprox_10.csv", status = "unknown")
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

!Problem 2 (Sample size 100):
n = 100
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "freal_100.csv", status = "unknown")
x = uniform(n)
do i=1,n
	y(i) = f(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "fapprox_100.csv", status = "unknown")
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

write (*,*) "Problem 2: Approximations to f(x) = 1/(1+10*x^2) printed to csv files"
write (*,*) "Plots will be printed to pdf after"

!Problem 3 (Sample size 100):
n = 100
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "greal_100.csv", status = "unknown")
x = uniform(n)
do i=1,n
	y(i) = g(x(i))
	write(1,*) x(i), y(i)
  y(i) = f(x(i))
end do
close(1)

open(1, file = "gapprox_100.csv", status = "unknown")
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,1)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

write (*,*) "Problem 3: Approximations to g(x) = f'(x) printed to csv files"
write (*,*) "Plots will be printed to pdf after"

!Problem 4 (f(x), Sample size 100):
allocate (x(m), y(m), xapprox(m), yapprox(m), error(m))
x = uniform(m)
do i=1,m
	y(i) = f(x(i))
end do
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
error(:) = abs(yapprox(:) - y(:))

write (*,*) "Problem 4a: Maximal error for f(x) approximation is ", maxval(error), " at x = ", x(maxloc(error))
deallocate (x, y, xapprox, yapprox, error)

!Problem 4 (g(x), Sample size 100):
allocate (x(m), y(m), xapprox(m), yapprox(m), error(m))
x = uniform(m)
do i=1,m
	y(i) = g(x(i))
end do
xapprox = uniform(m)
yapprox = fapprox(x,y,xapprox,n,m,1)
error(:) = abs(yapprox(:) - y(:))

write (*,*) "Problem 4b: Maximal error for g(x) approximation is ", maxval(error), " at x = ", x(maxloc(error))
deallocate (x, y, xapprox, yapprox, error)

!Problem 5
!A copy&paste of problems 2-4 but with chebyGrid(n) replacing all instances of uniform(n)
write (*,*) "Problem 5: Repeat of problems 2-4 with a non-uniform grid"
!Problem 5-2 (Sample size 10):
n = 10
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "freal_10_P5.csv", status = "unknown")
x = chebyGrid(n)
do i=1,n
	y(i) = f(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "fapprox_10_P5.csv", status = "unknown")
xapprox = chebyGrid(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

!Problem 5-2 (Sample size 100):
n = 100
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "freal_100_P5.csv", status = "unknown")
x = chebyGrid(n)
do i=1,n
	y(i) = f(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "fapprox_100_P5.csv", status = "unknown")
xapprox = chebyGrid(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

write (*,*) "Problem 5-2: Approximations to f(x) = 1/(1+10*x^2) printed to csv files"
write (*,*) "Plots will be printed to pdf after"

!Problem 5-3 (Sample size 100):
n = 100
m = 1000
allocate (x(n), y(n), xapprox(m), yapprox(m))
open(1, file = "greal_100_P5.csv", status = "unknown")
x = chebyGrid(n)
do i=1,n
	y(i) = g(x(i))
	write(1,*) x(i), y(i)
end do
close(1)

open(1, file = "gapprox_100_P5.csv", status = "unknown")
xapprox = chebyGrid(m)
yapprox = fapprox(x,y,xapprox,n,m,1)
do i=i,m
	write(1,*) xapprox(i), yapprox(i)
end do
close(1)
deallocate (x, y, xapprox, yapprox)

write (*,*) "Problem 5-3: Approximations to g(x) = f'(x) printed to csv files"
write (*,*) "Plots will be printed to pdf after"

!Problem 5-4a (f(x), Sample size 100):
allocate (x(m), y(m), xapprox(m), yapprox(m), error(m))
x = chebyGrid(m)
do i=1,m
	y(i) = f(x(i))
end do
xapprox = chebyGrid(m)
yapprox = fapprox(x,y,xapprox,n,m,0)
error(:) = abs(yapprox(:) - y(:))

write (*,*) "Problem 5-4a: Maximal error for f(x) approximation is ", maxval(error), " at x = ", x(maxloc(error))
deallocate (x, y, xapprox, yapprox, error)

!Problem 5-4b (g(x), Sample size 100):
allocate (x(m), y(m), xapprox(m), yapprox(m), error(m))
x = chebyGrid(m)
do i=1,m
	y(i) = g(x(i))
end do
xapprox = chebyGrid(m)
yapprox = fapprox(x,y,xapprox,n,m,1)
error(:) = abs(yapprox(:) - y(:))

write (*,*) "Problem 5-4b: Maximal error for g(x) approximation is ", maxval(error), " at x = ", x(maxloc(error))
deallocate (x, y, xapprox, yapprox, error)

end program
