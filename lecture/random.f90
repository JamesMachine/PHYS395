!random number generator
!gfortran -fdefault-real-8 random.f90 -o random && ./random

program random
  implicit none

  integer, parameter::n = 100
  real, parameter::pi = 3.14159265358979323846264338327950288419716939937510


  integer ::i
  real ::A(n), U(n), V(n), X(n)

  call random_seed()
  call random_number(U)
  call random_number(V)
  call random_number(X)


  A = sqrt(-2.0*log(U))*cos(2*pi*V) !fortran perform array operation automatically

  do i = 1,n
    write(*,*) X(i), (0.5*X(i)+1.0)+A(i)
  end do

end program

!gnuplot
!plot 'output' u 0:1
!plot 'output' -u 0:1, x/100
!plot 'output' u 1:($1-$0)/10000,erf(x)
!sort -k1g output > output_sorted


!
!?fit

!Box-muller transform


