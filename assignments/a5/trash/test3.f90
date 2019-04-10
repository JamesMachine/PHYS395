program test3; implicit none

  integer ::i
  integer,parameter::n=150
  real::pi = 3.14159
  real,parameter::ell=1.0
  real::theta(n),x(n)

  real::u(3)
  u = [1.0,1.0,1.0]
  call f(u)
  write(*,*) u



contains
  subroutine f(y)
    real::y(1)
    write(*,*) y
  end subroutine

end program
