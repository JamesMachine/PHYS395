program question1

implicit none

real :: x
real, dimension(9) :: x_array
integer :: i, j

! initial guesses
x_array = (/-10.0, -5.0, -1.0, -0.5, 0.0, 0.5, 1.0, 5.0, 10.0/)

do j = 1,7
    x = x_array(j)
    print *, "Initial guess is:", x
    ! Newton's iteration
    do i = 1,16
        x = x - f(x)/df(x)
       ! write (*,*) x, f(x)
    end do
    write (*,*) "Final answer:", x
end do

contains

! function to find a root of...
pure function f(x); 
    real, intent(in) ::  x
    real :: f
    f = x*(x*x - 1) + 0.25
end function

! derivative of a function to find a root of
pure function df(x); 
    real, intent(in) :: x
    real :: df
    df = 3*x*x - 1
end function

end program question1