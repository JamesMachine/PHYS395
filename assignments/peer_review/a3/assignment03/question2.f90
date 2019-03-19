program question2
implicit none

! define the initial triplet
real :: ax_init
real :: bx_init
real :: cx_init

! declare tolerance variable
real :: epsilon

! declare variable that holds the argmin of the function
real :: xmin
! declare variable that holds the minimum value of the function
real :: minval

! initial points bracket the left
ax_init = -10.0
bx_init = 0.0
cx_init = 10.0

epsilon = 0.0001 ! define the tolerance

minval = gss(ax_init, bx_init, cx_init, epsilon, xmin)

print *, "Consider the window defined by [", ax_init, bx_init, cx_init, "]"
print *, "Minimum value:", minval
print *, "Argmin value:", xmin
print *, "Tolerance:", epsilon

contains

! Define function to minimize
pure function f(x)
    real, intent(in) :: x
    real :: f
    f = x*(x*x*x-2*x+1)+1
end function f

! function that performs golden section search
function gss(ax, bx, cx, tol, xmin)
    real, intent(in) :: ax, bx, cx
    real, intent(in) :: tol ! defines the tolerance
    real, intent(out) :: xmin
    real :: gss

    ! declare variable to hold value of the fractional part of the golden ratio
    real, parameter :: gr=0.61803398874989485
    ! define a coefficient C that is one minus golden ratio
    real, parameter :: C=1-gr
    real :: f1, f2, x0, x1, x2, x3

    ! track the endpoints of the interval
    x0=ax
    x3=cx

    ! determine the midpoint
    if (abs(cx-bx) > abs(bx-ax)) then
        x1=bx
        x2 = bx+C*(cx-bx)
    else
        x2=bx
        x1=bx-C*(bx-ax)
    end if

    ! initial function evaluations
    f1=f(x1); f2=f(x2)

    do
        ! exit criterion
        if (abs(x3-x0) <= tol*(abs(x1) + abs(x2))) exit

        ! one possible outcome
        if (f2 < f1) then
            ! shift the bounds accordingly 
            call shft3(x0, x1, x2, gr*x2+C*x3)
            ! function evaluation
            call shft2(f1,f2,f(x2))
        else ! another possible outcome
            ! shift the bounds accordingly
            call shft3(x3, x2, x1, gr*x1+C*x0)
            ! function evaluation
            call shft2(f2, f1, f(x1))
        end if
    end do

    if (f1 > f2) then
        gss=f1
        xmin=x1
    else
        gss=f2
        xmin=x2
    end if
end function gss

! define subroutines that shift variables around
subroutine shft2(a, b, c)
    real, intent(out) :: a
    real, intent(inout) :: b
    real, intent(in) :: c
    a=b; b=c
end subroutine shft2

subroutine shft3(a, b, c, d)
    real, intent(out) :: a
    real, intent(inout) :: b, c
    real, intent(in) :: d
    a=b; b=c; c=d
end subroutine shft3

end program question2