program Q2
implicit none

! number of iterations for bisection
integer, parameter :: iterations = 1000

! value which we would like to aim for as psi -> inf
real, parameter :: tt = 0.0

! direction to shoot. +1.0 for postive, -1.0 for negative directions
real, parameter :: dir = 1.0

! initial guess for the wave function in terms of x, psi, d(psi)/dx
real, dimension(3) :: init=[0.0, 1.0, 0.0]

! declare spatial step size (suggested from collegue)
real, parameter :: dx = 6.2831853071795864769252867665590057684Q0/2**9

! declare which potential to use. 0 for quadratic, 1 for quartic
integer, parameter :: option = 0
real, parameter :: par = 1.0

integer h
real j

write(*,*) "The eigenvalues for the quadratic potential are:"
! use bisection method to compute the eigenvalues given the quadratic potential
do h = 0,9
	write(*,*) bisection(h*1.0, h*1.0 + 1.0, par)
end do

contains

! check the values of psi at both pos and neg inf
function check(E, par); intent(in) E, par
    real par, E, check
    real y(3), y_neg(3)
    integer l

    y = init
    y_neg = -1*init

    do l=1,2**10
        call gl10(y, dx, 1.0, E) ! pos direction
        call gl10(y_neg, dx, -1.0, E) ! neg direction
    end do

    ! check the boundary values of both pos and neg
    check = (y(2) + par*y_neg(2)) - tt
end function

! quadratic potential
function V(x); intent(in) x
    real V, x

    V = 0.5*(x*x)

end function

function bisection(a0, b0, par); intent(in) a0, b0, par
    real a0, b0, par, bisection
    real a, b, c, fa, fb, fc
    real, parameter :: eps=1e-14

    ! initial interval
    a = a0; fa = check(a, par)
    b = b0; fb = check(b, par)

    if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."

    ! bisect the interval
    do while (abs(b-a) > eps*1.0)
        c = (a+b)/2.0; fc = check(c, par); if (fc == 0.0) exit

        ! Ridder's variation on the basic method
        c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = check(c, par)

        if (fc == 0.0) exit
        if (fa*fc < 0.0) then; b = c; fb = fc; end if
        if (fc*fb < 0.0) then; a = c; fa = fc; end if
    end do

    bisection = c
end function

! evaluate derivatives
subroutine evalf(u, dudt, dir, E)
        real u(3), dudt(3), E, dir

        associate( x => u(1), psi => u(2), phi => u(3) )
        dudt(1)=dir*1; dudt(2)=phi; dudt(3)=2.0*(V(x)-E)*psi
        end associate
end subroutine

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt, dir, E)
        integer, parameter :: s = 5, n = 3
        real y(n), g(n,s), dt, E, dir; integer i, k

        ! Butcher tableau for 10th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1/), (/s,s/))
        real, parameter ::   b(s) = [ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1]

        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i), dir, E)
                end do
        end do

        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine

end program
