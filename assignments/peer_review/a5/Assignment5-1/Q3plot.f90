program soliton
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

real a, b, c, fa, fb, fc
integer i, varia

real energy(10,1)

real, dimension(10) :: lower
real, dimension(10) :: upper

integer h

! bracket the energies of quartic potential
lower=[0.38,1.47,2.92,4.50,6.40,8.32,10.48,12.60,14.90,17.41]
upper=[0.54,1.51,3.09,4.70,6.50,8.55,10.60,12.80,15.10,17.58]

! use biesction method to compute the eigenvalues given the quartic potential
do h=1,10
    ! initial guess for even potential
	if (mod((h-1), 2)==0) then
        init=[0.0, 1.0, 0.0]
    ! initial guess for odd potential
	else
        init=[0.0, 0.0, 1.0]
    end if
    energy(h,1) = bisect(lower(h), upper(h), par)
end do

do varia = 1,10
	if (mod((varia-1), 2)==0) then
		fa = integrate(0.0, 3.5, 1e-3, energy(varia,1), 1)
	else
		fa = integrate(1.0, 3.5, 1e-3, energy(varia,1), 0)
	end if
end do

contains

! potential
pure function V(x); intent(in) x
	real V, x

	V = 0.25*(x**4)
end function

! evaluate derivatives
subroutine evalf(u, dudt, num)
        real u(3), dudt(3), num

        associate( x => u(1), psi => u(2), phi => u(3) )
        dudt(1) = 1.0; dudt(2) = u(3); dudt(3) = 2.0*(V(x) - (energy(num,1)))*u(2)
        end associate
end subroutine

! evaluate derivatives
subroutine evalf1(u, dudt, dir, E)
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
                        call evalf1(y + g(:,i)*dt, g(:,i), dir, E)
                end do
        end do

        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine

! 10th order implicit Gauss-Legendre integrator
subroutine gl101(y, dt, num)
        integer, parameter :: s = 5, n = 3
        real y(n), g(n,s), dt, num; integer i, k

        ! Butcher tableau for 8th order Gauss-Legendre method
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
                        call evalf(y + g(:,i)*dt, g(:,i), num)
                end do
        end do

        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl101

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

function bisect(a0, b0, par); intent(in) a0, b0, par
    real a0, b0, par, bisect
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

    bisect = c
end function

! integrate equations of motion for a given amount of time
function integrate(phi, t, dt, num, choice); intent(in) phi, t, dt, num, choice
	real phi, t, dt, integrate, num
	real u(3), E0; integer i, n, choice

	! start from x=0
	if (choice==1) then
		u = [0.0, 1.0, phi]
	else
		u = [0.0, 0.0, phi]
	end if

	! number of time steps needed
	n = floor(t/dt)

	do i = 1,n
		call gl101(u, dt, num); if (abs(u(1)) > 100.0) exit
		if (mod(i, 100)==0) write (*,'(4g24.16)') u(1), u(2), u(3)
	end do

	call gl101(u,t-n*dt, num)

	! return state at time t
	integrate = u(1)
	write (*,*); write (*,*)
end function

end program