! harmonic_shooter.f90 integrate Schrodingerequation from x=0 to infinity
! using the shooting method

! outputs some files:
! fort.1 is x, psi, normalized psi^2 for a given E. Created in dump
! fort.2 is the path taken to get to the even eigenvalue. Created in bisect
! fort.3 is the path taken to get to the odd eigenvalue. Created in bisect
! fort.4 is the list of eigenvalues, and boundary value violations
! compile with: gfortran -O3 -fdefault-real-8 harmonic_shooter.f90
program harmonic_shooter
implicit none

! number of iterations to find solution
integer, parameter :: iterations = 100

! size of steps to take
real, parameter :: dt = 0.01

! counters
integer i, j, k, narg, cutoff

! variables and things we want to save
real a, b, c, fa, fb, fc, E, Elim(2), psi0, phi0, roots(2), violation(2)

integer :: status = 0

character(len=8000) :: arg

real :: u(4), psi (1001, 3)

roots = -1.0
violation = -1.0

! find eigenvalue by bisection for even solution, then odd solution
do j = 1, 2
	call initE()
	if (fa*fb <= 0.0) then
		! bisect the interval
		do i = 1,iterations
			call bisect()
			if (fc == 0.0) exit
		end do
		roots(j) = c
		cutoff = truncate(size(psi,1))
		violation(j) = psi(cutoff, 2)
		call dump(cutoff)
	end if
end do

write (*,*) 'An energy eigenvalue or boundary violation of -1.0 means an eigenvalue was not found in the interval'
write (4,'(4g24.16)') roots(1)
write (4,'(4g24.16)') roots(2)

contains

! ======================================================================
! initialize energy variables from command argument
! calls: f
! ======================================================================
subroutine initE ()
	write (*,*) 'Interval is:'
	call get_command_argument(1, arg)
	write (*,*) trim(arg)
	read (arg,*,iostat=status) a

	call get_command_argument(2, arg)
	write (*,*) trim(arg)
	read (arg,*,iostat=status) b

	psi0 = 2.0-j; phi0 = j-1.0
	fa = f(a)
	fb = f(b)
end subroutine initE

! ======================================================================
! perform a bisection on an interval
! calls: f
! ======================================================================
subroutine bisect ()
	c = (a+b)/2.0; fc = f(c)

	write (j+1,*) c, fc

	if (fa*fc < 0.0) then; b = c; fb = fc; end if
	if (fc*fb < 0.0) then; a = c; fa = fc; end if
end subroutine bisect

! ======================================================================
! integrate schrodinger equation from x=0 to x=infinity
! called by: f
! calls: gl10, V
! ======================================================================
subroutine integrate(t)
	real t
	integer l, n

	! number of time steps needed
	n = floor(t/dt)

!	write (1,'(4g24.16)') u(1), u(2), u(2)*u(2), u(4)
	psi = 0.0
	psi(1,:) = [u(1), u(2), u(4)]
	do l = 1,n
		call gl10(u, dt); if (abs(u(2)) > 1000.0) exit
		psi(l+1,:) = [u(1), u(2), u(4)]
!		write (1,'(4g24.16)') u(1), u(2), u(2)*u(2), u(4)
	end do

	call gl10(u,t-n*dt)
end subroutine integrate

! ======================================================================
! function to find a root of
! calls: integrate
! called by: bisect
! ======================================================================
function f(x); intent(in) x
	real f, x
	E = x
	u = [0.0, psi0, phi0, 0.0]
	call integrate(10.0)
	f = u(2)
end function f

! ======================================================================
! potential energy for the problem
! called by: integrate
! ======================================================================
pure function V(x); intent(in) x
	real V, x

	V = 0.5*x*x
end function V

! ======================================================================
! truncate psi at the appropriate spot. Removes weird stuff far outside
! the region we care about.
! ======================================================================
function truncate(sz)
	integer i, truncate, sz, imin
	real minimum

	minimum = 1000.0
	imin = sz

	do i = sz,1,-1
		if (abs(psi(i,2)) < minimum .and. abs(psi(i,2)) > 0.0) then
		minimum = abs(psi(i,2))
		imin = i
		end if
	end do

	truncate = imin
end function truncate

! ======================================================================
! dump results into fort.1
! x, psi, psi^2/2N
! ======================================================================
subroutine dump(cutoff)
	integer m, cutoff
	real norm

	norm = psi (cutoff,3)

	do m = 1,cutoff
		write (1,'(4g24.16)') psi(m,1), psi(m,2), &
		psi(m,2)*psi(m,2)/(norm*2.0)
	end do

	write(1,*) ''
	write(1,*) ''
end subroutine

! ======================================================================
! evaluate derivatives
! called by: gl10
! calls: V
! ======================================================================
subroutine evalf(u, dudx)
!dpsi/dx=pi, d2psi/dx2=2v-E, dx/dx=1.0, dN/dx=psi2
	real u(4), dudx(4)

	associate( x => u(1), psi => u(2), phi => u(3), norm => u(4) )
	dudx(1) = 1.0
	dudx(2) = phi
	dudx(3) = 2*(V(x)-E)*psi
	dudx(4) = psi*psi
	end associate
end subroutine evalf

! ======================================================================
! 10th order implicit Gauss-Legendre integrator
! calls: evalf
! called by: integrate
! ======================================================================
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 4
        real y(n), g(n,s), dt; integer i, k

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
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do

        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

end program
