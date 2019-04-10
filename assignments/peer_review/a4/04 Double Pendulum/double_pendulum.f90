! double_pendulum.f90 - model a chaotic double pendulum system, where each arm has
! equal mass and length
! compile with: gfortran -O3 -fdefault-real-8 double_pendulum.f90

program double_pendulum2; implicit none

! physical constants
real, parameter :: pi = 4.D0*DATAN(1.D0)
real, parameter :: g = 9.8

! parameters of the system
real, parameter :: m = 1.0
real, parameter :: l = 1.0

! characteristic timestep
real, parameter :: dt = sqrt(l/g)/100.0

real :: th1(2) = (/ pi/3.0, 0.0 /)
real :: th2(2) = (/ -pi/3.0, 0.0 /)
real :: data(4,1,1)
integer i



data(:,1,1) = integrate(th1(1), th2(1), 100.0*sqrt(l/g), dt)

contains

! integrate and output some information about the system
function integrate(th1, th2, t, dt)
	real th1, th2, t, dt, integrate(4)
	real u(4), E0, cartesian (2,2); integer n
	
	! start from a given position at rest
	u = [th1, th2, 0.0, 0.0]; E0 = Etot(u)
	
	! number of time steps needed
	n = floor(t/dt)
	
	do i = 1,n
		call gl10(u, dt)
		
		cartesian = polar2cart(u(1),u(2))
		write (1,'(4g24.16)') i*dt, cartesian(2,:), Etot(u)/E0 - 1.0
		call animation(u(1),u(2))
		
	end do

	call gl10(u,t-n*dt)
	
	! return state at time t
	integrate = u
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATIONS USED TO CHANGE THE STATE OF THE SYSTEM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! evaluate derivatives
subroutine evalf(u, dudt)
        real u(4), dudt(4)
        
        associate( th1 => u(1), th2 => u(2), thdot1 => u(3), thdot2 => u(4) )
        dudt(1) = thdot1
        dudt(2) = thdot2
        dudt(3) = (-g*3*m*sin(th1) - m*g*sin(th1-2*th2) - 2*m*l*sin(th1-th2)*(thdot2**2+thdot1**2*cos(th1-th2))) &
			/(l*m*(3-cos(2*th1-2*th2)))
        dudt(4) = (2*sin(th1-th2)*(thdot1**2*l*2*m + g*2*cos(th1) + thdot2**2*l*m*cos(th1-th2))) &
			/(l*m*(3 - cos(2*th1 - 2*th2)))
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 4
        real y(n), g(n,s), dt; integer i, k
        
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
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATIONS USED TO DETERMINE THE ENERGY OF THE SYSTEM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! potential for double pendulum
pure function Epot(th1,th2); intent(in) th1, th2
	real Epot, th1, th2
	
	Epot = -g*l*m*(2*cos(th1)+cos(th2))
end function Epot

! total energy for the double pendulum
pure function Ekin(u); intent(in) u
        real u(4), Ekin
        
        associate( th1 => u(1), th2 => u(2), thdot1 => u(3), thdot2 => u(4) )
        Ekin = 0.5*m*l*l*(2*thdot1**2 + thdot2**2 + 2*thdot1*thdot2*cos(th1-th2))
        end associate
end function Ekin

pure function Etot(u); intent(in) u
	real u(4), Etot
	
	Etot = Ekin(u) + Epot(u(1),u(2))
end function Etot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATIONS TO DETERMINE CARTESIAN POSITION OF SYSTEM FOR OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! converts given angles to cartesian coordinates and writes out the
! positions of the origin, and the tips of both rods to make a pretty
! animation
subroutine animation (th1, th2)
	real th1, th2, cartesian (2,2)
	write (2,*) 0.0, 0.0
	cartesian = polar2cart(th1, th2)
	write (2,*) cartesian (1,:)
	write (2,*) cartesian (2,:)
	write (2,*) ""
	write (2,*) ""
end subroutine

! converts polar coordinates to cartesian coordinates
function polar2cart (th1, th2); intent(in) th1, th2
	real th1, th2, polar2cart(2,2), x(2), y(2)
	
	x(1) = l*sin(th1)
	x(2) = x(1) + l*sin(th2)
	y(1) = -l*cos(th1)
	y(2) = y(1) - l*cos(th2)
	
	polar2cart = reshape((/ x(1), x(2), y(1), y(2) /), (/ 2,2 /))
end function

end program
