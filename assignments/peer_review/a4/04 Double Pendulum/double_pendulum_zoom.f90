! double_pendulum_flip.f90  -  simple chaotic evolution demo
! compile with: gfortran -O3 -fdefault-real-8 -fopenmp double_pendulum_zoom.f90 -lcfitsio

program double_pendulum_flip; implicit none

! physical constants
real, parameter :: pi = 4.D0*DATAN(1.D0)
real, parameter :: g = 9.8

! parameters of the system
real, parameter :: m = 1.0
real, parameter :: l = 1.0

! characteristic timestep
real, parameter :: dt = sqrt(l/g)/100.0
real thrange1(2), thrange2(2), hold(4)
integer, parameter :: nth = 2**4

real th1, th2
real(4) :: data(1,nth,nth) = 0.0
integer i, j
integer :: status = 0

character(len=8000) :: arg

do i = 1,4
	call get_command_argument(i, arg)
	write (*,*) trim(arg)
	read (arg,*,iostat=status) hold(i)
end do

thrange2 = (/ hold(1), hold(2) /)
thrange1 = (/ hold(3), hold(4) /)

write (*,*) thrange1
write (*,*) thrange2

! abuse the fact that the system has inversion symmetry. For speed, turn
! down the number of time steps and/or nth.
!$omp parallel do
do i = 1,nth
	th2 = thrange2(2) - (thrange2(2) - thrange2(1)) * (i-1)/(nth-1)
	do j = 1,nth
		th1 = thrange2(1) + (thrange2(2) - thrange2(1)) * (j-1)/(nth-1)
		data(1,i,j) = integrate(th1, th2, 10000.0*sqrt(l/g), dt)/100.0
		write(*,*) i,j, data(1,i,j)
	end do
end do

do i = 1,nth
!	write (2,*) data(1,i,:)
	write (1,*) log (1+data(1,i,:))
end do

! write out image to file
call write2fits('data.fit', log(1+data), thrange1, thrange2, ['iterations'], '(x,y)')

contains

function integrate(th1, th2, t, dt)
	real th1, th2, t, dt, integrate
	real u(4), E0; integer n, k
	
	! start from a given position at rest
	u = [th1, th2, 0.0, 0.0]; E0 = Etot(u)
	
	! maximal number of time steps needed
	n = floor(t/dt)
	
	! Track the number of steps required for a flip, bailing out early
	! if the initial energy isn't enough for a flip.
	if (E0 <= -g*l*m) then
		k = n
	else
		do k = 1,n
			call gl10(u, dt)
			if (abs(u(2)) > pi .or. abs(u(1)) > pi) exit
	!		write (*,'(4g24.16)') i*dt, u(1), u(2), Etot(u)/E0 - 1.0
	!		call animation(u(1),u(2))
		end do
	end if

!	call gl10(u,t-n*dt)
	
	integrate = k
end function

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

! converts given angles to cartesian coordinates and writes out the
! positions of the origin, and the tips of both rods to make a pretty
! animation
subroutine animation (th1, th2)
	real th1, th2, cartesian (2,2)
	write (*,*) 0.0, 0.0
	cartesian = polar2cart(th1, th2)
	write (*,*) cartesian (1,:)
	write (*,*) cartesian (2,:)
	write (*,*) ""
	write (*,*) ""
end subroutine

! potential for double pendulum
pure function Epot(th1,th2); intent(in) th1, th2
	real Epot, th1, th2
	
	Epot = -g*l*m*(2*cos(th1)+cos(th2))
end function Epot

! total energy of the dynamical system
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

function polar2cart (th1, th2); intent(in) th1, th2
	real th1, th2, polar2cart(2,2), x(2), y(2)
	
	x(1) = l*sin(th1)
	x(2) = x(1) + l*sin(th2)
	y(1) = -l*cos(th1)
	y(2) = y(1) - l*cos(th2)
	
	polar2cart = reshape((/ x(1), x(2), y(1), y(2) /), (/ 2,2 /))
end function

subroutine write2fits(file, array, xx, yy, vars, coords)
        character(len=*) file, vars(:), coords
        real(4) array(:,:,:); real(8) xx(2), yy(2)
        optional xx, yy, vars, coords
        
        integer i, j, status, unit
        integer :: hdus, naxis = 2, n(2), npix
        integer :: bitpix = -32, group = 1, blocksize = -1
        
        ! data dimansions
        hdus = size(array,1)
        n(1) = size(array,2)
        n(2) = size(array,3)
        npix = n(1)*n(2)
        
        ! delete file if it already exists
        open(unit=1234, iostat=status, file=file, status='old')
        if (status == 0) close(1234, status='delete'); status = 0
        
        ! initialize FITS file
        call ftgiou(unit, status)
        call ftinit(unit, file, blocksize, status)
        
        ! write image extensions
        do i = 1,hdus
                call ftiimg(unit, bitpix, naxis, n, status)
                call ftppre(unit, group, 1, npix, array(i,:,:), status)
                
                if (present(vars)) then
                        if (present(coords)) then
                                call ftpkys(unit, 'EXTNAME', trim(vars(i))//coords, 'variable stored in extension', status)
                        else
                                call ftpkys(unit, 'EXTNAME', trim(vars(i)), 'variable stored in extension', status)
                        end if
                end if
                if (present(xx)) then
                        call ftpkyj(unit, 'CRPIX1', 1, 'x-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL1', xx(1), 14, 'x-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT1', (xx(2)-xx(1))/(n(1)-1), 14, 'x-axis increment', status)
                end if
                if (present(yy)) then
                        call ftpkyj(unit, 'CRPIX2', 1, 'y-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL2', yy(1), 14, 'y-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT2', (yy(2)-yy(1))/(n(2)-1), 14, 'y-axis increment', status)
                end if
        end do
        
        ! clean up
        call ftclos(unit, status)
        call ftfiou(unit, status)
end subroutine write2fits

end program
