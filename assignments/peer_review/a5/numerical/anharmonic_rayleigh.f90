! anharmonic_rayleigh.f90 - solve Schrodinger's equation using dgeev()
! for a quartic potential
! fort.1 is eigenfunctions
! fort.4 is eigenvalues
! compile with: gfortran -O3 -fdefault-real-8 anharmonic_rayleigh.f90 -llapack

program anharmonic_rayleigh
implicit none

! order of the spectral scheme
integer, parameter :: n = 600

! scale of compactification
real, parameter :: ell = 5.0

! this is what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028842Q0

! collocation grids
real, dimension(n) :: x, theta

! second derivative operator
real, dimension (n,n) :: L, H

real psi(n,10)
real eigenvalue(10)

real lambda
integer i, k

! initialize spectral operators
call initg(); call initl()

! Hamiltonian in static Schrödinger equation
H = -L/2.0; forall (i=1:n) H(i,i) = -L(i,i)/2.0 + V(x(i))

! wavefunction
call esolve()

do i = 1,10
	call dump(psi(:,i))
end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! potential in the problem we are trying to solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! simple harmonic oscillator potential
elemental function V(x); intent(in) x
	real V, x, lambda
	
	lambda = 1.0
	
	V = x**4 * lambda /4.0
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spectral grid, basis functions and derivative operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize the collocation grid
subroutine initg()
	integer i
	
	forall (i=1:n) theta(i) = pi*(n-i+0.5)/n; x = ell/tan(theta)
end subroutine

! evaluate rational Chebyshev basis on a grid theta
subroutine evalb(n, pts, theta, Tn, Tnx, Tnxx)
        integer n, pts; real, dimension(pts), intent(in) :: theta
        real, dimension(pts), intent(out), optional :: Tn, Tnx, Tnxx
        
        ! Chebyshev basis and its derivatives
        if (present(Tn))   Tn = cos(n*theta)
        if (present(Tnx))  Tnx = n * sin(n*theta) * sin(theta)**2/ell
        if (present(Tnxx)) Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3/ell**2
end subroutine evalb

! initialize linear spectral derivative operator
subroutine initl()
	integer i, pivot(n), status; real A(n,n), B(n,n)
	
	! evaluate basis and differential operator values on collocation grid
	do i = 1,n
		call evalb(i-1, n, theta, Tn=A(i,:), Tnxx=B(i,:))
	end do
	
        ! find linear operator matrix
        status = 0; select case (kind(A))
                case(4); call sgesv(n, n, A, n, pivot, B, n, status)
                case(8); call dgesv(n, n, A, n, pivot, B, n, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        L = transpose(B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rayleigh itertation solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solve eigenvalue problem:
! [-1/2 d^2/dx^2 + V(x)] psi(x) = E psi(x)
subroutine esolve()
	real A(n,n), WR(n), WI(n), VL(n,n), VR(n, n), work(10*n)
	integer i, ind(10), status
	
	A = H
	
	! eigenvectors and eigenvalues
	status = 0; select case (kind(A))
		case(4); call sgeev('V', 'V', n, A, n, WR, WI, VL, n, VR, n, work, n*10, status)
		case(8); call dgeev('V', 'V', n, A, n, WR, WI, VL, n, VR, n, work, n*10, status)
		case default; call abort
	end select
	
	ind = lowest(WR,10)
	do i = 1,10
		eigenvalue(i) = WR(ind(i))
		write (4,'(1g24.16)') eigenvalue(i)
		psi(:,i) = VR(:,ind(i))
	end do
end subroutine esolve

! Get the position of the m lowest elements in a vector
function lowest(B,m)
	real B(n), large, minimum
	integer i, j, lowest(m), m
	
	minimum = 0.0
	
	do i = 1,10
		lowest(i) = minloc(B, DIM = 1, MASK = B > minimum)
		minimum = B(lowest(i))
	end do
end function lowest

! dump the solution and its residual
subroutine dump(psi)
	real psi(n), delta(n), norm; integer i, j
	
	norm = 0.0
	do j = 2,n
		norm = norm + (psi(j-1)**2 + psi(j)**2)/2.0 * (x(j)-x(j-1))
	end do
		
	do i = 1,n
		write (1,'(4g24.16)') x(i), psi(i), psi(i)*psi(i)/norm
	end do
	
	write (1,*) ""; write (1,*) ""
end subroutine

end program
