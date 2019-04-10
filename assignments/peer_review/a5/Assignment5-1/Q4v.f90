! compile with: gfortran -O3 -fdefault-real-8 Q4.f90 -llapack

program Q4
implicit none

! order of the spectral scheme as suggested
integer, parameter :: n = 150

! scale of compactification
real, parameter :: ell = 1.0

! this is what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028842Q0

! collocation grids
real, dimension(n) :: x, theta, psi

! second derivative operator
real, dimension (n,n) :: L, H

real :: lambda(n), lambda_im(n), left(n,n), right(n,n)
real :: lambda_sorted(10), temp, psi_sorted(n,10)
integer i, j, k



! initialize spectral operators
call initg(); call initl()

! Hamiltonian in static Schrödinger equation
H = -L/2.0; forall (i=1:n) H(i,i) = -L(i,i)/2.0 + V(x(i))

! directly solve for eigenvalues of H
call lsolve_eg(H,lambda,lambda_im,left,right)

! need to sort the eigenvalues. Do first one seperately
Write(*,*) "The eigenvaues for question 4 harmonic oscillator are:"
lambda_sorted(1) = minval(lambda)
psi_sorted(:,1) = right(:,minloc(lambda,1))
temp = lambda_sorted(1)

write(*,*) lambda_sorted(1)

do i = 2,10
	lambda_sorted(i) = minval(lambda,1,lambda > temp)
	psi_sorted(:,i) = right(:,minloc(lambda,1,lambda > temp))
	temp = lambda_sorted(i)
	write(*,*) lambda_sorted(i)
end do
!write(*,*) ""; write(*,*) ""

! output psi sorted from low to high order even eigenfunctions are negative, fixed manually.
do i = 1,10
	if (mod(i,2) == 1) then
		do j = 1,n
			!write(*,*) x(j), -psi_sorted(j,i)
		end do
	else
		do j = 1,n
			!write(*,*) x(j), psi_sorted(j,i)
		end do
	end if
	!write(*,*) ""; write(*,*) ""
end do



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! potential in the problem we are trying to solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! simple harmonic oscillator potential
elemental function V(x); intent(in) x
	real V, x

	V = x*x/2.0
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
! Solves for eigenvalues of Hamiltonian using sgeev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! [V(x) - 1/2 d^2/dx^2] psi(x) = lambda psi(x) is Schrodinger Eq
! => H*psi(x) = lambda*psi(x)
subroutine lsolve_eg(A,lambda,lambda_im,vl,vr)
	real :: A(n,n), lambda(n), lambda_im(n)
	real :: vl(n,n), vr(n,n), work(6*n)
	integer :: i, status

	! find eigenvalues of A (= H)
   	status = 0; select case (kind(A))
            case(4); call sgeev('V', 'V', n, A, n, lambda, lambda_im, vl, n, vr, n, work, 6*n, status)
            case(8); call dgeev('V', 'V', n, A, n, lambda, lambda_im, vl, n, vr, n, work, 6*n, status)
            case default; call abort
    end select

    ! bail at first sign of trouble
    if (status /= 0) call abort
end subroutine

! dump the solution and its residual
subroutine dump(psi)
	real :: psi(n)
	integer :: i

	do i = 1,n
		write (*,'(3g24.16)') x(i), -psi(i)
	end do

	write (*,*) ""; write (*,*) ""
end subroutine

end program
