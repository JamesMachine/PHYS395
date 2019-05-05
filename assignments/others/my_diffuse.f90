! diffuse.f90  -  one dimension diffusion equation using forward time finite difference scheme
! compile with: gfortran -O3 -fdefault-real-8 my_diffuse.f90 -o my_diffuse && ./my_dffuse

program diffusion
implicit none

!constants
integer, parameter :: fields = 1	! total number of fields to evolve
integer, parameter ::  n = 2**8	! total number of spatial points
integer, parameter :: tt = 2**12	! total number of time steps
integer, parameter :: nt = 2**3		! output solution every nt steps
integer, parameter :: nx = 2**3		! output solution every nx points

real, parameter :: alpha = 0.50		! parameter setting dt/dx^2
real, parameter :: dx = 1.0/(n-1)	! spacing between grid points
real, parameter :: dt = alpha*dx**2	! time step size


! declarations
real :: smp(fields,n,2) !smp(fields, x position, next and now)
integer :: l
open(1, file="./output/diffuse_phi.dat") ! for dump destinatio file

! initialize the problem; get initial shape (function phi)
call init(smp(:,:,1))
call dump(0.0, smp(:,:,1)) !dump (t=0.0, phi = smp(:,:,1))


! time evolution loop
do l=1,tt,2 !1,3,5,7,9 .. tt
	call step(smp(:,:,1), smp(:,:,2)); if (mod(l+0,nt) == 0) call dump((l+0)*dt, smp(:,:,2))
	call step(smp(:,:,2), smp(:,:,1)); if (mod(l+1,nt) == 0) call dump((l+1)*dt, smp(:,:,1))
end do


close(1)

contains

! initial field profile
subroutine init(phi)
	real phi(fields,n); integer i

  ! this is at t=0, phi = this shape
	forall (i=1:n) phi(:,i) = exp( -1000.0*((i-n/2)*dx)**2 )
end subroutine

! step solution forward
subroutine step(phi,phinext)
	real :: phi(fields,n)
  real :: phinext(fields,n)
  integer :: i

	! use equations of motion for interior of the grid
	do i=2,n-1
		phinext(:,i) = alpha*(phi(:,i-1)+phi(:,i+1)) + (1.0-2.0*alpha)*phi(:,i)
	end do

	! use boundary conditions at boundary points
	phi(:,1) = 0.0		! Dirichlet boundary conditions
	phi(:,n) = phi(:,n-1)	! Neumann boundary conditions
end subroutine

! dump solution at timestep t
subroutine dump(t, phi)
	real ::t
  real :: phi(fields,n)
  integer :: i

	do i=1,n,nx
		write(1,'(8g32.16e3)') t, (i-1)*dx, phi(:,i)
	end do

	write (1,*) ""
  write (1,*) ""
end subroutine

end program
