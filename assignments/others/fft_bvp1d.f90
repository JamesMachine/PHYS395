! FFTW demo, compile with

program fft_bvp1d; implicit none

include "fftw3.f"

integer, parameter :: N = 2**22, NN = N/2 + 1
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0, dk = pi

real(8) in(N), out(N), dst(N), x(N)
integer*8 plan

integer i, k

!in = (0,0)
!in(3*N/8:5*N/8) = 1.0
forall (i=1:N) in(i) = exp( -1000.0*((i-N/2)*(1.0/N))**2 )
forall (i=1:N) x(i) = (i-1)*(1.0/N)


call dfftw_plan_r2r_1d(plan,N,in,dst,FFTW_RODFT00,FFTW_ESTIMATE)
call dfftw_execute_r2r(plan, in, dst)
call dfftw_destroy_plan(plan)

do k = 1,N
	dst(k) = -dst(k)/(k*dk)**2	! inverse Laplacian
end do

call dfftw_plan_r2r_1d(plan,N,dst,out,FFTW_RODFT00,FFTW_ESTIMATE)
call dfftw_execute_r2r(plan, dst, out)
call dfftw_destroy_plan(plan)

do i = 1,N,1024
	write (*,*) x(i), in(i), out(i)/(2*(N+1))
end do

end program
