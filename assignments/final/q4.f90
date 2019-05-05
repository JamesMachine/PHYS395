
! compile and run with: gfortran -O3 -fdefault-real-8 q4.f90 -lcfitsio -o q4 && ./q4

program q3; implicit none

  !declarations
  real,parameter :: m = 1
  real,parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
  real,parameter :: a = 1
  real,parameter :: x0 = a
  real,parameter :: V0 = 12*m
  real,parameter :: eps = 10e-002
  real,parameter :: dt = 0.025
  real,parameter :: tlength = 100
  integer,parameter :: n = floor(tlength/dt)



  ! non-parameter variables
  integer:: ii, j,i
  integer :: flag
  real :: u(2), dx(n),t(n), period, Energy

  open(1,file = "./output/output_q4/data.dat")

  write(*,*) "q4 started ..."

  ! get t(n)
  do i=1,n
    t(i) = (i-1)*dt
  end do

  flag = 0
  do ii=1,50
    u = [ (ii/10)*a, 0.0 ] !initial values for IVP
    dx = integrate(u, tlength, dt)

    ! find the first v = 0 to find period
    do j=1,n

      if (abs(dx(j)) > 0.5) then
        flag = 1
      end if

      if (flag == 1 .and. abs(dx(j)) < eps) then
          period = 2*t(j)
          Energy = E(u)
          write(*,*) "Found period: ",Energy, period
          write(1,*) Energy, period
          exit
      end if

    end do
    write(*,*) "index",ii

  end do



  write(*,*) ""
  write(*,*) ""

  close(1)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! dynamical system and Gauss-Legendre integrator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! potential
  pure function V(x); intent(in) x
    real ::V,x

    V = -V0 / (cosh(x/a))**2
  end function

  ! total energy of the dynamical system
  pure function E(u); intent(in) u
      real u(2), E

      E = 0.5*m*u(2)**2 + V(u(1))
  end function

  ! evaluate derivatives
  subroutine evalf(u, dudt)
          real u(2), dudt(2)

          associate( x => u(1), dx => u(2) )

          dudt(1) = dx
          dudt(2) = -2*V0/(m*a) * ( sinh(x/a) / (cosh(x/a))**3 )

          end associate
  end subroutine evalf

  ! 10th order implicit Gauss-Legendre integrator
  subroutine gl10(y, dt)
          integer, parameter :: s = 5, n = 2
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


  ! integrate equations of motion for a given amount of time
  function integrate(u, t, dt)
    real:: t, dt, integrate(n)
    real:: u(2), E0, Evio, period
    integer:: n
    logical::writeFlag

    ! start from a given position at rest
    E0 = E(u)

    ! number of time steps needed
    n = floor(t/dt)

    do i = 1,n

      !Evio = E(u)/E0-1

      call gl10(u,dt)
      integrate(i) = u(2)

    end do
    call gl10(u,t-n*dt)

  end function



end program
