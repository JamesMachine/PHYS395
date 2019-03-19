! Assignment 4

!Running instruction for this specific file:
!gfortran -g -O3 -fdefault-real-8 -fopenmp double_pendulum.f90 -llapack -lcfitsio  -o double_pendulum && ./double_pendulum

program double_pendulum; implicit none

  !declarations
  real,parameter::m=1
  real,parameter::l=1
  real,parameter::grav=9.80665
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
  real,parameter :: x0 = 0, y0 = 0
  real,parameter :: eps = 10e-12
  real,parameter :: dt = 0.1

  integer ::i
  real :: u(4) = (/ pi/3.0, -pi/3.0, 0.0, 0.0 /) !initial values for IVP
  real:: x1,x2,y1,y2
  real::period, E0, Evio

  period = 100 * sqrt(l/grav)
  E0 = E(u)

  u = integrate(u(1),u(2),period,dt)


contains

  ! potential
  pure function V(th1,th2); intent(in) th1, th2
    real:: V, th1, th2,y1,y2

    y1 = -l/2 * cos(th1)
    y2 = - l * (cos(th1) + 0.5 * cos(th2))

    V = - 0.5 * m * grav * l * (3*cos(th1) + cos(th2))
    !V = m*grav*(y1+y2) !alternatively
  end function

  ! total energy
  pure function E(u); intent(in) u
          real u(4), E

          associate( th1 => u(1), th2 => u(2), thd1 => u(3), thd2 => u(4) )
          E = 1.0/6.0*m*grav*l**2 * (thd2**2 + 4*thd1**2 + 3*thd1*thd2*cos(th1-th2))+ V(th1,th2)
          end associate
  end function


  ! evaluate derivatives
  subroutine evalf(u, dudt)
          real u(4), dudt(4)

          !th: theta
          !thd: thetadot
          associate( th1 => u(1), th2 => u(2), thd1 => u(3), thd2 => u(4) )
          dudt(1) = thd1
          dudt(2) = thd2
          dudt(3) = -1.0 / 2.0 * l ** 2.0 * (thd1 * thd2 * sin(th1-th2) + 3 * grav / l * sin(th1))
          dudt(4) = -1.0 / 2.0 * l ** 2.0 * ( - thd1 * thd2 * sin(th1-th2) + grav / l * sin(th2))
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

  ! integrate equations of motion for a given amount of time
  function integrate(th1,th2, t, dt)
    real th1,th2,t, dt, integrate(4)
    real u(4), E0; integer n

    ! start from a given position at rest
    u = [th1, th2, 0.0, 0.0]; E0 = E(u)

    ! number of time steps needed
    n = floor(t/dt)

    !write it
    open(1,file="./output/results")
    do i = 1,n
      associate( th1 => u(1), th2 => u(2), thd1 => u(3), thd2 => u(4) )
        x1 = l/2 * sin(th1)
        y1 = -l/2 * cos(th1)

        x2 = l * (sin(th1) + 0.5 * sin(th2))
        y2 = - l * (cos(th1) + 0.5 * cos(th2))

        Evio = E(u)/E0-1

        write (1,*) (i-1)*dt,x0,y0,x1,y1,x2,y2,Evio
        call gl10(u, dt)
      end associate
    end do
    close(1)


    call gl10(u,t-n*dt)

    ! return state at time t
    integrate = u
  end function


end program
