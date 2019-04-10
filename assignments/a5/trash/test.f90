! For Assignment 1 -  Q1 -  solve BVP using shooting method
! compile with: gfortran -O3 -fdefault-real-8 shrodinger_shooting_q1.f90 -o shrodinger_shooting_q1 && ./shrodinger_shooting_q1

program shrodinger_shooting_q1
implicit none

!define constants
real,parameter :: eps = 1.0e-5
real,parameter :: m = 1.0
real,parameter :: hbar = 1.0
real,parameter :: omega = 1.0
real,parameter :: xrange = 4
real,parameter :: dx = 0.001
integer,parameter :: n = floor(xrange/dx)


!declarations
real :: psi_pos(2,n)
real :: psi_neg(2,n)
real :: psi(2,2*n-1) !including negative and postive x range; 2 dim - time and actual value
real :: psi_array(10,2,2*n-1)
real :: E
real :: E_array(10)
!real :: psiend,x
real :: initCond(3)
integer :: i,ii
real :: real_i



!Step1: finding inital conditions for given BVP and intial condtions
write(*,* ) "Q2 Started..."

!Step1: find the first 10 eigenvalues for even wave function for iteration
do i=1,10
  !real_i = real(i)
  if (mod(i-1,2) == 0) then !when the index is even. Note I subtract 1 because fortran starts index at 1
    !find the eigenvalues for even wave function
    initCond = [0.0, 1.0, 0.0]
    E_array(i) = bisect2(real(i-1), real(i), initCond) ! store at even index 0,2,4,6,8
  else
    !find the eigenvalues for odd wave function
    initCond = [0.0, 0.0, 1.0]
    E_array(i) = bisect2(real(i-1), real(i), initCond) !store at odd index 1,3,5,7,9
  end if
end do



!Step2: Get wave function in both postive and negative direction
initCond = [0.0, 1.0, 0.0]
psi_pos = integrate2(initCond, E_array(1), xrange, dx)
psi_neg = integrate2(initCond, E_array(1), xrange, -dx)


!Step3: Construct the psi for both postive and negative direction
do i=1,2*n-1
  ! for postive range
  if (i <= n) then
    psi(1,i) = psi_neg(1,n-i+1)
    psi(2,i) = psi_neg(2,n-i+1)
  else
    psi(1,i) = psi_pos(1,i-n)
    psi(2,i) = psi_pos(2,i-n)
  end if
end do

open(1, file = "./output/q2_psi")
do i = 1, 2*n-1
  write(*,*) psi(1,i),psi(2,i)
  write(1,*) psi(1,i),psi(2,i)
end do
close(1)

write(*,* ) ""
write(*,* ) ""

write(*,*) sum(psi(2,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!Functions & subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  ! potential
  pure function V(x); intent(in) x
    real V, x

    V = 0.5 * m * omega**2 * x**2
  end function


  ! evaluate derivatives
  subroutine evalf(u, dudx, E) ! need to include parameter x?
          real u(3), dudx(3), E

          associate(x => u(1), psi => u(2), dpsi => u(3))
            dudx(1) = 1.0
            dudx(2) = dpsi
            dudx(3) = 2*m/hbar**2 * (V(x) - E)*psi
          end associate
  end subroutine evalf

  ! 10th order implicit Gauss-Legendre integrator
  subroutine gl10(y, dt, E)
          integer, parameter :: s = 5, n = 3
          real y(n), g(n,s), dt, E; integer i, k

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
                          call evalf(y + g(:,i)*dt, g(:,i), E)
                  end do
          end do

          ! update the solution
          y = y + matmul(g,b)*dt
  end subroutine gl10

  ! integrate for finding eigenvalues
  function integrate(initCond, E, t, dt)
    real :: initCond(3)
    real :: E
    real :: t
    real :: dt
    real :: integrate
    real :: u(3)
    integer :: i
    integer :: n

    ! start from zero at a given gradient
    u = [initCond(1), initCond(2), initCond(3)]

    ! number of time steps needed
    n = floor(t/abs(dt))

    do i = 1,n
      call gl10(u, dt, E)
    end do

    integrate = u(3) !make sure you return last derivative --> end of psi value

  end function


  ! integrate for returning wave function
  function integrate2(initCond, E, t, dt)
    real :: initCond(3)
    real :: t
    real :: dt
    real :: integrate2(2,n)
    real :: u(3)
    real :: E
    integer :: i
    integer :: n

    ! start from zero at a given gradient
    u = [initCond(1), initCond(2), initCond(3)]

    ! number of time steps needed
    n = floor(t/abs(dt))

    do i = 1,n
      call gl10(u, dt, E)
      integrate2(1,i) = u(1)
      integrate2(2,i) = u(2)
    end do

  end function



  ! function to find a root of even function...
  function f(x, initCond); intent(in) x, initCond
    real :: f
    real :: x
    real :: initCond(3)
    f = integrate(initCond, x, xrange, dx)
  end function


  function bisect2(a0,b0,initCond); intent(in) ::a0,b0,initCond
    real ::a0, b0, a, b, c, fa, fb, fc
    integer ::i
    real :: initCond(3)
    real ::bisect2


    ! initial interval
    a = a0; fa = f(a,initCond)
    b = b0; fb = f(b,initCond)
    !write(*,*) fa, fb

    if (fa*fb > 0.0) then
      write(*,*) "The root is not bracketed, bailing out..."
      return
    end if

    fc = 1 ! initialize fc to 1, because fortran contains small value in fc before initializing it
    ! bisect the interval
    do while(abs(fc) > eps)
      c = (a+b)/2.0; fc = f(c,initCond); if (fc == 0.0) exit

      ! Ridder's variation on the basic method
      c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(c,initCond)

      if (fa*fc < 0.0) then; b = c; fb = fc; end if
      if (fc*fb < 0.0) then; a = c; fa = fc; end if
    end do

    bisect2 = c
  end function


end program
