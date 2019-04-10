! For Assignment 5 -  Q1 -  solve BVP using shooting method
! compile with: gfortran -O3 -fdefault-real-8 shrodinger_shooting_q1.f90 -o shrodinger_shooting_q1 && ./shrodinger_shooting_q1

program shrodinger_shooting_q1
implicit none

!define constants
real,parameter :: eps = 1.0e-12
real,parameter :: m = 1.0
real,parameter :: hbar = 1.0
real,parameter :: omega = 1.0
real,parameter :: E = 4.2
real,parameter :: xrange = 4.0
real,parameter :: dx = 0.001
integer,parameter :: n = floor(xrange/dx)


!declarations
real :: x0
real :: psi0
real :: dpsi0
real :: initCond(3)
real :: psi_pos(2,n)
real :: psi_neg(2,n)
real :: psi(2,2*n-1) !including negative and postive x range; 2 dim - time and actual value
real :: psi_even(2,2*n-1)
real :: psi_odd(2,2*n-1)
real :: norm
integer :: i



! Step1: Finding inital conditions for given BVP and E
write(*,* ) "Q1 Started..."
dpsi0 = bisect(-10.0,10.0)
write(*,*) "Found initial condition parameter dpsi(0)/dx = ",dpsi0
write(*,*) "For E = ",E
write(*,*) "With psi(0) = 0.0"
write(*,*) ""


! Step2: Get wave function in both postive and negative direction
x0 = 0.0
psi0 = 1.0 !boundary condtion psi(0) = 1.0
initCond = [x0, psi0, dpsi0]
psi_pos = integrate2(initCond, xrange, dx)
psi_neg = integrate2(initCond, xrange, -dx)

write(*,*) "Wave function plotted with:"
write(*,*) "x0 = ", x0
write(*,*) "psi(0) = ", psi0
write(*,*) "dpsi(0) = ", dpsi0

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

!Step3.5: Normalization!!
norm = sum(psi(2,:)*psi(2,:))*dx
psi(2,:) = psi(2,:) / sqrt(norm)

!Step4: construct the psi_even and psi_odd wave function from psi
psi_even(1,:) = psi(1,:)
psi_odd(1,:) = psi(1,:)
do i=1,2*n-1
  psi_even(2,i) = 0.5 * ( psi(2,(2*n-1)-i+1) + psi(2,i) )
  psi_odd(2,i) = 0.5 * ( psi(2,(2*n-1)-i+1) - psi(2,i) )
end do

open(1, file = "./output/output_q1/output_q1_psi")
do i = 1, 2*n-1
  write(1,*) psi(1,i),psi_even(2,i),psi_odd(2,i),psi(2,i)
end do
close(1)


!write(*,*) sum(psi(2,:)*psi(2,:))*dx
write(*,* ) ""
write(*,* ) ""

contains

  ! potential
  pure function V(x); intent(in) x
    real V, x

    V = 0.5 * m * omega**2 * x**2
  end function


  ! evaluate derivatives
  subroutine evalf(u, dudx) ! need to include parameter x?
          real u(3), dudx(3)

          associate(x => u(1), psi => u(2), dpsi => u(3))
            dudx(1) = 1.0
            dudx(2) = dpsi
            dudx(3) = 2*m/hbar**2 * (V(x) - E)*psi
          end associate
  end subroutine evalf

  ! 10th order implicit Gauss-Legendre integrator
  subroutine gl10(y, dt)
          integer, parameter :: s = 5, n = 3
          real y(n), g(n,s), dt, x; integer i, k

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

  ! integrate for finding initial condtions
  function integrate(initCond, t, dt)
    real :: initCond(3)
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
      call gl10(u, dt)
    end do

    integrate = u(3)

  end function



  ! integrate for returning wave function
  function integrate2(initCond, t, dt)
    real :: initCond(3)
    real :: t
    real :: dt
    real :: integrate2(2,n)
    real :: u(3)
    integer :: i
    integer :: n

    ! start from zero at a given gradient
    u = [initCond(1), initCond(2), initCond(3)]

    ! number of time steps needed
    n = floor(t/abs(dt))

    do i = 1,n
      call gl10(u, dt)
      integrate2(1,i) = u(1)
      integrate2(2,i) = u(2)
    end do

  end function



  ! function to find a root of...
  function f(x); intent(in) x
    real :: initCond(3)
    real :: f
    real :: x !this is dpsi0

    initCond = [0.0, 1.0, x]
    f = integrate(initCond, xrange, dx)
  end function


  function bisect(a0,b0); intent(in) ::a0,b0
    real ::a0, b0, a, b, c, fa, fb, fc
    integer ::i
    real ::bisect


    ! initial interval
    a = a0; fa = f(a)
    b = b0; fb = f(b)

    if (fa*fb > 0.0) then
      write(*,*) "The root is not bracketed, bailing out..."
      return
    end if

    fc = 1 ! initialize fc to 1, because fortran contains small value in fc before initializing it
    ! bisect the interval
    do while(abs(fc) > eps)
      c = (a+b)/2.0; fc = f(c); if (fc == 0.0) exit

      ! Ridder's variation on the basic method
      c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(c)

      if (fa*fc < 0.0) then; b = c; fb = fc; end if
      if (fc*fb < 0.0) then; a = c; fa = fc; end if
    end do

    bisect = c
  end function


end program
