! Only for problem 3, scanning phase space

! compile and run with: gfortran -O3 -fdefault-real-8 -fopenmp dbpened_longrun.f90 -lcfitsio -o dbpened_longrun && ./dbpened_longrun

program dbpened_longrun; implicit none

  !declarations
  real,parameter :: m1 = 1, m2 = 1
  real,parameter::  l1 = 1, l2 = 1
  real,parameter::  grav = 9.80665
  real,parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0
  real,parameter :: x0 = 0, y0 = 0
  real,parameter :: eps = 10e-12
  real,parameter :: dt = 0.1

  ! image size and scan bounds
  integer, parameter :: nx = 2**8, ny = 2**8
  real(4), parameter :: xx(2) = [-3, 3], yy(2) = [-3, 3]

  ! non-parameter variables
  integer:: i, j
  real :: u(4)
  real::period = sqrt(l1/grav)
  real::tflip, tlength

  ! data array (allocatable to avoid problems with system resources)
  real(4), allocatable :: data(:,:,:)
  allocate(data(1,nx,ny))

  !stop !de-comment this if you dont want to phase scan:

  write(*,*) "Q3: phase space scan starts..."
  tlength = 10001*period !this should be 100001*period but taking too much time (will take > 24hours..)
  !$omp parallel do
  do j = 1,ny
    do i = 1,nx
      write(*,*) 'th1(0)=',xx(1) + (xx(2)-xx(1))*(i-1)/(nx-1),"th2(0)=",yy(1) + (yy(2)-yy(1))*(j-1)/(ny-1)
      tflip = integrate2(xx(1) + (xx(2)-xx(1))*(i-1)/(nx-1), yy(1) + (yy(2)-yy(1))*(j-1)/(ny-1), tlength, dt)
      data(:,i,j) = log10(tflip) / tlength
    end do
  end do
  !$omp end parallel do
  write(*,*) "phase space scan ends..."

  ! write out image to file
  call write2fits('./output_longrun/data.fit', data, xx, yy, ['th1 ','th2 ','thd1','thd2'], '(th1(0),th2(0))')

  deallocate(data)


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! dynamical system and Gauss-Legendre integrator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! potential for 2D particle motion
  pure function V(th1,th2); intent(in) th1, th2
    real ::V,th1, th2

    V = - (m1 + m2)*grav*l1*cos(th1) - m2*grav*l2*cos(th2)
  end function

  ! total energy of the dynamical system
  pure function E(u); intent(in) u
          real u(4), E

          associate( th1 => u(1), th2 => u(2), thd1 => u(3), thd2=>u(4) )
          E = 0.5*m1*l1**2*thd1**2 + 0.5*m2*(l1**2*thd1**2+l2**2*thd2**2+2*l1*l2*thd1*thd2*cos(th1-th2)) + V(th1,th2)
          end associate
  end function

  ! evaluate derivatives
  subroutine evalf(u, dudt)
          real u(4), dudt(4), p1,p2, c1, c2

          associate( th1 => u(1), th2 => u(2), thd1 => u(3), thd2 => u(4) )

          dudt(1) = thd1
          dudt(2) = thd2

          dudt(3) =(-grav*(2*m1+m2)*sin(th1) - m2*grav*sin(th1-2*th2)-2*sin(th1-th2)*m2*(thd2**2*l2+thd1**2*l1*cos(th1-th2))) &
          / ( l1*(2*m1+m2-m2*cos(2*th1-2*th2)) )

          dudt(4) = (2*sin(th1-th2)*(thd1**2*l1*(m1+m2)+grav*(m1+m2)*cos(th1)+thd2**2*l2*m2*cos(th1-th2))) &
          / ( l2*(2*m1+m2-m2*cos(2*th1 - 2*th2)) )

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
  function integrate(th1, th2, t, dt, writeFlag)
    real:: th1, th2, t, dt, integrate(4)
    real:: u(4), E0, Evio
    real:: x1,x2,y1,y2
    integer:: n
    logical::writeFlag

    ! start from a given position at rest
    u = [th1, th2, 0.0, 0.0]; E0 = V(th1,th2)

    ! number of time steps needed
    n = floor(t/dt)

    open(1,file = "./output/results")
    do i = 1,n
      associate( th1 => u(1), th2 => u(2), thd1 => u(3), thd2 => u(4) )
          x1 = l1 * sin(th1)
          y1 = -l1 * cos(th1)

          x2 = l1*sin(th1) + l2*sin(th2)
          y2 = -l1*cos(th1) -l2*cos(th2)

          Evio = E(u)/E0-1

          call gl10(u,t-n*dt)

          !if writeFlag is true, write the result into results file
          if(writeFlag) then
            write (1,*) (i-1)*dt,x0,y0,x1,y1,x2,y2,Evio
          end if
      end associate
    end do
    close(1)

    call gl10(u,t-n*dt)

    ! return state at time t
    integrate = u
  end function


  ! integrate equations of motion for a given amount of time
  function integrate2(th1, th2, t, dt)
    real:: th1, th2, t, dt,integrate2
    real:: u(4)
    integer:: n

    ! start from a given position at rest
    u = [th1, th2, 0.0, 0.0]

    ! number of time steps needed
    n = floor(t/dt)

    do i = 1,n
      associate( th1 => u(1), th2 => u(2), thd1 => u(3), thd2 => u(4) )

        if (abs(th1)>pi .or. abs(th2) > pi) then
          integrate2 = (i-1)*dt
          return
        end if
        call gl10(u,t-n*dt)

      end associate
    end do
    call gl10(u,t-n*dt)

    ! return state at time t
    integrate2 = (i-1)*dt
  end function


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write array data into FITS file as sequence of image extensions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !using FITSIO library
  subroutine write2fits(file, array, xx, yy, vars, coords)
          character(len=*) file, vars(:), coords
          real(4) array(:,:,:); real(4) xx(2), yy(2)
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
