! Assignment 5 - Q4
!gfortran -O3 -fdefault-real-8 shrodinger_dgeev_q4.f90 -o shrodinger_dgeev_q4 && ./shrodinger_dgeev_q4

program shrodinger_dgeev_q4; implicit none


  ! define constants
  integer, parameter :: n = 150   ! order of the spectral scheme
  real, parameter :: ell = 1.0   ! scale of compactification
  real, parameter :: pi = 3.1415926535897932384626433832795028842Q0


  ! declarations
  real, dimension(n) :: x, theta   ! collocation grids
  real, dimension (n,n) :: L, H   ! second derivative operator
  integer i, k,min_i
  real::WR(n), VR(n,n)
  !WR: eigenvalue array
  !VR: psi array


  write(*,*) "Q4 Started..."

  ! initialize spectral operators L
  call initg()
  call initl()

  ! Get Hamiltonian in static Schrödinger equation
  H = -L/2.0; forall (i=1:n) H(i,i) = -L(i,i)/2.0 + V(x(i))

  ! solve eigenvalue problem using dgeev
  call esolve(H, WR, VR)

  ! write into files
  call dump(x,WR,VR)

  write(*,*) "Found Eigenvalue E:"
  do i=1,10
    write(*,*) WR(i)
  end do


  write(*,*) ""
  write(*,*) ""


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
  ! lapack eigen solver
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Solving eigenvalue problem using dgeev lapack library:
  ! H * psi(x) = lambda * psi(x)
  subroutine esolve(A,WR,VR)
    character::JOBVL,JOBVR
    real::A(n,n),WR(n), WI(n),VL(n,n),VR(n,n),WORK(6*n)
    integer::LDA,LDVL,LDVR,LWORK,status
    real::temp, temp_array(n)

    JOBVL='N'
    JOBVR='V'
    LDA=n
    LDVL = n
    LDVR = n
    LWORK = 6*n
    status = 0

    select case (kind(A))
      case(4); call dgeev(JOBVL, JOBVR, n, A, LDA, WR, WI,VL,LDVL,VR,LDVR,WORK,LWORK,status)
      case(8); call dgeev(JOBVL, JOBVR, n, A, LDA, WR, WI,VL,LDVL,VR,LDVR,WORK,LWORK,status)
      case default; call abort
    end select

    ! bail at first sign of trouble
    if (status /= 0) call abort

    ! Originally,each column corresponds to each psi,
    ! but make VR as each row corresponds to each psi by taking transpose
    VR=transpose(VR)



    !sort WR and VR from lower value to larger value
    do i = 1,n
      ! find the index at minimum value for WR
      min_i = minloc(WR(i:n), dim=1) + i - 1

      !swap this with the first row for both WR and VR
      temp = WR(i)
      WR(i) = WR(min_i)
      WR(min_i) = temp

      temp_array(:) = VR(i,:)
      VR(i,:) = VR(min_i,:)
      VR(min_i,:) = temp_array(:)
    end do

    !Global phase correction of -1. Lapack gives phase error in giving eigenstate
    do i=1,n
      if (VR(i,n/2+1) < 0) then
        VR(i,:)=-VR(i,:)
      end if
    end do

  end subroutine



  subroutine dump(x,WR,VR)
    real::x(n),WR(n), VR(n,n)
    integer::i,j

    open(1,file="./output/output_q4/q4_E")
    ! write all the eigenvalues
    do i = 1,n
      write (1,'(4g24.16)') WR(i)
    end do
    close(1)

    ! write all the psi
    open(1,file="./output/output_q4/q4_psi")
    do i=1,n
      do j=1,n
        write (1,'(4g24.16)') x(j),VR(i,j), VR(i,j)**2
      end do
      !write(*,*) sum(VR(i,:)*VR(i,:)) !normalization factor
      write (1,*) ""; write (1,*) ""
    end do
    close(1)
  end subroutine


end program
