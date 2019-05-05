! compile with: gfortran -O3 -fdefault-real-8 prob4.f90 -llapack -o a && ./a

program problem4
    implicit none
    
    ! order of the spectral scheme
    integer, parameter :: n = 100
    
    ! scale of compactification
    real, parameter :: ell = .6
    
    ! this is what you think it is...
    real, parameter :: pi = 3.1415926535897932384626433832795028842Q0
    
    ! collocation grids
    real, dimension(n) :: x, theta
    
    ! second derivative operator
    real, dimension (n,n) :: L, H
    
    real eigen_vector(n,10)
    real eigen_value(10)
    
    real lambda
    integer i, j
    
    ! initialize spectral operators
    call initg(); call initl()
    
    ! Hamiltonian in static Schrödinger equation
    H = -L/2.0; forall (i=1:n) H(i,i) = -L(i,i)/2.0 + V(x(i))
    
    !solving for psi eigenvector and energy eigenvalue 
    write(*,*) "------------------------------------------------"
    write(*,*) "Eigenvalue of the energy for harmionic potential"
    write(*,*) "------------------------------------------------"
    call lsolve()
    ! Evolve untill 10th energy eigen state and save the data
    
    open(unit=4,file="data_q4_1.dat")

    do j = 1,10
        write (4,'(4g24.16)')
        call dump(eigen_vector(:,j))
    end do
    
    close(4)


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
    ! Rayleigh itertation solver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Rayleigh's iteration solving eigenvalue problem:
    ! [-1/2 d^2/dx^2 + V(x)] psi(x) = lambda psi(x) 
    subroutine lsolve()
        real A(n,n), WR(n), WI(n), VL(n,n), VR(n, n), work(10*n)
        integer i, ind(10), status
        
     	! linear improvement inflating eigenvector
        A = H; !forall (i=1:n) A(i,i) = H(i,i) - lambda
        ! B = psi
	
        ! eigenvectors and eigenvalues
        status = 0; select case (kind(A))
            case(4); call sgeev('V', 'V', n, A, n, WR, WI, VL, n, VR, n, work, n*10, status)
            case(8); call dgeev('V', 'V', n, A, n, WR, WI, VL, n, VR, n, work, n*10, status)
            case default; call abort
        end select
        
        ind = lowest(WR,10)
        ! write(*,*) ind(2)
        

        do i = 1,10
            eigen_value(i) = WR(ind(i))
            eigen_vector(:,i) = VR(:,ind(i))
            write (*,'(1g24.16)')   eigen_value(i)
        end do
    end subroutine lsolve
    
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
    
    ! dump the solution 
    subroutine dump(eigen_vector)
        real eigen_vector(n);integer i, j
        ! open(unit=4,file="data_q4_1.dat")
        do i = 1,n
            write (4,'(4g24.16)') x(i), eigen_vector(i)
        end do
        ! close(4)
        ! write (1,*) ""; write (1,*) ""
    end subroutine
    
    end program
    