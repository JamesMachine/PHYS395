module utilities
implicit none

real, parameter :: PI = 3.141592653589793238462643383279502884197169399375
real, parameter :: e = 2.7182818284590452353602874713527
real, parameter :: e_charge = 1.6021766208e-19 ! in coulombs
real, parameter :: hbar = 1.054571800e-34 ! in J.s
real, parameter :: G = 6.67408e-11 ! in m^3 kg^{-1} s^{-2}
real, parameter :: kb = 1.38064852e-23 ! in J K^{-1}

contains

! ============================================================================
! ==============
! ==============
! ============== LAPACK WRAPPERS
! ==============
! ==============
! ============================================================================

subroutine lsolve(n, A, B)
    integer, intent(in) :: n
    real, dimension(n,n), intent(inout) :: A
    real, dimension(n,1), intent(inout) :: B
    integer, dimension(n) :: pivot
    integer :: status

    ! initialize status to all clear
    status = 0

    ! call appropriate version of xGESV
    select case (kind(A))
        case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
        case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
        case default; call abort
    end select

    ! abort at first sign of trouble
    if (status /= 0) stop "singular matrix in lsolve()"
end subroutine

subroutine svd(n, A, B)
    integer, intent(in) :: n
    real, dimension(n,n), intent(inout) :: A
    real, dimension(n,1), intent(inout) :: B
    real, dimension(n) :: pivot
    integer :: status

    ! initialize status to all clear
    status = 0

    ! call appropriate version of xGESV
    select case (kind(A))
        case(4); call sgesvd(n, 1, A, n, pivot, B, n, status)
        case(8); call dgesvd(n, 1, A, n, pivot, B, n, status)
        case default; call abort
    end select

    ! abort at first sign of trouble
    if (status /= 0) stop "singular matrix in lsolve()"
end subroutine svd

subroutine svdsolve(n, A, B, epsilon)
    ! solve A.x = B using LAPACK xGESVD
    ! A gets destroyed, answer is returned in B
    integer, intent(in) :: n
    real, dimension(n,n), intent(inout) :: A
    real, dimension(n), intent(inout) :: B
    real, intent(in) :: epsilon
    integer :: status

    ! SVD matrices
    real, dimension(n,n) :: U, V
    real, dimension(n) :: S
    real, dimension(6*n) :: W

    ! initialize status to all clear
    status = 0

    ! call appropriate version of xGESV
    select case (kind(A))
        case(4); call sgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
        case(8); call dgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
        case default; call abort
    end select

    ! abort at first sign of trouble
    if (status /= 0) stop "singular matrix in svdsolve()"

    ! compute the solution using pseudo-inverse
    B = matmul(transpose(U),B)
    where (S > epsilon*S(1))
        B = B/S
    elsewhere
        B = 0.0
    end where
    B = matmul(transpose(V),B)
end subroutine svdsolve

subroutine svdsolve_with_vals(n, S, A, B, epsilon)
    ! variation on svdsolve that returns the singular values as well
    integer, intent(in) :: n
    real, dimension(n), intent(inout) :: S
    real, dimension(n, n), intent(inout) :: A
    real, dimension(n), intent(inout) :: B
    real :: epsilon
    integer :: status

    ! SVD matrices
    real, dimension(n,n) ::  U
    real, dimension(n,n) :: V
    real, dimension(6*n) :: W

    ! initialize status to all clear
    status = 0

    ! call appropriate version of xGESV
    select case (kind(A))
        case(4); call sgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
        case(8); call dgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
        case default; call abort
    end select

    ! abort at first sign of trouble
    if (status /= 0) stop "singular matrix in svdsolve_with_vals()"

    ! compute the solution using pseudo-inverse
    B = matmul(transpose(U),B)
    where (S > epsilon*S(1))
        B = B/S
    elsewhere
        B = 0.0
    end where
    B = matmul(transpose(V),B)
end subroutine svdsolve_with_vals

! minimize |A.x - B|^2 using LAPACK canned SVD routine (A gets destroyed, the answer is returned in B)
! rcond determines the effective rank of A as described in LAPACK docs. Pass -1.0 for machine precision
subroutine lstsq(num_grid_pts, max_order, A, B)
    integer, intent(in) :: num_grid_pts
    integer, intent(in) :: max_order
    real, dimension(num_grid_pts, max_order), intent(inout) :: A
    real, dimension(num_grid_pts), intent(inout):: B
    integer :: l
    real, dimension(max_order) :: S
    real, dimension(6*num_grid_pts) :: W
    integer :: rank, status
    real, parameter :: rcond = -1.0


    l = 3*num_grid_pts + max_order

    ! find solution by singular value decomposition
    status = 0; select case (kind(A))
            case(4); call sgelss(num_grid_pts, max_order, 1, A, num_grid_pts, B, num_grid_pts, S, rcond, rank, W, l, status)
            case(8); call dgelss(num_grid_pts, max_order, 1, A, num_grid_pts, B, num_grid_pts, S, rcond, rank, W, l, status)
            case default; call abort
    end select

    ! bail at first sign of trouble
    if (status /= 0) call abort
end subroutine lstsq

! ============================================================================
! ==============
! ==============
! ============== TEMPLATES
! ==============
! ==============
! ============================================================================

! basis functions we are fitting
! subroutine evalb(x, n, B)
!     integer, intent(in) :: n
!     real, intent(in) :: x
!     real, dimension(n), intent(inout) :: B
!     integer :: i

!     forall (i=1:n) B(i) = cos(2*PI*x*i) ! set basis function here
! end subroutine evalb

! evaluate polynomial expansion
! TEMPLATE FOR DUMPING FIT FUNCTIONS
! subroutine dump_expression(n, x_min, x_max, m, U)
!     integer, intent(in) :: n
!     integer, intent(in) :: m
!     real, intent(in) :: x_min, x_max
!     real, dimension(n) :: U
!     real, dimension(m) :: x_pts
!     real, dimension(n) :: basis
!     integer :: i

!     call linspace(x_pts, x_min, x_max, m)

!     do i = 1,m
!         basis = 0.0
!         call evalb(x_pts(i), n, basis)
!         write (*,*) x_pts(i), sum(U*basis)
!     end do
! end subroutine dump_expression

! ============================================================================
! ==============
! ==============
! ============== SELF MADE MODULES
! ==============
! ==============
! ============================================================================


elemental function ChebyshevT(x, n)
    real :: ChebyshevT
    real, intent(in) :: x
    integer, intent(in) :: n
    ChebyshevT = cos(n*acos(x))
end function ChebyshevT

elemental function DChebyshevT(x, n)
    real :: DChebyshevT
    real, intent(in) :: x
    integer, intent(in) :: n
    DChebyshevT = n*sin(n*acos(x)) / sqrt(1-x*x)
end function DChebyshevT

subroutine linspace(array, start_point, stop_point, n)
    real, dimension(n, 1), intent(inout) :: array
    real, intent(in) :: start_point, stop_point
    real :: step
    integer, intent(in) :: n
    integer :: i
    ! set up interval, sampled uniformly from [-1.0, 1.0]
    step = (stop_point-start_point) / (n-1.0)
    do i=1,n
        array(i,1) = start_point + (i-1)*step
    end do
end subroutine linspace

subroutine swap(a,b)
    ! subroutine to swap two vectors with one another
    real, dimension(:), intent(inout) :: a,b
    real, dimension(size(a)) :: dummy
    dummy=a
    a=b
    b=dummy
end subroutine swap


subroutine gaussj(n, a, b)
    ! Gauss-Jordan integrator
    integer, intent(in) :: n ! number of elements in b
    integer :: i, j ! iterator variables
    real, dimension(n,n), intent(inout) :: a
    real :: scale_factor
    real :: subtract_factor ! initialize a factor
    real, dimension(n,1), intent(inout) :: b

    integer, dimension(n) :: row, col ! pivot indices
    integer,target :: pivot(2) ! store max pivot index
    ! declare pointers that point to the elements of the pivot vector
    integer, pointer :: pivot_1,pivot_2
    ! initialize the pointers
    pivot_1 => pivot(1)
    pivot_2 => pivot(2)

    forall (i=1:n) row(i) = i
    forall (j=1:n) col(j) = j

    ! ========================================================================
    ! ==========
    ! ==========
    ! ========== BEGIN FORWARD ELIMINATION PHASE
    ! ==========
    ! ==========
    ! ========================================================================

    ! iterate through all of the columns save the last
    do j = 1,n-1

        ! find the pivot in the sub-matrix and store them in the pivot array
        pivot_1 = maxloc(abs(A(j:,j)),1) + (j-1);
        pivot_2 = j;

        ! if pivot_1 /= pivot_2 there needs to be a row swap to get the pivot
        ! on the diagonal
        if (pivot_1 /= pivot_2) then
            call swap(A(pivot_1,:),A(pivot_2,:))
            call swap(B(pivot_1,:),B(pivot_2,:))
        end if

    ! zero out the jth column entries that are under the diagonal
    do i = j+1,n
        scale_factor = (A(i,j)/A(j,j))
        A(i,:) = A(i,:) - scale_factor * A(j,:)
        B(i,1) = B(i,1) - scale_factor * B(j,1)
    end do

    end do ! end of iteration through column

    ! check here to make sure that the matrix is in upper triangular form
    ! by checking the lower triangular part of the matrix to make sure all
    ! the entries there are zero
    do i=1,n
        do j=i+1,n
            ! matrix is not in upper triangular form. Bail out.
            if (A(j,i) >= 1.1920929E-07) call abort
        end do
    end do
    ! ========================================================================
    ! ==========
    ! ==========
    ! ========== BEGIN BACK SUBSTITUTION PHASE
    ! ==========
    ! ==========
    ! ========================================================================

    do i = n,1,-1
        subtract_factor = 0.0
        do j = i+1,n
            subtract_factor = subtract_factor + A(i,j)*B(j,1)
        end do
        B(i,1) = (B(i,1) - subtract_factor)/A(i,i)
    end do

end subroutine gaussj

end module utilities
