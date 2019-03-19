program question4
use utilities
implicit none

! number of points in the data set
integer, parameter :: n = 9130

! number of basis functions
integer, parameter :: m = 3

real :: lambda

real :: x_var, y_var
real :: xmin, xmax
real :: criterion
integer, parameter :: num_grid_fine = 10000
real, dimension(num_grid_fine) :: fine_grid
real, dimension(num_grid_fine) :: fit_vals
real, dimension(n) :: x
real, dimension(n) :: y
real, dimension(m+1) :: c
real, dimension(m+1) :: delta
real, dimension(n, m) :: basis_matrix
real, dimension(n, m+1) :: J
real, dimension(m+1, n) :: JT
real, dimension(m+1, m+1) :: JTJ
real, dimension(m+1, m+1) :: JTJ_diagonal

real :: chi_sq_old, chi_sq_new
real :: epsilon
integer :: i
integer :: jj
integer :: index1, iteration
integer :: status

status = 0
index1=1

! read in data from std out
do while (status == 0)
    read (*,*,iostat=status) x_var, y_var
    if (status < 0) exit

    x(index1) = x_var; y(index1) = y_var

    index1 = index1 + 1
end do

xmin = minval(x)
xmax = maxval(x)

epsilon = 0.000001

! initial guess
c = (/0.5, 0.75, 3.2, 0.0/)

lambda = 0.0001

basis_matrix = basis(n, m, x)
iteration=1

! implement Levenberg-Marquardt step
do
    J = Jacobian(n, m, c, basis_matrix)
    JT = transpose(J)
    JTJ = matmul(JT, J)
    JTJ_diagonal(:,:) = 0.0
    do jj=1,m
        JTJ_diagonal(jj,jj) = JTJ(jj,jj)
    end do

    delta = matmul(inv(JTJ + lambda*JTJ_diagonal), matmul(JT, y-f(n,m,c,basis_matrix)))

    chi_sq_old = chi_sq(n, m, y, c, basis_matrix)
    chi_sq_new = chi_sq(n, m, y, c + delta, basis_matrix)

    write (*,*) "Iteration:", iteration, " Chi-Square(t):",chi_sq_old,"Chi_square(t+1):",chi_sq_new

    criterion = chi_sq_new - chi_sq_old

    ! check convergence criterion
    if (abs(criterion/chi_sq_old) < epsilon) exit

    ! update the Marquardt parameter
    if (abs(criterion) >= 0.0) then
        lambda = 10.0*lambda
    else
        lambda = lambda/10.0
    end if

    c = c + delta

    iteration = iteration + 1

end do

write (*,*) "Final coefficients:", c
write (*,*) "Final Chi-sq value:", chi_sq(n, m, y, c, basis_matrix)

call linspace(fine_grid, xmin, xmax, num_grid_fine)

call dump_expression(m, xmin, xmax, num_grid_fine, c, fit_vals)

open(unit=1,file="q5_output.dat")
    do i=1,num_grid_fine
        write(1,*) fine_grid(i), fit_vals(i)
    end do
close(1)

contains

! evaluate the chi_squared function
pure function chi_sq(n, m, y, c, basis_matrix)
    integer, intent(in) :: n, m
    real, dimension(n), intent(in) :: y
    real, dimension(m), intent(in) :: c
    real, dimension(n,m), intent(in) :: basis_matrix
    real :: chi_sq

    chi_sq = sum((y-f(n,m,c,basis_matrix))**2)
end function chi_sq

! fit function
pure function f(n, m, c, basis_matrix)
    integer, intent(in) :: n, m
    real, dimension(m+1), intent(in) :: c
    real, dimension(n,m), intent(in) :: basis_matrix
    real, dimension(n) :: f

    f = exp(matmul(basis_matrix, c(:m))) + c(m+1)
end function f

! derivative of a function to find a minimum of
function Jacobian(n, m, c, basis_matrix)
    integer, intent(in) :: n, m
    real, dimension(n, m), intent(in) :: basis_matrix
    real, dimension(m), intent(in) :: c
    real, dimension(n, m+1) :: Jacobian
    real, dimension(n) :: temp

    integer :: j

    do j=1,m
        temp = matmul(basis_matrix, c)
        ! write (*,*) temp
        Jacobian(:,j) = basis_matrix(:,j)*exp(matmul(basis_matrix,c))
    end do
    Jacobian(:,m+1) = 1.0
end function Jacobian

! Returns the inverse of a matrix calculated by finding the LU
! decomposition
function inv(A) result(Ainv)
    real, dimension(:,:), intent(in) :: A
    real, dimension(size(A,1),size(A,2)) :: Ainv

    real, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
        stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
        stop 'Matrix inversion failed!'
    end if
end function inv

! construct basis matrix
pure function basis(n, m, x)
    integer, intent(in) :: n, m
    real, dimension(n), intent(in) :: x
    real, dimension(n, m) :: basis
    integer :: i, j

    do j=1,m
        do i=1,n
            basis(i,j) = cos(2*PI*(j-1)*x(i))
        end do
    end do
end function basis

! basis functions we are fitting
subroutine evalb(x, n, B)
    real, intent(in) :: x
    integer, intent(in) :: n
    real, dimension(n), intent(inout) :: B
    integer :: i

    forall (i=1:n) B(i) = cos(2*PI*x*(i-1))
end subroutine evalb

subroutine dump_expression(m, x_min, x_max, num_grid_fine, U, fit_vals)
    integer, intent(in) :: m
    integer, intent(in) :: num_grid_fine
    real, intent(in) :: x_min, x_max
    real, dimension(num_grid_fine), intent(inout) :: fit_vals
    real, dimension(m+1) :: U
    real, dimension(num_grid_fine) :: x_pts
    real, dimension(m) :: basis
    integer :: i

    call linspace(x_pts, x_min, x_max, num_grid_fine)

    do i = 1,num_grid_fine
        basis = 0.0
        call evalb(x_pts(i), m, basis)
        fit_vals(i) = exp(sum(U(:m)*basis)) + U(m+1)
    end do
end subroutine dump_expression

end program question4