module TypesAndDefs
  implicit none
  !Parameters
  integer, parameter :: pad = 4
  logical :: plan_made
  integer :: N
  integer :: N2
  integer :: r_end
  integer, parameter :: xrange = 2**8
  real, parameter :: w=2.0**(-3)
  real, parameter :: kappa = 2.0**(-3)
  real :: dt
  real :: dx
  complex, parameter :: i = (0.0,1.0)
  integer :: t_end
  integer :: ii, j_s
  real, allocatable, dimension(:) :: k
  real, allocatable, dimension(:) :: r
  real, parameter :: dtmin = 1.e-15
  real, parameter :: dtsafe = 1.e-7
  real, parameter :: eps = 1.e-7
  !Temporary fft stuff
  !integer*8:: plan

  real, allocatable, dimension(:) :: ones
  complex, allocatable, dimension(:,:) :: cap_B
  complex, dimension(3,3) :: identity

  complex, allocatable, dimension(:,:) :: Grad
  complex, allocatable, dimension(:) :: Grad2
  complex, dimension(3,3) :: E_mat
  
  interface operator (.x.)
     module procedure cross_product
  end interface operator (.x.)

  interface operator (.dot.)
     module procedure dot_prod_vectorized
  end interface operator (.dot.)

  interface operator (.dot.)
     module procedure dot_prod
  end interface operator (.dot.)
contains
  function cross_product(x,y) result(z)
    !Define the cross product between 2 vectors x(r) and y(r)
    implicit none
    complex,intent(in), dimension(:,:) :: x, y
    complex, dimension(size(x,1),size(x,2)) ::z
    if(size(x,1) /= size(y,1))then
       write(*,*) 'Error: vectors not the same length'
       stop
    end if
    !z_x = x_y*y_z - x_z*y_y
    !z_y = x_z*y_x - x_x*y_z
    !z_z = x_x*y_y - x_y*y_x
    z(:,1) = x(:,2)*y(:,3) - x(:,3)*y(:,2)
    z(:,2) = x(:,3)*y(:,1) - x(:,1)*y(:,3)
    z(:,3) = x(:,1)*y(:,2) - x(:,2)*y(:,1)
    return
  end function cross_product

  function dot_prod_vectorized(x,y) result(z)
    implicit none
    complex, dimension(:,:), intent(in) :: x, y
    complex :: z(size(x,1))
    z = x(:,1)*y(:,1) + x(:,2)*y(:,2) + x(:,3)*y(:,3)
  end function dot_prod_vectorized

  function dot_prod(x,y) result(z)
    implicit none
    complex, dimension(:), intent(in) :: x, y
    complex :: z
    integer :: j
    do j=1,size(x)
       z = z + x(j)*y(j)
    end do
  end function dot_prod

  function outer_product(x,y) result(z)
    implicit none
    !For vectors
    complex, dimension(:), intent(in) :: x, y
    complex, dimension(size(x),size(y)) :: z
    integer :: i1,j1
    
    do j1=1,3
       do i1=1,3
          z(j1,i1) = x(i1)*y(j1)
       end do
    end do
  end function outer_product
end module TypesAndDefs

