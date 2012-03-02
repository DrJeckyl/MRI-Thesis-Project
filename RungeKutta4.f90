subroutine RungeKutta4(x,y,h,derivs)
  use TypesAndDefs
  implicit none

  complex, dimension(:), intent(inout) :: y
  real, intent(inout) :: x
  real, intent(in) :: h

  interface
     subroutine derivs(x,y,dydx)
       use TypesAndDefs
       implicit none
       real, intent(in) :: x
       complex, dimension(:), intent(in) :: y
       complex, dimension(:), intent(out) :: dydx
     end subroutine derivs
  end interface

  complex, dimension(size(y)) :: k1,k2,k3,k4
  complex, dimension(size(y)) :: dydx

  call derivs(x,y,dydx)
  k1 = h*dydx
  call derivs(x+.5*h,y+k1/2.,dydx)
  k2 = h*dydx
  call derivs(x+.5*h,y+k2/2.,dydx)
  k3 = h*dydx
  call derivs(x+h,y+k3,dydx)
  k4 = h*dydx

  y = y + k1/6. +k2/3. + k3/3. + k4/6.

end subroutine RungeKutta4
