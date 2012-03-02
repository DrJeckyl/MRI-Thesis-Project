subroutine EquationsOfMotion(x,y,dydx)
  use TypesAndDefs
  implicit none
  real, intent(in) :: x
  complex, intent(in), dimension(:) :: y
  complex, intent(out), dimension(:) :: dydx
  
  !Algorithm Variables
  complex :: dummy
  complex, dimension(3) :: GradQ
  complex, dimension(3) :: Nabla
  complex :: Nabla2
  complex, dimension(3,3) :: Pi_k, Pi_k2
  integer :: j
  real :: xt

  !clear dydx before doing anything
  dydx = (0.,0.)

  xt = x + k(j_s)
  
  Nabla(1) = 1.5*xt
  Nabla(2) = w
  Nabla(3) = w/kappa
  Nabla2 = Nabla .dot. Conjg(Nabla)

!!$  dummy = 9./4.*xt**2 + w**2 + w**2/kappa**2
!!$
!!$  gradQ(2) = (3.*xt*y(5) - 2.*y(4)*w)/dummy
!!$  gradQ(1) = 1.5 * xt / w * gradQ(2)
!!$  gradQ(3) = gradQ(2)/kappa

  Pi_k = outer_product(Nabla,Conjg(Nabla))/Nabla2      
  Pi_k2 = outer_product(Nabla,[i*cmplx(1.5),cmplx(0.),cmplx(0.)])/Nabla2
  
  GradQ = MATMUL(Pi_k,MATMUL(E_mat,y(4:6))) - MATMUL(Pi_k2,y(4:6)) + MATMUL(Pi_k,i*y(1:3))

!!$  if(GradQ(1) .ne. GradQ2(1))then
!!$     print *, 'GradQ(1) not the same, GradQ(1)=',GradQ(1), GradQ2(1)
!!$     !stop
!!$  end if
!!$  if(GradQ(2) .ne. GradQ2(2))then
!!$     print *, 'GradQ(2) not the same, GradQ(2)=',GradQ(2), GradQ2(2)
!!$     !stop
!!$  end if
!!$  if(GradQ(3) .ne. GradQ2(3))then
!!$     print *, 'GradQ(3) not the same, GradQ(3)=',GradQ(3), GradQ2(3)
!!$     !stop
!!$  end if

  dydx(1) =                    + i*y(4)
  dydx(2) =    -1.5*y(1)/w + i*y(5)
  dydx(3) =                    + i*y(6)
  dydx(4) =     2.*y(5)/w  + i*y(1) -gradQ(1)
  dydx(5) =    -.5*y(4)/w  + i*y(2) -gradQ(2)
  dydx(6) =                    + i*y(3) -gradQ(3)
  
  return
  
end subroutine EquationsOfMotion
