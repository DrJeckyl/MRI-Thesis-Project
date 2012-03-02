Program Main
  use TypesAndDefs
  use Funcs
  use nr, ONLY: EquationsOfMotion,RungeKutta4,omp_set_num_threads,omp_set_dynamic
  implicit none
  include "fftw3.f"

  !Buffer variables to read in dt and dx
  integer :: arg1, arg2
  character *100 buffer

  !PGPLOT Variables
  integer :: PGOPEN
  real(kind=4), allocatable, dimension(:) :: varx   !Single Precision to make PGPLOT happy
  real(kind=4), allocatable, dimension(:) :: vary   !
  real(kind=4), allocatable, dimension(:) :: vart   !
  real(kind=4), allocatable, dimension(:) :: varz   !
  real(kind=4), allocatable, dimension(:) :: vark   !
  real(kind=4), allocatable, dimension(:,:) :: ampu !Variable to store log of the amplitude, single precision

  !Solution Variables
  complex, allocatable, dimension(:,:) :: u !Real Space
  complex, allocatable, dimension(:,:) :: u_!K-space
  complex, allocatable, dimension(:) :: div_b, div_u !divergence checking
  complex, allocatable, dimension(:,:) :: u_ic !Initial conditions
  complex, allocatable, dimension(:,:) :: u_ic_
  complex, allocatable, dimension(:,:) :: u_dumb 

  !Quadradic Quantities and dependancy variables
  complex, allocatable, dimension(:,:) :: gradQ  !gradient of the pressure
  complex, allocatable, dimension(:) :: phi !scalar potential
  complex, allocatable, dimension(:) :: phi_
  complex, allocatable, dimension(:) :: del_phi !the other scalar potential
  complex, allocatable, dimension(:) :: del_phi_
  complex, allocatable, dimension(:,:) :: jc  !current density
  complex, allocatable, dimension(:,:) :: a !Vector Potentials
  complex, allocatable, dimension(:,:) :: a_
  !complex, allocatable, dimension(:) :: Grad2 !Laplace operator in k-space
  complex, allocatable, dimension(:,:) :: jh !magnetic helicity
  complex, allocatable, dimension(:,:) :: jh_
  complex, allocatable, dimension(:,:) :: L !Angular Momentum
  complex, allocatable, dimension(:,:) :: L_
  real, allocatable, dimension(:) :: Jh_flux
  real, allocatable, dimension(:) :: L_flux
  real, allocatable, dimension(:) :: dummy
  !complex, allocatable, dimension(:,:) :: Grad !Gradient in k-space
  complex, allocatable, dimension(:,:) :: Grad_temp
  complex, allocatable, dimension(:) :: phi_source !Source term for the laplace equation
  
  !Time and space variables
  integer :: t
  real :: sigma  !Initial condition variables
  real :: ag
  real :: t_current !Actual time
  real, allocatable, dimension(:) :: kt        !Space grid--->k + t

  !Wall time calculations
  integer, dimension(8) :: wall_time_start, wall_time_end, time_diff
  character*10 :: string_time_start, string_time_end

  !Code variables, counters etc.
  integer :: j
  logical :: FFTIC
  integer :: omp_get_num_threads, omp_get_num_procs, omp_get_max_threads
  complex, dimension(3,3) :: dummy_matrix
  complex(kind=8), allocatable, dimension(:) :: k_c
!!!!!!!!!!!!
  ! Program Block
!!!!!!!!!!

  complex(kind=8) :: array(11)


  !Get the arguments from the main program
  call getarg(1,buffer)
  read(buffer,*) arg1
  call getarg(2,buffer)
  read(buffer,*) arg2

  N = arg1   !Sets dx
  N2 = arg2  !Sets dt
  
  dx = xrange/2.0**(N)  !(x-N), x sets the domain to 0..2**x
  dt = 2.0**(-N2)
  
  r_end = 2**(N)  !0 to r_end*dx
  t_end = int(16./dt)  !int(x/dt) -- x is the time to integrate up to
  
  FFTIC = .false.
  FFTIC = .true.

  !Plan is defaulted to false
  plan_made = .false.

  !3x3 identity matrix
  identity = 0.
  identity(1,1) = 1.
  identity(2,2) = 1.
  identity(3,3) = 1.

  !The matrix operator to pick out the E12 and E21 components of u
  !Corresponds to the E matrix from the WC method
  E_mat = 0.
  E_mat(1,2) = 2./w
  E_mat(2,1) = -.5/w


  !Allocate the arrays
  !Allocate all the arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(k(r_end),r(r_end))
  allocate(kt(r_end))
  allocate(ones(r_end))
  allocate(cap_B(r_end,3))
  allocate(u(r_end,6))
  allocate(u_(r_end,6))
  allocate(u_ic(r_end,6))
  allocate(u_ic_(r_end,6))
  allocate(u_dumb(r_end,3))
  allocate(ampu(6,t_end))
  allocate(gradQ(r_end,3))
  allocate(phi(r_end))
  allocate(del_phi(r_end))
  allocate(phi_(r_end))
  allocate(del_phi_(r_end))
  allocate(jc(r_end,3))
  allocate(a(r_end,3))
  allocate(a_(r_end,3))
  allocate(Grad2(r_end))
  allocate(jh(r_end,3))
  allocate(jh_(r_end,3))
  allocate(L(r_end,3))
  allocate(L_(r_end,3))
  allocate(Jh_flux(t_end))
  allocate(L_flux(t_end))
  allocate(dummy(r_end))
  allocate(vart(t_end))
  allocate(varz(t_end))
  allocate(varx(r_end))
  allocate(vary(r_end))
  allocate(Grad(r_end,3))
  allocate(phi_source(r_end))
  allocate(div_u(r_end),div_b(r_end))
  allocate(Grad_temp(r_end,3))
  allocate(k_c(r_end))

  !Set the grid up
  do j=1,r_end
     r(j) = real((j-1)*dx)
  end do
  
  !do j=1,r_end
  !   k(j) = real((j-1)*dx - xrange/2.)
  !end do

  !k = r - xrange/2.

  !k_c = k
  !call fft_shift(k_c)
  !k = real(k_c)

  do j=1,r_end/2
     k(j) = real((j-1)/(r_end*1.)/2./dx)
     k(r_end/2 + j) = real((-r_end + j - 1 + r_end/2.)/2./(r_end*1.)/dx)
  end do    
  
!!$  print *, 'k:'
!!$  print *, k(1:5)
!!$  print *, k(r_end/2-2:r_end/2+2)
!!$  print *, k(r_end-4:r_end)
!!$  print *, 'r:'
!!$  print *, r(1:5)
!!$  print *, r(r_end/2-2:r_end/2+2)
!!$  print *, r(r_end-4:r_end)
!!$  stop

  !Define some module variables to use later
  ones = 1.0
  cap_B(:,1) = 0.
  cap_B(:,2) = 1.
  cap_B(:,3) = 0.

  !Print out the parameters of the run
  print *, 'dt:',dt
  print *, 'dx:',dx
  print *, 'N:',N
  print *, 'w_A',w
  print *, 'kappa',kappa
  print *, 'k ranges from:',k(1),'to:',k(r_end)

!!$
!!$  array(1:5) = 1.0
!!$  array(6:11) = 0.0
!!$
!!$  print *, 'Testing shift routine'
!!$  print *, 'array before shift'
!!$  print *, array
!!$  print *, 'array after shift'
!!$  call fft_shift(array)
!!$  print *, array
!!$  print *, 'array after ishift'
!!$  call ifft_shift(array)
!!$  print *, array
!!$
!!$  stop


  !Set up initial conditions
  sigma = 10.0/2.0/sqrt(2.0*log(2.0))
  ag = 1.0/2.0/(sigma)**2

  !u is the solution vector containing both u and b
  !u(:,1:3) = bx,by,bz
  !u(:,4:6) = ux,uy,uz
  u=0.0
  u_=0.0

  if(FFTIC)then
     u(:,1) = exp(-ag*(r-xrange/2.)**2)
     !u(:,2) = -3.*ag*r/i/w*u(:,1)
     u(:,2) = -1.5*xder(u(:,1))/i/w
  else
     !Do an analytic Fourier Transform of the initial conditions
     u(:,1) = exp(-1./ag*k**2)
     u(:,2) = -1.5*k/w*u(:,1)
  end if

  !Make a snapshot of the initial conditions to plot later
  u_ic = u 


!!$!Test the FT
!!$  print *, 'bx'
!!$  print *, u(1:10,1)
!!$  print *, 'by'
!!$  print *, u(1:10,2)

  if(FFTIC)then
     !FFT the initial conditions and centre them
     !call fftw_init_threads
     
     !call date_and_time(values=wall_time_start)
     
     !Use the new fft routines I wrote
     !call dfftw_plan_dft_1d(plan,r_end,u(:,1),u_(:,1),FFTW_FORWARD,FFTW_ESTIMATE)
     call fft(u(:,1),u_(:,1))
     call fft(u(:,2),u_(:,2))
     !call dfftw_destroy_plan(plan)
     !plan_made = .false.


!!$  !Test the FT
!!$  print *, 'bx'
!!$  print *, u_(1:10,1)
!!$  print *, 'by'
!!$  print *, u_(1:10,2)
!!$  stop
  else
     u_ = u
  end if
  

  !Grab a snapshot of the FT'ed initial conditions
  u_ic_ = u_

  !Start Timing the code
  call date_and_time(values=wall_time_start)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*********************** MAIN LOOP ******************************************
  !Solving the equations
  print *, 'Using Runge Kutta 4 routine'
  print *, 'Start Time:'
  time_diff = wall_time_start
  write(*,"(F8.5,a,I2.2,a,I2.2,a,I2.2,a,I3.3,a)") real(100.0*t)/real(t_end),'% ',time_diff(5),'h:',time_diff(6),'m:',time_diff(7),'.',time_diff(8),'s'

  do t=1,t_end
     !Set the current time and space variables
     t_current = real((t-1)*dt)
     kt = k + t_current

     !Define the Gradient operator in Fourier dpace
     Grad(:,1) = 1.5*kt
     Grad(:,2) = w
     Grad(:,3) = w/kappa

     !Laplace operator in Fourier Dpace
     Grad2 = Grad .dot. conjg(Grad)

     !Print out code progress every 1000 time steps
     if(mod(t,1000)==0)then
        call date_and_time(values=wall_time_end)
        time_diff = t_diff(wall_time_start,wall_time_end)
        write(*,"(F8.5,a,I2.2,a,I2.2,a,I2.2,a,I3.3,a)") real(100.0*t)/real(t_end),'% ',time_diff(5),'h:',time_diff(6),'m:',time_diff(7),'.',time_diff(8),'s'
     endif

     !Solving the equations at every space step
     do j_s=1,r_end
        if(sum(abs(u_(j_s,:))) .ne. 0.)then
           call RungeKutta4(t_current,u_(j_s,:),dt,EquationsOfMotion)
        end if
     end do

     do j=1,r_end
        if(abs(Grad2(j)) > 1.*dx)then
           dummy_matrix = outer_product(Grad(j,:),Conjg(Grad(j,:)))
           dummy_matrix = identity-dummy_matrix/Grad2(j)
        end if
        u_(j,1:3) = MATMUL(dummy_matrix,u_(j,1:3)) 
     end do

     !Calculate the vector potential from \nabla^2 a = -j, jc = \nabla x b
     ! \nabla^2 a = - (\nabla x b)
     ! a = -(\nabla x b) * \nabla^-2
     jc = (Grad .x. u_(:,1:3))
     a = -jc/Grad2(1)

     !Calulate the scalar potentials phi, del_phi
     !phi = \nabla^-2 \nabla . (u x B)
     phi_ = (Grad .dot. (u_(:,4:6) .x. cap_B))/Grad2(1)

     !Calculate b(del_phi - a.V)
     ! b(del_phi - a.V) = -2b \int dx'/4pi * G(x,x')(d_x j_i)S_ik
     ! \nabla^4 b(del_phi - a.V) = -2b * (d_x j_i)S_ik
     ! \nabla^4 (del_phi - a.V) = -2 * (d_x j_i) * S_ik
     ! The code uses del_phi to represent (del_phi - a.V)
     ! del_phi = -2 * \nabla^-4 * (d_x j_i) * S_ik
     ! del_phi = -2 * \nabla^-4 * (-9./8. * \partial_x j_y - 3./4.*i*w*j_x)
     !                                    * 3./2. * i * kt * j_y - 3./4. * i * w * j_x
     del_phi_ = -2. / Grad2(1)**2 * (-9./8. * 1.5 * i * kt * jc(:,2) - 0.75 * i * w * jc(:,1))

     !Put everything together to find the magnetic helicity flux
     !Calculate jh as (a.B)u + b\phi - b(a.V + b\varphi)
     do j=1,3
        jh(:,j) =  u_(:,j+3)*a(:,2) + u_(:,j)*phi_ + u_(:,j)*del_phi_
     end do
     
     !calculate divb and divu in fourier space
     div_b = Grad.dot.u_(:,1:3)!i*1.5*k*u_(:,1) + i*w*u_(:,2) + i*w/kappa*u_(:,3)
     div_u = Grad.dot.u_(:,4:6)!i*1.5*k*u_(:,4) + i*w*u_(:,5) + i*w/kappa*u_(:,6)

     !Shift and inverse transform everything
     !call date_and_time(values=wall_time_start)

     call ifft(u_(:,1),u(:,1))  
     call ifft(u_(:,2),u(:,2))
     call ifft(u_(:,3),u(:,3))
     call ifft(u_(:,4),u(:,4))
     call ifft(u_(:,5),u(:,5))
     call ifft(u_(:,6),u(:,6))

     do j=1,3
        call ifft(jh(:,j),jh(:,j))
        call ifft(div_b,div_b)
        call ifft(div_u,div_u)
     end do

!!$     call ifft_shift(div_b)
!!$     call ifft_shift(div_u)
!!$     call ifft_shift(u_(:,1))
!!$     call ifft_shift(u_(:,2))
!!$     call ifft_shift(u_(:,3))
!!$     call ifft_shift(u_(:,4))
!!$     call ifft_shift(u_(:,5))
!!$     call ifft_shift(u_(:,6))
!!$     
     !call dfftw_destroy_plan(plan)
     !plan_made = .false.

     !Calculate divb and divu
     !div_b = 1.5*xder(u(:,1)) + i*w*u(:,2) + i*w/kappa*u(:,3)
     !div_u = 1.5*xder(u(:,4)) + i*w*u(:,5) + i*w/kappa*u(:,6)

     !Find the maximum of the variables to log versus the time
     !ampu(1,t) = real(log(maxval(abs(u_(:,1)))+1.0E-30),kind=4)
     ampu(1,t) = real(log(maxval(abs(u(:,1)))+1.0E-30),kind=4)
     ampu(2,t) = real(log(maxval(abs(u(:,2)))+1.0E-30),kind=4)
     ampu(3,t) = real(log(maxval(abs(u(:,3)))+1.0E-30),kind=4)
     ampu(4,t) = real(log(maxval(abs(u(:,4)))+1.0E-30),kind=4)
     ampu(5,t) = real(log(maxval(abs(u(:,5)))+1.0E-30),kind=4)
     ampu(6,t) = real(log(maxval(abs(u(:,6)))+1.0E-30),kind=4)

     !Romberg integration routine, syntax: Romberg(input_array, lbound, ubound)
     L_flux(t) = real(Romberg((conjg(u(:,4))*u(:,5)-conjg(u(:,1))*u(:,2)),0.,real(xrange)))

     !Jh_flux(t) = real(Romberg(sqrt(conjg(jh(:,1))*jh(:,2)),-64.,64.))
     Jh_flux(t) = real(Romberg(sqrt(conjg(jh(:,3))*jh(:,3)),0.,real(xrange)))
  end do


     !Book keeping for the PGPLOT routines(required)
     call PGSLCT(PGOPEN('/ps'))
     call PGASK(.false.)
     !if(PGOPEN('?') < 1)stop
     do j=1,t_end
        vart(j) = real((j-1)*dt,kind=4)
     end do

     varx = real(r,kind=4)
     vark = real(r - xrange/2.,kind=4)


     !Plot the solution variables
     !syntax is plot(x,y,'xaxis','yaxis','title')

     vary = real(u(:,1),kind=4)
     call plot(varx,vary,'x','bx','bx vs x')

     vary = real(u(:,2),kind=4)
     call plot(varx,vary,'x','by','by vs x')

     vary = real(u(:,3),kind=4)
     call plot(varx,vary,'x','bz','bz vs x')

     vary = real(u(:,4),kind=4)
     call plot(varx,vary,'x','ux','ux vs x')

     varz = real(ampu(1,:),kind=4)
     call plot(vart,varz,'t','|bx(t)|','amplitude of bx vs. time')

     vary = real(abs(u_ic(:,1)),kind=4)
     call plot(varx,vary,'x','b','bx IC')
  
     vary = real(abs(u_ic_(:,1)),kind=4)
     call plot(varx,vary,'x','b','bx IC FTed')

     vary = real(jh(:,1),kind=4)
     call plot(varx,vary,'t','jh','jh(x)')

     varz = real(Jh_flux,kind=4)
     call plot(vart,varz,'t','Jh','Jh(t)')

     varz = real(L_flux,kind=4)
     call plot(vart,varz,'t','L_flux','L(t)')

     varz = log(abs(Jh_flux/L_flux))
     call plot(vart,varz,'t','log(J/L)','Raito of J to L')

     vary = real(abs(div_b),kind=4)
     call plot(varx,vary,'x','divb','div b vs x')

     vary = real(abs(div_u),kind=4)
     call plot(varx,vary,'x','divu','div u vs x')

     call PGCLOS

     !Write some files out to work with later, binary format
     
     open(60,file='bx.out',status='unknown',form='formatted')
     do j=1,t_end
        write(60,*) ampu(1,j)
     end do
     close(60)

     open(60,file='J.out',status='unknown',form='formatted')
     do j=1,t_end
        write(60,*) Jh_flux(j)
     end do

     close(60)

     open(60,file='L.out',status='unknown',form='formatted')
     do j=1,t_end
        write(60,*) L_flux(j)
     end do
     close(60)

   End Program Main
