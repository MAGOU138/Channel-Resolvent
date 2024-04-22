module global
   implicit none
   include 'mpif.h'

   ! mathematic constant
   real*8, parameter::pi = 4d0*datan(1d0)
   complex*16, parameter::IC = (0.d0, 1.d0)

   !=================================================================================!
   ! These are read from cfg file
   ! physical parameters
   real*8::Re_tau = 180.d0    ! Friction Reynolds number
   real*8::Karman = 0.426d0    ! Karman constant
   real*8::B_wall = 0.41d0    ! wall law constant
   real*8::B1_wall = 0.2d0   ! wall law constant

   ! discrete parameters
   integer::Nx = 2, Ny = 33, Nz = 2, Nt = 2
   real*8::dkx = 1.d0, dkz = 1.d0, domega = 2.d0*PI

   ! model choices
   logical::if_resolv_eddy = .false., if_resolv_sweep = .false.
   logical::if_force_composite = .false.
   character(len=16)::force_spectra_type = 'plain'
   character(len=8)::eddy_model_name = 'Cess'
   ! end vars read from cfg file
   !=================================================================================!

   ! variables shared
   ! discrete size info.
   integer::N, num_3var, num_4var

   real*8::Re                       ! Reynolds number based on half channel height and bulk velocity
   real*8::ratio_u_tau_2_U_bulk     ! Friction velocity

   ! Chebyshev
   real*8, allocatable, dimension(:)::cheby_points
   real*8, allocatable, dimension(:, :)::cheby_derivative_1, cheby_derivative_2
   complex*16, allocatable, dimension(:, :)::cheby_integrate
   real*8, allocatable, dimension(:)::integral_coeff

   ! mean velocity
   real*8, allocatable, dimension(:)::du_dy, u

   ! rms fluctutaing velocity
   real*8, allocatable, dimension(:)::m_uu, m_vv, m_ww, d_m_v_dy

   ! scale used for sweep
   real*8, allocatable, dimension(:, :, :)::lambda_kx_y_kz
   real*8, allocatable, dimension(:)::lambda

   ! effective viscosity
   real*8, allocatable, dimension(:)::nu_e, nu_t, d_nu_t_dy

   ! wave and frequency
   real*8, allocatable, dimension(:)::wave_x, wave_z, freque

   ! Resolvent matrix
   real*8, allocatable, dimension(:, :)::Resolvent_B
   complex*16, allocatable, dimension(:, :)::H, H_spatial
   complex*16, allocatable, dimension(:, :)::H_bc_1, H_bc_2
   real*8, allocatable, dimension(:, :)::Force_spectra
   complex*16, allocatable, dimension(:, :)::Force_spectra_comp
   complex*16, allocatable, dimension(:, :)::Resolvent_H, Resolvent_R, Resolvent_F, Resolvent_P
   complex*16, allocatable, dimension(:, :)::Resolvent_R_trans
   complex*16, allocatable, dimension(:, :)::Resolvent_BF
   complex*16, allocatable, dimension(:, :)::Oper_nonl, Oper_prod, Oper_pres, Oper_diss
   complex*16, allocatable, dimension(:, :)::Res_nonl, Res_prod, Res_pres, Res_diss

   ! tke
   real*8, allocatable, dimension(:)::tke_nonl, tke_prod, tke_pres, tke_diss

   ! output
   ! variables among all processes
   real*8, allocatable, dimension(:, :)::auto_corr_uu, auto_corr_vv, auto_corr_ww, auto_corr_pp
   complex*8, allocatable, dimension(:, :)::cros_corr_uu, cros_corr_vv, cros_corr_ww, cros_corr_pp

   ! mpi parallel
   integer::ierr, num_process, my_id
   integer::size_x, size_z
   integer, allocatable, dimension(:)::index_x, index_z

   ! Cess parameter
   real*8, parameter:: Cess_A = 25.4d0

contains
   subroutine allocate_var
      implicit none

      ! discrete size info.
      N = Ny - 1         ! degree of Chebyshev polynomials
      num_4var = 4*Ny
      num_3var = 3*Ny

      ! wave and frequency
      allocate (wave_x(0:Nx - 1), wave_z(0:Nz - 1), freque(0:Nt - 1))

      ! Chebyshev
      ! wall-normal discrete points-Chebyshev collocation, differential, and integral matrix
      allocate (cheby_points(0:N))  ! coordinates of Chebyshev extreme points
      allocate (cheby_derivative_1(0:N, 0:N), cheby_derivative_2(0:N, 0:N))   ! 1st, 2nd derivativ
      allocate (cheby_integrate(0:N, 0:N), integral_coeff(0:N))   ! integral matrix and integral coefficients
      allocate (du_dy(0:N), u(0:N))
      allocate (m_uu(0:N), m_vv(0:N), m_ww(0:N), d_m_v_dy(0:N))
      allocate (lambda_kx_y_kz(0:Nx - 1, 0:N, 0:Nz - 1))
      allocate (lambda(0:N))
      allocate (nu_e(0:N), nu_t(0:N), d_nu_t_dy(0:N))
      allocate (Resolvent_B(0:num_4var - 1, 0:num_3var - 1))
      allocate (H(0:num_4var - 1, 0:num_4var - 1), H_spatial(0:num_4var - 1, 0:num_4var - 1))
      allocate (H_bc_1(0:3, 0:num_4var - 1))
      allocate (H_bc_2(0:3, 0:num_4var - 1))
      allocate (Force_spectra(0:num_3var - 1, 0:num_3var - 1))
      allocate (Force_spectra_comp(0:num_3var - 1, 0:num_3var - 1))
      allocate (Resolvent_H(0:num_4var - 1, 0:num_4var - 1))
      allocate (Resolvent_R(0:num_4var - 1, 0:num_3var - 1))
      allocate (Resolvent_F(0:num_3var - 1, 0:num_4var - 1))
      allocate (Resolvent_P(0:num_4var - 1, 0:num_4var - 1))
      allocate (Resolvent_R_trans(0:num_3var - 1, 0:num_4var - 1))
      allocate (Resolvent_BF(0:num_4var - 1, 0:num_4var - 1))
      allocate (Oper_nonl(0:num_4var - 1, 0:num_4var - 1))
      allocate (Oper_prod(0:num_4var - 1, 0:num_4var - 1))
      allocate (Oper_pres(0:num_4var - 1, 0:num_4var - 1))
      allocate (Oper_diss(0:num_4var - 1, 0:num_4var - 1))
      allocate (Res_nonl(0:num_4var - 1, 0:num_4var - 1))
      allocate (Res_prod(0:num_4var - 1, 0:num_4var - 1))
      allocate (Res_pres(0:num_4var - 1, 0:num_4var - 1))
      allocate (Res_diss(0:num_4var - 1, 0:num_4var - 1))
      allocate (tke_nonl(0:num_4var - 1))
      allocate (tke_prod(0:num_4var - 1))
      allocate (tke_pres(0:num_4var - 1))
      allocate (tke_diss(0:num_4var - 1))
      allocate (auto_corr_uu(0:N, 0:Nt - 1))
      allocate (auto_corr_vv(0:N, 0:Nt - 1))
      allocate (auto_corr_ww(0:N, 0:Nt - 1))
      allocate (auto_corr_pp(0:N, 0:Nt - 1))
      allocate (cros_corr_uu(0:N, 0:N))
      allocate (cros_corr_vv(0:N, 0:N))
      allocate (cros_corr_ww(0:N, 0:N))
      allocate (cros_corr_pp(0:N, 0:N))

      Re = Re_tau*((-1.d0 + log(Re_tau))/Karman + B_wall + B1_wall)
      ratio_u_tau_2_U_bulk = Re_tau/Re

   end subroutine allocate_var

   subroutine deallocate_var
      implicit none

      deallocate (wave_x, wave_z, freque)

      deallocate (cheby_points)
      deallocate (cheby_derivative_1, cheby_derivative_2)
      deallocate (cheby_integrate, integral_coeff)
      deallocate (du_dy, u)
      deallocate (m_uu, m_vv, m_ww, d_m_v_dy)
      deallocate (lambda_kx_y_kz, lambda)
      deallocate (nu_e, nu_t, d_nu_t_dy)
      deallocate (Resolvent_B)
      deallocate (H, H_spatial)
      deallocate (Resolvent_H)
      deallocate (Resolvent_R)
      deallocate (Resolvent_R_trans)
      deallocate (Resolvent_F, Resolvent_BF)
      deallocate (Resolvent_P)
      deallocate (Res_nonl, Res_prod, Res_pres, Res_diss)
      deallocate (H_bc_1)
      deallocate (H_bc_2)
      deallocate (Oper_nonl, Oper_prod, Oper_pres, Oper_diss)
      deallocate (tke_nonl, tke_prod, tke_pres, tke_diss)
      deallocate (Force_spectra, Force_spectra_comp)
      deallocate (auto_corr_uu, auto_corr_vv, auto_corr_ww, auto_corr_pp)
      deallocate (cros_corr_uu, cros_corr_vv, cros_corr_ww, cros_corr_pp)

      deallocate (index_x, index_z)

   end subroutine deallocate_var

end module global

!=========================================================================================!
subroutine parallel_init
   use global
   implicit none
   call mpi_init(ierr)
   call mpi_comm_size(mpi_comm_world, num_process, ierr)
   call mpi_comm_rank(mpi_comm_world, my_id, ierr)
   if (my_id == 0) write (*, *) 'Initializ MPI and get number of process'
   if (my_id == 0) write (*, *) 'Total process', num_process, 'are ready'
end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine load_cfg
   use global
   implicit none
   real*8::Nomega, f_temp
   integer::file_exist, i
   character(len=4)::resolv_eddy, resolv_sweep, force_composite
   open (20, file='set-up.ini')
   read (20, *)
   read (20, *) Re_tau              !  Friction Reynolds number
   read (20, *)
   read (20, *) Karman              !  Karman constant
   read (20, *)
   read (20, *) B_wall, B1_wall     !  log-wall constants
   read (20, *)
   read (20, *) Nx, Nz, Ny, Nt      !  wave numbers, Chebyshev point numbers, discrete frequency numbers
   read (20, *)
   read (20, *) dkx, dkz, Nomega    !  Resolution
   read (20, *)
   read (20, *) resolv_eddy         !  if eddy-enhanced resolvent
   read (20, *)
   read (20, *) eddy_model_name     !  which eddy model
   read (20, *)
   read (20, *) resolv_sweep        !  if sweep-enhanced resolvent
   read (20, *)
   read (20, *) force_spectra_type  !  which force spectra type
   read (20, *)
   read (20, *) force_composite     !  if force spectra composite
   close (20)
   ! assignment variables
   domega = 2*PI/Nomega
   if (trim(resolv_eddy) == 'true') if_resolv_eddy = .true.
   if (trim(resolv_sweep) == 'true') if_resolv_sweep = .true.
   if (trim(force_composite) == 'true') if_force_composite = .true.

   ! check input
   if (my_id == 0) then
      write (*, '(A, F8.3)') 'Friction Reynolds number, ', Re_tau
      write (*, '(A,3F8.3)') 'Wall constants, Karman, B, B1', Karman, B_wall, B1_wall
      write (*, '(A,4I8)') 'Discret size in x, z, y, t', Nx, Nz, Ny, Nt
      write (*, '(A,4F8.3)') 'Resolution dkx, dkz, domega', dkx, dkz, domega
      write (*, '(A,L2)') 'If resolvent eddy-enhanced: ', if_resolv_eddy
      write (*, '(A,L2)') 'If resolvent sweep-enhanced: ', if_resolv_sweep
      write (*, '(A)') 'Which eddy-enhanced model: '//eddy_model_name
      write (*, '(A)') 'Force spectra type: '//force_spectra_type
   end if

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine init_spectra
   use global
   implicit none
   integer::i
   wave_x = [(dkx*DBLE(i), i=0, Nx - 1)]
   wave_z = [(dkz*DBLE(i), i=0, Nz - 1)]
   freque = [(domega*DBLE(i - Nt/2), i=0, Nt - 1)]
end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine Chebyshev
   use global
   implicit none
   complex*16, allocatable, dimension(:, :)::cheby_derivative
   integer::i, j, status
   real*8::c

   do i = 0, N
      cheby_points(i) = cos(i*PI/N)
   end do

   do i = 0, N
      do j = 0, N
         ! i/=j
         if (i /= j) then
            if (mod(i, N) == 0 .and. mod(j, N) /= 0) then
               c = 2.d0
            else if (mod(i, N) /= 0 .and. mod(j, N) == 0) then
               c = 0.5d0
            else
               c = 1.d0
            end if
            cheby_derivative_1(i, j) = c*(-1)**(i + j)/(cheby_points(i) - cheby_points(j))
         end if
         ! i==j
         if (i == j) then
            cheby_derivative_1(i, j) = -cheby_points(i)/(2.d0*(1.d0 - cheby_points(i)**2))
         end if
         ! first and last diagonal
         cheby_derivative_1(0, 0) = (2.d0*N**2 + 1.d0)/6.d0
         cheby_derivative_1(N, N) = -cheby_derivative_1(0, 0)
      end do
   end do

   cheby_derivative_2 = matmul(cheby_derivative_1, cheby_derivative_1)

   allocate (cheby_derivative(0:N, 0:N), stat=status)
   cheby_derivative = cheby_derivative_1
   cheby_derivative(0, :) = 0.d0
   cheby_derivative(:, 0) = 0.d0
   cheby_derivative(0, 0) = 1.d0
   call INV_MAT(N + 1, cheby_derivative, cheby_integrate)
   integral_coeff = -real(cheby_integrate(N, :))
   deallocate (cheby_derivative)

   ! write out Cheby integral coefficients
   block
      character(len=128) ::filename
      filename = 'set-up/integral-coefficients.dat'
      open (10, file=filename, status='replace', action='write')
      do i = 0, N
         write (10, '(E24.12)') integral_coeff(i)
      end do
      close (10)
   end block

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine mean_velocity
   use global
   implicit none

   character(len=100)::fileName
   integer::i
   real*8::y_plus, y

   ! use Cess model to calculate total eddy and mean profile
   do i = 0, N
      y = cheby_points(i)
      y_plus = Re_tau*(1.d0 - abs(y))
      nu_e(i) = 0.5d0*sqrt(1.d0 + (Karman*Re_tau/3.d0*(1 + y**2 - 2.d0*y**4)*(1.d0 - exp(-y_plus/Cess_A)))**2) - 0.5d0
      nu_t(i) = nu_e(i) + 1.d0
   end do

   ! get mean shear from Cess viscosity
   do i = 0, N
      du_dy(i) = -(Re_tau)*cheby_points(i)/nu_t(i)
   end do
   ! get mean velocity profile by integration
   u = real(matmul(cheby_integrate, du_dy))
   u(0) = 0.d0

   if (my_id == 0) then

      fileName = 'set-up/mean-velocity.dat'
      open (20, file=fileName, status='replace', action='write')
      write (20, *) 'y yplus U U_plus'
      do i = 0, N
         y = cheby_points(i)
         y_plus = Re_tau*(1.d0 - abs(y))
         write (20, '(4E16.6)') y, y_plus, u(i)*ratio_u_tau_2_U_bulk, u(i)
      end do
      close (20)

      fileName = 'set-up/wave_x.dat'
      open (20, file=fileName, status='replace', action='write')
      write (20, *) 'wave_x'
      do i = 0, Nx - 1
         write (20, '(E16.6)') wave_x(i)
      end do
      close (20)

      fileName = 'set-up/wave_z.dat'
      open (20, file=fileName, status='replace', action='write')
      write (20, *) 'wave_z'
      do i = 0, Nz - 1
         write (20, '(E16.6)') wave_z(i)
      end do
      close (20)

      fileName = 'set-up/frequency.dat'
      open (20, file=fileName, status='replace', action='write')
      write (20, *) 'frequency'
      do i = 0, Nt - 1
         write (20, '(E16.6)') freque(i)
      end do
      close (20)

      fileName = 'set-up/Reynolds.dat'
      open (20, file=fileName, status='replace', action='write')
      write (20, *) 'Re_tau,    Re,    u_tao_2_u_bulk'
      write (20, '(3E16.6)') Re_tau, Re, ratio_u_tau_2_U_bulk
      close (20)

   end if

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine load_fluc_velocity
   use global
   implicit none
   integer::i, j, k, file_exist
   real*8::f_temp
   character(len=124)::filename

   m_uu = 0.d0
   m_vv = 0.d0
   m_ww = 0.d0
   d_m_v_dy = 0.d0
   if (my_id == 0) then
      ! processor id 0 reads mean square velocity
      open (3, file='./inputs/DNS-rms.dat', status='old', action='read', iostat=file_exist)
      if (file_exist /= 0) then
         write (*, *) 'NO DNS-rms.dat file provided'
         if (if_resolv_sweep) then
            stop
         else
            return
         end if
      end if
      read (3, *)! header of file
      do i = 0, N
         read (3, *) f_temp, m_uu(i), m_vv(i), m_ww(i)
      end do
      close (3)
      write (*, *) 'read DNS mean square of fluctuating velocity'

      ! processor id 0 reads scale file
      write (filename, '(A,I2,A,I2,A,I3,A)') 'lambda_kx_kz_y_', Nx, '_', Nz, '_', Ny, '.dat'
      open (4, file=filename, status='old', action='read', iostat=file_exist)
      if (file_exist /= 0) then
         write (*, *) 'NO scale file provided'
         if (if_resolv_sweep) then
            stop
         else
            return
         end if
      end if
      read (4, *)
      read (4, *)
      do i = 0, Nx - 1
         do k = 0, Nz - 1
            do j = 0, N
               read (4, *) f_temp, f_temp, f_temp, lambda_kx_y_kz(i, j, k)
            end do
         end do
      end do
      close (4)
      write (*, *) 'read scale file lamda finished'

      d_m_v_dy = matmul(cheby_derivative_1, sqrt(m_vv))
      ! broadcast to all other processors
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(m_uu, N + 1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(m_vv, N + 1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(m_ww, N + 1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(d_m_v_dy, N + 1, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(lambda_kx_y_kz, Nx*Ny*Nz, mpi_double_precision, 0, mpi_comm_world, ierr)
      call mpi_barrier(mpi_comm_world, ierr)
   end if

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine assemble_B
   use global
   implicit none
   integer::i
   Resolvent_B = 0.d0
   ! leading and trailing block (4*3) are zero
   ! this is the no-slip boundary condition
   ! Resolvent_B is a block-diagonal matrix
   do i = 1, N - 1
      Resolvent_B(4*i, 3*i) = 1.d0
      Resolvent_B(4*i + 1, 3*i + 1) = 1.d0
      Resolvent_B(4*i + 2, 3*i + 2) = 1.d0
   end do
end subroutine assemble_B
!=========================================================================================!

!=========================================================================================!
subroutine parallel_decompose
   use global
   implicit none
   integer::i, t1, t2
   integer::num_process_x, num_process_z

   if (mod(Nx*Nz, num_process) /= 0) then
      write (*, *) 'Error! number of waves must be divided evenly by number of processes'
      stop
   end if
   if (Nx*Nz < num_process) then
      write (*, *) 'Error! number of waves must larger than number of processes'
      stop
   end if

   num_process_x = 1
   do i = 1, max(num_process, Nx)
      if (mod(Nx, i) == 0 .and. mod(num_process, i) == 0) then
         num_process_x = max(num_process_x, i)
      end if
   end do
   num_process_z = num_process/num_process_x

   t1 = mod(my_id, num_process_x)
   t2 = my_id/num_process_x

   size_x = Nx/num_process_x
   size_z = Nz/num_process_z

   ! mpi parallel
   allocate (index_x(0:size_x - 1), index_z(0:size_z - 1))

   do i = 0, size_x - 1
      index_x(i) = size_x*t1 + i
   end do
   do i = 0, size_z - 1
      index_z(i) = size_z*t2 + i
   end do
   write (*, '(A,5I4)') 'process_id, index_x(start, end), index_z(start, end)', &
      my_id, index_x(0), index_x(size_x - 1), index_z(0), index_z(size_z - 1)

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine viscosity(kx_, kz_)
   use global
   implicit none
   real*8, intent(in)::kx_, kz_
   integer::i
   real*8::y_plus, y, lambda_G, lambda_m

   if (.not. if_resolv_eddy) then
      nu_e = 0.d0
      nu_t = 1.d0
      d_nu_t_dy = 0.d0
      return
   end if

   if (trim(eddy_model_name) == 'Cess') then
      ! get total viscosity derivative
      ! do nothing
      d_nu_t_dy = matmul(cheby_derivative_1, nu_t)
   end if

   if (trim(eddy_model_name) == 'Gupta') then
      ! get total viscosity
      lambda_G = 2*PI/sqrt(kx_*kx_ + kz_*kz_)
      do i = 0, N
         y = cheby_points(i)
         y_plus = Re_tau*(1.d0 - abs(y))
         nu_e(i) = 0.5d0*sqrt(1.d0 + (Karman*Re_tau/3.d0*(1 + y**2 - 2.d0*y**4)*(1.d0 - exp(-y_plus/Cess_A)))**2) - 0.5d0
         lambda_m = 50.d0/Re_tau + (2.d0 - 50.d0/Re_tau)*tanh(6.d0 - 6.d0*abs(y))
         nu_e(i) = nu_e(i)*lambda_G/(lambda_G + lambda_m)
      end do
      nu_t = nu_e + 1.d0
      ! get total viscosity derivative
      d_nu_t_dy = matmul(cheby_derivative_1, nu_t)
   end if

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine assemble_F(kx_, kz_)
   use global
   implicit none
   real*8, intent(in)::kx_, kz_
   integer::i
   real*8::Force_weight

   Force_spectra = 0.d0
   do i = 1, N - 1
      Force_spectra(3*i + 0, 3*i + 0) = 1.d0/integral_coeff(i)
      Force_spectra(3*i + 1, 3*i + 1) = 1.d0/integral_coeff(i)
      Force_spectra(3*i + 2, 3*i + 2) = 1.d0/integral_coeff(i)
   end do

   if (trim(force_spectra_type) == 'plain') then
      return
   end if

   if (trim(force_spectra_type) == 'eddy') then
      do i = 1, N - 1
         Force_weight = abs(nu_e(i))
         Force_spectra(3*i + 0, 3*i + 0) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 1, 3*i + 1) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 2, 3*i + 2) = 1.d0/integral_coeff(i)*Force_weight
      end do
      return
   end if

   if (trim(force_spectra_type) == 'eddy_sqr') then
      do i = 1, N - 1
         Force_weight = nu_e(i)*nu_e(i)
         Force_spectra(3*i + 0, 3*i + 0) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 1, 3*i + 1) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 2, 3*i + 2) = 1.d0/integral_coeff(i)*Force_weight
      end do
      return
   end if

   if (trim(force_spectra_type) == 'eddy_sqr_exp') then
      do i = 1, N - 1
         Force_weight = nu_e(i)*nu_e(i)*(1.d0 - exp(-abs(cheby_points(i))))
         Force_spectra(3*i + 0, 3*i + 0) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 1, 3*i + 1) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 2, 3*i + 2) = 1.d0/integral_coeff(i)*Force_weight
      end do
      return
   end if

   if (trim(force_spectra_type) == 'sweep') then
      do i = 1, N - 1
         Force_weight = m_vv(i)*(kx_*kx_*m_uu(i) + kz_*kz_*m_ww(i))**1.5
         Force_spectra(3*i + 0, 3*i + 0) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 1, 3*i + 1) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 2, 3*i + 2) = 1.d0/integral_coeff(i)*Force_weight
      end do
      return
   end if

   if (trim(force_spectra_type) == 'tke') THEN
      do i = 1, N - 1
         Force_weight = 0.5d0*(m_uu(i) + m_vv(i) + m_ww(i))
         Force_spectra(3*i + 0, 3*i + 0) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 1, 3*i + 1) = 1.d0/integral_coeff(i)*Force_weight
         Force_spectra(3*i + 2, 3*i + 2) = 1.d0/integral_coeff(i)*Force_weight
      end do
      return
   end if

   if (trim(force_spectra_type) == 'tke_component') then
      do i = 1, N - 1
         Force_spectra(3*i + 0, 3*i + 0) = 1.d0/integral_coeff(i)*m_uu(i)
         Force_spectra(3*i + 1, 3*i + 1) = 1.d0/integral_coeff(i)*m_vv(i)
         Force_spectra(3*i + 2, 3*i + 2) = 1.d0/integral_coeff(i)*m_ww(i)
      end do
      return
   end if

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine assemble_H_spatial(kx_, kz_)
   use global
   implicit none
   real*8, intent(in)::kx_, kz_
   real*8::k2
   complex*16::H_block(0:3, 0:3)
   integer::i, j

   H_spatial = 0.d0
   ! H - block-diagonal part - mean shear, mean convection
   do i = 0, N
      H_block = 0.d0
      H_block(0, 0) = IC*kx_*u(i)  ! H(1,1)
      H_block(1, 1) = IC*kx_*u(i)  ! H(2,2)
      H_block(2, 2) = IC*kx_*u(i)  ! H(3,3)
      H_block(0, 3) = IC*kx_       ! H(1,4)
      H_block(3, 0) = IC*kx_       ! H(4,1)
      H_block(2, 3) = IC*kz_       ! H(3,4)
      H_block(3, 2) = IC*kz_       ! H(4,3)
      H_block(0, 1) = du_dy(i)     ! H(1,2)
      H_spatial(4*i:4*i + 3, 4*i:4*i + 3) = H_block
   end do

   ! H - block-diagonal part - total viscosity

   k2 = kx_*kx_ + kz_*kz_
   do i = 0, N
      H_block = 0.d0
      H_block(0, 0) = k2*nu_t(i)/Re_tau + sqrt(kx_*kx_*m_uu(i) + kz_*kz_*m_ww(i))   ! H(1,1)
      H_block(1, 1) = k2*nu_t(i)/Re_tau + sqrt(kx_*kx_*m_uu(i) + kz_*kz_*m_ww(i))   ! H(2,2)
      H_block(2, 2) = k2*nu_t(i)/Re_tau + sqrt(kx_*kx_*m_uu(i) + kz_*kz_*m_ww(i))   ! H(3,3)
      H_block(0, 1) = -IC*kx_*d_nu_t_dy(i)/Re_tau  ! H(1,2)
      H_block(2, 1) = -IC*kz_*d_nu_t_dy(i)/Re_tau  ! H(3,2)
      H_spatial(4*i:4*i + 3, 4*i:4*i + 3) = H_spatial(4*i:4*i + 3, 4*i:4*i + 3) + H_block
   end do

   ! Non-diagonal part
   ! H - continuity and pressure gradient part
   do i = 0, N
      do j = 0, N
         H_spatial(4*i + 3, 4*j + 1) = H_spatial(4*i + 3, 4*j + 1) + cheby_derivative_1(i, j)
         H_spatial(4*i + 1, 4*j + 3) = H_spatial(4*i + 1, 4*j + 3) + cheby_derivative_1(i, j)
      end do
   end do

   ! H - total viscosity
   do i = 0, N
      do j = 0, N
         H_spatial(4*i + 0, 4*j + 0) = H_spatial(4*i + 0, 4*j + 0) &
                                       - 1.d0*d_nu_t_dy(i)*cheby_derivative_1(i, j)/Re_tau &
                                       - nu_t(i)*cheby_derivative_2(i, j)/Re_tau &
                                       - lambda(i)*(d_m_v_dy(j)*cheby_derivative_1(i, j) + sqrt(m_vv(i))*cheby_derivative_2(i, j))
         H_spatial(4*i + 1, 4*j + 1) = H_spatial(4*i + 1, 4*j + 1) &
                                       - 2.d0*d_nu_t_dy(i)*cheby_derivative_1(i, j)/Re_tau &
                                       - nu_t(i)*cheby_derivative_2(i, j)/Re_tau &
                                       - lambda(i)*(d_m_v_dy(j)*cheby_derivative_1(i, j) + sqrt(m_vv(i))*cheby_derivative_2(i, j))
         H_spatial(4*i + 2, 4*j + 2) = H_spatial(4*i + 2, 4*j + 2) &
                                       - 1.d0*d_nu_t_dy(i)*cheby_derivative_1(i, j)/Re_tau &
                                       - nu_t(i)*cheby_derivative_2(i, j)/Re_tau &
                                       - lambda(i)*(d_m_v_dy(j)*cheby_derivative_1(i, j) + sqrt(m_vv(i))*cheby_derivative_2(i, j))
      end do
   end do

   ! boundary condition
   ! velocity bc
   H_bc_1 = 0.d0
   H_bc_1(0:2, :) = 0.d0
   H_bc_1(3, :) = H_spatial(3, :)
   H_bc_1(0, 0) = 1.d0
   H_bc_1(1, 1) = 1.d0
   H_bc_1(2, 2) = 1.d0

   H_bc_2 = 0.d0
   H_bc_2(0:2, :) = 0.d0
   H_bc_2(3, :) = H_spatial(num_4var - 1, :)
   H_bc_2(0, num_4var - 4) = 1.d0
   H_bc_2(1, num_4var - 3) = 1.d0
   H_bc_2(2, num_4var - 2) = 1.d0

   ! pressure bc
   do i = 0, N
      H_bc_1(3, 4*i + 3) = -cheby_derivative_1(0, i)
      H_bc_1(3, 4*i + 1) = cheby_derivative_2(0, i)/Re_tau

      H_bc_2(3, 4*i + 3) = -cheby_derivative_1(N, i)
      H_bc_2(3, 4*i + 1) = cheby_derivative_2(N, i)/Re_tau
   end do

   ! apply boundary condition on resolvent matrix
   H_spatial(0:3, :) = H_bc_1
   H_spatial(num_4var - 4:num_4var - 1, :) = H_bc_2

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine assemble_H_temporal(omega_)
   use global
   implicit none
   real*8, intent(in):: omega_
   integer::i

   H = H_spatial
   ! H - time rate
   do i = 1, N - 1
      H(4*i + 0, 4*i + 0) = H(4*i + 0, 4*i + 0) - IC*omega_
      H(4*i + 1, 4*i + 1) = H(4*i + 1, 4*i + 1) - IC*omega_
      H(4*i + 2, 4*i + 2) = H(4*i + 2, 4*i + 2) - IC*omega_
   end do

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine assemble_operator_nonlinear(kx_, kz_)
   use global
   implicit none
   real*8, intent(in)::kx_, kz_
   real*8::k2
   complex*16::O_block(0:3, 0:3)
   integer::i, j

   if (if_resolv_eddy .eqv. .false.) then
      Oper_nonl = 0.d0
      return
   end if

   Oper_nonl = 0.d0
   ! block-diagonal part - eddy viscosity
   k2 = kx_*kx_ + kz_*kz_
   do i = 0, N
      O_block = 0.d0
      O_block(0, 0) = k2*nu_e(i)/Re_tau + sqrt(kx_*kx_*m_uu(i) + kz_*kz_*m_ww(i))   ! H(1,1)
      O_block(1, 1) = k2*nu_e(i)/Re_tau + sqrt(kx_*kx_*m_uu(i) + kz_*kz_*m_ww(i))   ! H(2,2)
      O_block(2, 2) = k2*nu_e(i)/Re_tau + sqrt(kx_*kx_*m_uu(i) + kz_*kz_*m_ww(i))   ! H(3,3)
      O_block(0, 1) = -IC*kx_*d_nu_t_dy(i)/Re_tau  ! H(1,2)
      O_block(2, 1) = -IC*kz_*d_nu_t_dy(i)/Re_tau  ! H(3,2)
      Oper_nonl(4*i:4*i + 3, 4*i:4*i + 3) = O_block
   end do

   ! Non-diagonal part
   ! eddy viscosity
   do i = 0, N
      do j = 0, N
         Oper_nonl(4*i + 0, 4*j + 0) = Oper_nonl(4*i + 0, 4*j + 0) &
                                       - 1.d0*d_nu_t_dy(i)*cheby_derivative_1(i, j)/Re_tau &
                                       - nu_e(i)*cheby_derivative_2(i, j)/Re_tau &
                                       - lambda(i)*(d_m_v_dy(j)*cheby_derivative_1(i, j) + sqrt(m_vv(i))*cheby_derivative_2(i, j))
         Oper_nonl(4*i + 1, 4*j + 1) = Oper_nonl(4*i + 1, 4*j + 1) &
                                       - 2.d0*d_nu_t_dy(i)*cheby_derivative_1(i, j)/Re_tau &
                                       - nu_e(i)*cheby_derivative_2(i, j)/Re_tau &
                                       - lambda(i)*(d_m_v_dy(j)*cheby_derivative_1(i, j) + sqrt(m_vv(i))*cheby_derivative_2(i, j))
         Oper_nonl(4*i + 2, 4*j + 2) = Oper_nonl(4*i + 2, 4*j + 2) &
                                       - 1.d0*d_nu_t_dy(i)*cheby_derivative_1(i, j)/Re_tau &
                                       - nu_e(i)*cheby_derivative_2(i, j)/Re_tau &
                                       - lambda(i)*(d_m_v_dy(j)*cheby_derivative_1(i, j) + sqrt(m_vv(i))*cheby_derivative_2(i, j))
      end do
   end do

   Oper_nonl(0:3, :) = 0.d0
   Oper_nonl(num_4var - 4:num_4var - 1, :) = 0.d0

end subroutine
!=========================================================================================!

!=========================================================================================!
subroutine assemble_operator_linear(kx_, kz_)
   use global
   implicit none
   real*8, intent(in)::kx_, kz_
   real*8::k2
   integer::i, j

   Oper_prod = 0.d0
   Oper_pres = 0.d0
   Oper_diss = 0.d0
   k2 = kx_*kx_ + kz_*kz_
   do i = 0, N
      ! advection
      Oper_prod(4*i + 0, 4*i + 0) = -IC*kx_*u(i)
      Oper_prod(4*i + 1, 4*i + 1) = -IC*kx_*u(i)
      Oper_prod(4*i + 2, 4*i + 2) = -IC*kx_*u(i)
      Oper_prod(4*i + 0, 4*i + 1) = -du_dy(i)

      ! pressure
      Oper_pres(4*i + 0, 4*i + 3) = -IC*kx_
      Oper_pres(4*i + 2, 4*i + 3) = -IC*kz_

      ! molecular viscosity
      Oper_diss(4*i + 0, 4*i + 0) = -k2/Re_tau
      Oper_diss(4*i + 1, 4*i + 1) = -k2/Re_tau
      Oper_diss(4*i + 2, 4*i + 2) = -k2/Re_tau
   end do

   ! eddy viscosity
   do i = 0, N
      do j = 0, N
         Oper_pres(4*i + 1, 4*j + 3) = Oper_pres(4*i + 1, 4*j + 3) - cheby_derivative_1(i, j)
         Oper_diss(4*i + 0, 4*j + 0) = Oper_diss(4*i + 0, 4*j + 0) + cheby_derivative_2(i, j)/Re_tau
         Oper_diss(4*i + 1, 4*j + 1) = Oper_diss(4*i + 1, 4*j + 1) + cheby_derivative_2(i, j)/Re_tau
         Oper_diss(4*i + 2, 4*j + 2) = Oper_diss(4*i + 2, 4*j + 2) + cheby_derivative_2(i, j)/Re_tau
      end do
   end do

   Oper_prod(0:3, :) = 0.d0
   Oper_pres(0:3, :) = 0.d0
   Oper_diss(0:3, :) = 0.d0
   Oper_prod(num_4var - 4:num_4var - 1, :) = 0.d0
   Oper_pres(num_4var - 4:num_4var - 1, :) = 0.d0
   Oper_diss(num_4var - 4:num_4var - 1, :) = 0.d0

end subroutine
!=========================================================================================!

!=========================================================================================!
SUBROUTINE INV_MAT(Nn, A, INV_A)
   INTEGER :: Nn, LDA, IPIV(Nn), INFO, LWORK
   COMPLEX*16 :: A(Nn, Nn), INV_A(Nn, Nn), WORK(64*Nn)

   INV_A = A
   LWORK = 64*Nn
   LDA = Nn

   CALL zgetrf(Nn, Nn, INV_A, LDA, IPIV, INFO)
   IF (INFO .NE. 0) THEN
      WRITE (0, *) 'Error occured in zgetrf!'
      STOP
   END IF
   CALL zgetri(Nn, INV_A, LDA, IPIV, WORK, LWORK, INFO)
   IF (INFO .NE. 0) THEN
      WRITE (0, *) 'Error occured in zgetri!'
      STOP
   END IF
END SUBROUTINE INV_MAT
!=========================================================================================!

!=========================================================================================!
subroutine write_results(kx_, kz_)
   use global
   implicit none
   real*8::kx_, kz_, y_plus
   character(len=125)::filename
   integer::i, j

   write (*, '(A,I4,A)') 'process', my_id, ' write down results'
   ! auto correlation uu
   write (filename, '(A,F6.2,A,F6.2,A)') 'results-auto/kx-', kx_, '-kz-', kz_, '-uu.dat'
   open (30, file=filename, status='replace', action='write')
   do j = 0, Nt - 1
      write (30, '(4096E16.6)') (auto_corr_uu(i, j), i=0, N)
   end do
   close (30)

   ! auto correlation vv
   write (filename, '(A,F6.2,A,F6.2,A)') 'results-auto/kx-', kx_, '-kz-', kz_, '-vv.dat'
   open (31, file=filename, status='replace', action='write')
   do j = 0, Nt - 1
      write (31, '(4096E16.6)') (auto_corr_vv(i, j), i=0, N)
   end do
   close (31)

   ! auto correlation ww
   write (filename, '(A,F6.2,A,F6.2,A)') 'results-auto/kx-', kx_, '-kz-', kz_, '-ww.dat'
   open (32, file=filename, status='replace', action='write')
   do j = 0, Nt - 1
      write (32, '(4096E16.6)') (auto_corr_ww(i, j), i=0, N)
   end do
   close (32)

   ! frequency-averaged cross correlation uu
   write (filename, '(A,F6.2,A,F6.2,A)') 'results-cross/kx-', kx_, '-kz-', kz_, '-uu.dat'
   open (40, file=filename, status='replace', action='write')
   do i = 0, N
      write (40, '(4096E16.6)') (abs(cros_corr_uu(i, j)), j=0, N)
   end do
   close (40)

   ! frequency-averaged cross correlation vv
   write (filename, '(A,F6.2,A,F6.2,A)') 'results-cross/kx-', kx_, '-kz-', kz_, '-vv.dat'
   open (41, file=filename, status='replace', action='write')
   do i = 0, N
      write (41, '(4096E16.6)') (abs(cros_corr_vv(i, j)), j=0, N)
   end do
   close (41)

   ! frequency-averaged cross correlation ww
   write (filename, '(A,F6.2,A,F6.2,A)') 'results-cross/kx-', kx_, '-kz-', kz_, '-ww.dat'
   open (42, file=filename, status='replace', action='write')
   do i = 0, N
      write (42, '(4096E16.6)') (abs(cros_corr_ww(i, j)), j=0, N)
   end do
   close (42)

   ! frequency-integrated tke production
   write (filename, '(A,F6.2,A,F6.2,A)') 'results-tke/kx-', kx_, '-kz-', kz_, '-budget.dat'
   open (50, file=filename, status='replace', action='write')
   write (50, *) 'y     y_plus      tke-nonl    tke-prod    tke-pres    tke-diss'
   do i = 0, N
      y_plus = Re_tau*(1.0 - abs(cheby_points(i)))
      write (50, '(6E16.6)') cheby_points(i), y_plus, tke_nonl(i), tke_prod(i), tke_pres(i), tke_diss(i)
   end do
   close (50)

end subroutine
!=========================================================================================!

!=========================================================================================!
program main
   use global
   implicit none

   ! dum variable
   integer::id_x, id_z, id_t, i_dum, j_dum

   call parallel_init
   call load_cfg
   call allocate_var
   call init_spectra
   call Chebyshev
   call mean_velocity
   call load_fluc_velocity
   call assemble_B
   call parallel_decompose
   call flush
   call mpi_barrier(mpi_comm_world, ierr)

   ! loop over each wave vector and frequency
   do id_x = index_x(0), index_x(size_x - 1)
      do id_z = index_z(0), index_z(size_z - 1)
         ! loop over all frequences for non-zero wave numbers
         if ((wave_x(id_x) == 0.d0) .and. (wave_z(id_z) == 0.d0)) cycle

         lambda = 0.d0
         if (if_resolv_sweep) then
            lambda = lambda_kx_y_kz(id_x, :, id_z)
         end if
         call viscosity(wave_x(id_x), wave_z(id_z))
         call assemble_F(wave_x(id_x), wave_z(id_z))
         call assemble_H_spatial(wave_x(id_x), wave_z(id_z))
         ! initial
         Resolvent_P = 0.d0
         auto_corr_uu = 0.d0
         auto_corr_vv = 0.d0
         auto_corr_ww = 0.d0
         auto_corr_pp = 0.d0
         cros_corr_uu = 0.d0
         cros_corr_vv = 0.d0
         cros_corr_ww = 0.d0
         cros_corr_pp = 0.d0
         tke_nonl = 0.d0
         tke_prod = 0.d0
         tke_pres = 0.d0
         tke_diss = 0.d0

         do id_t = 0, Nt - 1
            if (my_id == 1) then
               write (*, *) id_x, id_z, id_t
            end if
            call assemble_H_temporal(freque(id_t))

            !===================================================================================!
            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'eddy viscosity, total viscosity'
            ! do i_dum = 0, N
            !    write (*, '(I4,2E16.6)') i_dum, nu_e(i_dum), nu_t(i_dum)
            ! end do

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'H real part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (real(H(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do
            ! write (*, *) 'H imag part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (imag(H(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Force spectra'
            ! write (*, '(4096E16.6)') (real(Force_spectra(i_dum*3, i_dum*3)), i_dum=0, N)
            !===================================================================================!

            call INV_MAT(num_4var, H, Resolvent_H)

            Resolvent_R = matmul(Resolvent_H, Resolvent_B)
            Resolvent_R_trans = dconjg(transpose(Resolvent_R))

            Resolvent_F = matmul(Force_spectra, Resolvent_R_trans)
            Resolvent_P = matmul(Resolvent_R, Resolvent_F)

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Resolvent_H real part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (real(Resolvent_H(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do
            ! write (*, *) 'Resolvent_H imag part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (imag(Resolvent_H(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Resolvent P real part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (real(Resolvent_P(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do
            ! write (*, *) 'Resolvent P imag part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (imag(Resolvent_P(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do

            ! composite force spectra
            if (if_force_composite) then
               do i_dum = 1, N - 1
                  do j_dum = 1, N - 1
                     Force_spectra_comp(3*i_dum + 0, 3*j_dum + 0) = Resolvent_P(4*i_dum + 0, 4*j_dum + 0)
                     Force_spectra_comp(3*i_dum + 1, 3*j_dum + 1) = Resolvent_P(4*i_dum + 1, 4*j_dum + 1)
                     Force_spectra_comp(3*i_dum + 2, 3*j_dum + 2) = Resolvent_P(4*i_dum + 2, 4*j_dum + 2)
                  end do
               end do
               Resolvent_F = matmul(Force_spectra_comp, Resolvent_R_trans)
               Resolvent_P = matmul(Resolvent_R, Resolvent_F)
            end if

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Force spectra composite'
            ! write (*, '(4096E16.6)') (real(Force_spectra_comp(i_dum*3, i_dum*3)), i_dum=0, N)

            Resolvent_BF = matmul(Resolvent_B, Resolvent_F)
            call assemble_operator_nonlinear(wave_x(id_x), wave_z(id_z))
            call assemble_operator_linear(wave_x(id_x), wave_z(id_z))

            Res_nonl = -matmul(Oper_nonl, Resolvent_P)
            Res_prod = matmul(Oper_prod, Resolvent_P)
            Res_pres = matmul(Oper_pres, Resolvent_P)
            Res_diss = matmul(Oper_diss, Resolvent_P)

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Res_nonl real part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (real(Res_nonl(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do
            ! write (*, *) 'Res_nonl imag part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (imag(Res_nonl(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Res_prod real part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (real(Res_prod(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do
            ! write (*, *) 'Res_prod imag part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (imag(Res_prod(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Res_pres real part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (real(Res_pres(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do
            ! write (*, *) 'Res_pres imag part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (imag(Res_pres(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do

            ! write (*, *) 'wave_x, wave_z, omega', wave_x(id_x), wave_z(id_z), freque(id_t)
            ! write (*, *) 'Res_diss real part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (real(Res_diss(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do
            ! write (*, *) 'Res_diss imag part'
            ! do i_dum = 0, num_4var - 1
            !    write (*, '(4096E16.6)') (imag(Res_diss(i_dum, j_dum)), j_dum=0, num_4var - 1)
            ! end do

            ! postProcess
            block
               integer::i, j
               ! get auto correlation
               do j = 0, N
                  auto_corr_uu(j, id_t) = real(Resolvent_P(4*j + 0, 4*j + 0))
                  auto_corr_vv(j, id_t) = real(Resolvent_P(4*j + 1, 4*j + 1))
                  auto_corr_ww(j, id_t) = real(Resolvent_P(4*j + 2, 4*j + 2))
                  auto_corr_pp(j, id_t) = real(Resolvent_P(4*j + 3, 4*j + 3))
               end do

               ! get cross correlation
               do i = 0, N
                  do j = 0, N
                     cros_corr_uu(i, j) = cros_corr_uu(i, j) + Resolvent_P(4*i + 0, 4*j + 0)
                     cros_corr_vv(i, j) = cros_corr_vv(i, j) + Resolvent_P(4*i + 1, 4*j + 1)
                     cros_corr_ww(i, j) = cros_corr_ww(i, j) + Resolvent_P(4*i + 2, 4*j + 2)
                     cros_corr_pp(i, j) = cros_corr_pp(i, j) + Resolvent_P(4*i + 3, 4*j + 3)
                  end do
               end do

               ! tke production
               do i = 0, N
                  tke_nonl(i) = tke_nonl(i) + real(Res_nonl(4*i, 4*i) + Res_nonl(4*i + 1, 4*i + 1) + Res_nonl(4*i + 2, 4*i + 2)) &
                                + real(Resolvent_BF(4*i, 4*i) + Resolvent_BF(4*i + 1, 4*i + 1) + Resolvent_BF(4*i + 2, 4*i + 2))
                  tke_prod(i) = tke_prod(i) + real(Res_prod(4*i, 4*i) + Res_prod(4*i + 1, 4*i + 1) + Res_prod(4*i + 2, 4*i + 2))
                  tke_pres(i) = tke_pres(i) + real(Res_pres(4*i, 4*i) + Res_pres(4*i + 1, 4*i + 1) + Res_pres(4*i + 2, 4*i + 2) &
                                                   + Res_pres(4*i + 3, 4*i + 3))
                  tke_diss(i) = tke_diss(i) + real(Res_diss(4*i, 4*i) + Res_diss(4*i + 1, 4*i + 1) + Res_diss(4*i + 2, 4*i + 2))
               end do
            end block
         end do

         call write_results(wave_x(id_x), wave_z(id_z))

      end do
   end do
   call mpi_barrier(mpi_comm_world, ierr)

   write (*, *) 'process', my_id, 'finished'

   call deallocate_var

   call mpi_finalize(ierr)

end program main
!=========================================================================================!