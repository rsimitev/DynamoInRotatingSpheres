! -----------
!> Generate the initial conditions relevant to the 1st benchmark exercise based on a given state.
! -----------
program Benchmarkv1
#include "drsDefs.F90"
   use drs_params
   use drs_dims
   use drs_mpi
   use drs_fftw3
   use drs_legendre
   use drs_radial
   use drs_flow
   use drs_field
   use drs_temp
   use drs_io_state
   use drs_time
   implicit none

   double precision, allocatable:: vec_r(:,:,:)
   double precision, allocatable:: vec_t(:,:,:)
   double precision, allocatable:: vec_p(:,:,:)
   double precision, allocatable:: scal(:,:,:)
   double precision, allocatable:: x(:), r(:)
   double precision:: dp
   integer:: error, l, j, i
   double precision, parameter:: A=0.1d0
   
   error = 0
   call init(error)
   if(error.ne.0) call drs_abort(error)

   r  = rcoll
   dp = 2.0d0*pi/Np
   x  = 2*r - ri - ro

   ! Construct the temperature in real space
   ! Used for cases 0 and 1
   ! TODO: Check that this corresponds to our adimensionalization
   do i=1, Nr
      do j=1, Np
         do l=0, Nt
            scal(l,j,i) = (1.0d0 - 3*x(i)**2 + 3*x(i)**4 - x(i)**6)*sintheta(l)**4*dcos(4*(l-1)*dp)
            scal(l,j,i) = ro*ri/r(i) - ri + 210*A/(17920*pi)*scal(l,j,i)
         enddo
      enddo
   enddo
   ! Transform to spectral space
   call ylmt_3D(scal, temp)
   ! The initial flow is zero.
   flow_pol = 0.0d0
   flow_tor = 0.0d0
   ! Construct the field in real space
   ! only for benchmark case 1
   do i=1, Nr
      do j=1, Np
         do l=0, Nt
            vec_r(l,j,i) = 5.0d0/8.0d0*( 8*ro   - 6*r(i) - 2*ri**4/r(i)**3 )*costheta(l)
            vec_t(l,j,i) = 5.0d0/8.0d0*( 9*r(i) - 8*ro -     ri**4/r(i)**3 )*sintheta(l)
            vec_p(l,j,i) = 5.0d0*dsin(pi*(r(i)-ri))*2*costheta(l)*sintheta(l)
         enddo
      enddo
   enddo
   ! Transform to pol/tor
   call vectorField2PolTor_common(vec_r,vec_t,vec_p,field_pol,field_tor)
   call PolTor_common2PolTor_field(field_pol, field_tor)

   ! Save the state
   call drs_save_state()

contains

   subroutine init(error)
      implicit none
      integer, intent(inout):: error
      character(len=100):: deflate, inflate
      error = 0
      deflate = 'compressdata '
      inflate = 'uncompressdata '

      ! TODO: Rename this file to be more parameter descriptive.
      io_calc_file_out = 'benchmark-christensen-init'

      ! Initialise everything
      drs_want_hypDiff = .FALSE.

      ! Set the resolution to use in the calculations.
      m0 = 1
      Nr = 33
      Nt = 64
      Np = 129
      eta = 0.35
      Nr_s = 33
      Nt_s = 64
      Np_s = 129
      lform = 1
      drift = 0.0d0
      comment = 'Initial condition for benchmark Christensen et al. (2001)'
      drs_calc_type = 4
      ! Update the simulation parameters
      Pt   = 1.0
      Ta   = 4.0e6
      Ra_t = 100
      Pm   = 5

      call drs_mpi_init(error)
      if(mpi_size.ne.1) then
         spew 'This program should be ran on a single cpu'
         call drs_abort(1)
      endif
      ! Start the initializations
      call drs_params_init()
      call drs_dims_init(error)

      call mpi_dims_init(Nt, Np_s, m0, error)
      call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)

      call drs_legendre_allocation()
      call drs_flow_allocation()
      call drs_field_allocation()
      call drs_temp_allocation()

      call drs_legendre_init()
      call drs_radial_init(eta)

      ! Initialize derivatives

      allocate(vec_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(vec_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(vec_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(scal(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(x(Nr), r(Nr))

   end subroutine init
end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
