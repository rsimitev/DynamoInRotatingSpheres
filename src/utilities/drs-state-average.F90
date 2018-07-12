! -----------
!> Takes an  Average in  time of all fields
! -----------
program StateAverage
#include "drsDefs.F90"
   use drs_time
   use drs_params
   use drs_dims
   use drs_mpi
   use drs_fftw3
   use drs_legendre
   use drs_radial
   use drs_hypDiff
   use drs_flow
   use drs_field
   use drs_temp
#ifdef COMP
   use drs_comp
#endif
   use drs_io_state
   use drs_probes
   use drs_real_space
   implicit none

   integer:: error, nstates
   character:: first
   io_calc_file_out = "averaged.0"

   nstates = 0
   call init(error)

   !--main loop
   do
      read (700,'(A50)', iostat=error) io_calc_file_in
      if (error.ne.0) exit
      io_calc_file_in = adjustl(io_calc_file_in)
      first = io_calc_file_in(1:1)
      if((first=='#').or.(first=='*')) cycle
      nstates = nstates + 1

      call system (inflate_state//trim(io_calc_file_in))
      call drs_load_state(error)
      if(error>0) call drs_abort(error)
      call system (deflate_state//trim(io_calc_file_in))

      flow_pol_avg = flow_pol_avg + flow_pol
      flow_tor_avg = flow_tor_avg + flow_tor
      if(field_evolves) then
         field_pol_avg = field_pol_avg + field_pol
         field_tor_avg = field_tor_avg + field_tor
      endif
      if(temp_evolves) then
         temp_avg = temp_avg + temp
      endif
   enddo
   close(700)
   flow_pol_avg = flow_pol_avg/nstates
   flow_tor_avg = flow_tor_avg/nstates
   if(field_evolves) then
      field_pol_avg = field_pol_avg/nstates
      field_tor_avg = field_tor_avg/nstates
   endif
   if(temp_evolves) then
      temp_avg = temp_avg/nstates
   endif

   call drs_save_state()
   call system (deflate_state//trim(io_calc_file_out))
contains
   subroutine init(error)
      implicit none
      integer, intent(inout):: error
      character:: first
      error = 0
      open(700,file='state-average.in',status="OLD")
      do
         read(700,'(A)', iostat=error) io_calc_file_in
         io_calc_file_in = adjustl(io_calc_file_in)
         first = io_calc_file_in(1:1)
         if((first.ne.'#').and.(first.ne.'*')) exit
      enddo
      rewind(700)
      ! Read parameters in
      open(900, file=trim(io_calc_file_in)//'.par', status='OLD')
      call drs_read_state_par(900) ! This is needed before initialising the modules
      close(900)

      ! Initialise everything
      drs_want_hypDiff = .FALSE.

      ! Set the resolution to use in the calculations.
      m0 = m0i
      Nr = Nri
      Nt = Nti
      Np = Npi
      eta = etai
      Nr_s = Nri_s
      Nt_s = Nti_s
      Np_s = Npi_s
      lform = lformi
      drift = drifti
      comment = commenti
      drs_calc_type = drs_calc_typei
      ! Update the simulation parameters
      Pt   = Pti
      Ta   = Tai
      Ra_t = Ra_ti
      Pm   = Pmi

      call drs_mpi_init(error)
      if (error.gt.0) call drs_abort(error)
      if(mpi_size.ne.1) then
         spew 'This program should be ran on a single cpu'
         call drs_abort(1)
      endif
      call select_calc_type(error)
      if (error.gt.0) call drs_abort(error)
      ! Start the initializations
      call drs_time_init()
      call drs_params_init()
      call drs_dims_init(error)
      if (error.gt.0) call drs_abort(error)

      call mpi_dims_init(Nt, Np_s, m0, error)
      call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)

      call drs_legendre_allocation()
      call drs_flow_allocation()
      call drs_field_allocation()
      call drs_temp_allocation()
      call drs_probes_allocation()

      call drs_legendre_init()
      call drs_radial_init(eta)
      call drs_hypDiff_init(Nt)
      call drs_probes_init(time)

   end subroutine
end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
