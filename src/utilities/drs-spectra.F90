! -----------
!> Computes the spectra of all quantities for a speciffic state.
! -----------
program spectra
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
   use drs_io
   use drs_probes
   use drs_real_space
   implicit none

   integer:: error, nargs

   nargs=command_argument_count()
   if(nargs.ne.1) then
      Write(*,*) '* Invalid number of states.'
      Write(*,*) '* Pass only one state at a time.'
      stop 1
   endif
   ! TODO: Add usage instructions and version number
   call get_command_argument(1,io_calc_file_in)
   io_calc_file_in  = adjustl(io_calc_file_in)
   io_calc_file_out = io_calc_file_in
   call init(error)

   call system (inflate_state//trim(io_calc_file_in))
   call drs_load_state(error)
   if(error>0) call drs_abort(error)
   call system (deflate_state//trim(io_calc_file_in))

   call save_l_spec()
   call save_m_spec()
   call save_n_spec()

contains
   subroutine init(error)
      implicit none
      integer, intent(inout):: error
      error = 0
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
      if(mpi_size.ne.1) then
         spew 'This program should be ran on a single cpu'
         call drs_abort(1)
      endif
      call select_calc_type(error)
      ! Start the initializations
      call drs_time_init()
      call drs_params_init()
      call drs_dims_init(error)

      call mpi_dims_init(Nt, Np_s, m0, error)
      call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)

      call drs_legendre_allocation()
      call drs_flow_allocation()
      call drs_field_allocation()
      call drs_temp_allocation()
#ifdef COMP
      call drs_comp_allocation()
#endif
      call drs_probes_allocation()

      call drs_legendre_init()
      call drs_radial_init(eta)
      call drs_hypDiff_init(Nt)
      call drs_probes_init(time)

   end subroutine
end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
