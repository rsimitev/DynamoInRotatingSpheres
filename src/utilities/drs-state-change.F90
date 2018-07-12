! -----------
!> Computes the spectra of all quantities for a speciffic state.
! -----------
program change_state
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

   character(len=2):: what_char
   integer:: error, nargs, what

   nargs=command_argument_count()
   if(nargs.ne.2) then
      Write(*,*) 'drs_change_state <state> <what>'
      Write(*,*) 'Changes the structure of a state.'
      Write(*,*) 'Takes two mandatory arguments:'
      Write(*,*) ' <state> is the base name of the state without extension;'
      Write(*,*) ' <what> is the action on the state.'
      Write(*,*) 'Will output a new state called <state>_new.'
      stop 1
   endif
   ! TODO: Add usage instructions and version number
   call get_command_argument(1,io_calc_file_in)
   call get_command_argument(2,what_char)
   read(what_char,*) what
   io_calc_file_in  = adjustl(io_calc_file_in)
   io_calc_file_out = trim(io_calc_file_in)//'_new'
   call init(error)

   call system (inflate_state//trim(io_calc_file_in))
   call drs_load_state(error)
   if(error>0) call drs_abort(error)
   call system (deflate_state//trim(io_calc_file_in))

   select case(what)
      case(0)
         call kill_differential_rotation()
   end select

   call drs_save_state()
   !if(error>0) call drs_abort(error)
   call system (deflate_state//trim(io_calc_file_out))

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
      Pm   = Pmi
      Ta   = Tai
      Ra_t = Ra_ti
#ifdef COMP
      Ra_c = Ra_ci
      Pc   = Pci
      compBC = compBCi
#endif
      tempBC = tempBCi
      flowBC = flowBCi
      magBC  = magBCi
      h = hi
      stepmax = stepmaxi
      transient = transienti
      sample_rate_steps = sample_ratei
      drift = drifti

      call drs_mpi_init(error)
      if(error.ne.0) call drs_abort(1)
      if(mpi_size.ne.1) then
         spew 'This program should be ran on a single cpu'
         call drs_abort(1)
      endif
      call select_calc_type(error)
      if(error.ne.0) call drs_abort(1)
      ! Start the initializations
      call drs_time_init()
      call drs_params_init()
      call drs_dims_init(error)
      if(error.ne.0) call drs_abort(1)

      call mpi_dims_init(Nt, Np_s, m0, error)
      if(error.ne.0) call drs_abort(1)
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

   !-----------------------------------------------------------------
   !> Sets all toroidal zonal components of the flow to zero.
   subroutine kill_differential_rotation()
      implicit none
      integer:: l, j, m, i

      do i=1, Nr
         jlm_do(j,l,m)
            if (m==0) flow_tor(l,j,i)=0.0d0
         jlm_enddo
      enddo
   end subroutine
end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
