! -----------
!> Computes the plots require for the Yokoi paper
! -----------
program drs2dx
#include "drsDefs.F90"
   use drs_time
   use drs_params
   use drs_dims
   use drs_mpi
   use drs_radial
   use drs_flow
   use drs_field
   use drs_temp
#ifdef COMP
   use drs_comp
#endif
   use drs_io_state
   use drs_probes
   use drs_real_space
   use drs_renderers
   use drs_io_DX
   use parser
   implicit none

   integer:: error, what

   error = 0
   call init(error)

   spew "Done init"

   call render(what)
   spew "Done rendering"
   ! TODO Document what this does
   if(what/100.eq.5) then
      ! Redefine the radial coordinates
      call redefine_radial_coordinate()
      spew "Done redefining the radial coordinates"
   endif
   if((what-10*(what/10)).eq.4) then
      call save2DX(XX,YY,ZZ,io_calc_file_in)
   else
      call save2DX(render_out,io_calc_file_in)
   endif
   spew "Done saving"

contains
   subroutine init(error)
      implicit none
      integer, intent(inout):: error
      error = 0
      call parse_drs2dx()
      ! Read parameters in
      open(900, file=trim(io_calc_file_in)//'.par', status='OLD')
      call drs_read_state_par(900) ! This is needed before initialising the modules
      close(900)

      ! Initialise everything
      drs_want_hypDiff = .FALSE.

      ! The harmonic resolution is the same as in the files.
      m0 = m0i
      eta = etai
      Nr_s = Nr
      Nt_s = Nt
      Np_s = Np
      lform = lformi
      drift = drifti
      ! The simulation parameters are the ones in the file
      Pt   = Pti
      Ta   = Tai
      Ra_t = Ra_ti
      Pm   = Pmi
#ifdef COMP
      Ra_c = Ra_ci
      Pc   = Pci
#endif
      drs_calc_type = drs_calc_typei
      call check_dims(error)
      if (error.ne.0) call drs_abort(error)

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
      call drs_renderers_allocation(what)

      ! Read a new state
      spew "Initialized."
      call system (inflate_state//trim(io_calc_file_in))
      call drs_load_state(error)
      if(error>0) call drs_abort(error)
      call system (deflate_state//trim(io_calc_file_in))
      spew "Done loading state."
      call drs_temp_init(error)
#ifdef COMP
      call drs_comp_init(error)
#endif
   end subroutine init

   subroutine parse_drs2dx()
      implicit none
      integer:: error
      character(len=256):: line
      character(len=40):: varname

      error=0
      open(unit=444, file='drs2dx.in', status='OLD', iostat=error)
      spew error
      if (error.ne.0) return

      ! TODO: Update to new parser infrastructure
      do while(error.eq.0)
         call parse(444, varname, line, error)
         select case(varname)
            case('io_calc_file_in')
               call read_val(line, io_calc_file_in)
               spew 'Input: ', io_calc_file_in
            case('comment')
               call read_val(line, comment)
            case('Nr')
               call read_val(line, Nr)
            case('Nt')
               call read_val(line, Nt)
            case('Np')
               call read_val(line, Np)
            case('what')
               call read_val(line, what)
            case('render_type')
               call read_val(line, cut_type)
            case('where')
               call read_val(line, where_to_cut)
            case default
               cycle
         end select
      enddo
   end subroutine parse_drs2dx

   subroutine redefine_radial_coordinate()
      implicit none
      integer:: i
      double precision:: rout, rin
      rin  = rcoll(Nr)
      rout = rcoll(1)
      rcoll = rout + 2.0d0*(rcoll - rin)
      do i=1, Nr
         rcoll2(i) = rcoll(i)**2
         if (i.gt.1) drcoll(i-1) = rcoll(i-1) - rcoll(i)
      enddo
      drcoll(Nr) = drcoll(Nr-1)
   end subroutine redefine_radial_coordinate
end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
