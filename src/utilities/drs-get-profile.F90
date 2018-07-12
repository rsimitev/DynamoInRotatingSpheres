! -----------
!> Computes the horizontally integrated radial profile for the
!! given quantity.
! -----------

program getProfile
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
#ifdef COMP
   use drs_comp
#endif
   use drs_temp
   use drs_io_state
   use drs_probes
   use drs_error_codes
   implicit none

   integer:: error, what, i
   character(len=64):: whatName
   double precision, allocatable:: profile(:)

   call init(error)

   call system (inflate_state//trim(io_calc_file_in))
   call drs_load_state(error)
   if(error>0) call drs_abort(error)
   call system (deflate_state//trim(io_calc_file_in))

   select case(what)
      case (1)
         profile(:) = temp(0,1,:)
      case (2)
         call drs_temp_init(error)
         profile(:) = temp(0,1,:)+ temp_profile(:)
#ifdef COMP
      case (3)
         profile(:) = comp(0,1,:)
      case (4)
         call drs_comp_init(error)
         profile(:) = comp(0,1,:)+ comp_profile(:)
#endif
      case default
         profile(:) = temp(0,1,:)
   end select

   open(900, file=trim(io_calc_file_in)//'.rprof', status='NEW')
   Write(900,*) '# ',trim(whatName)
   do i=1, Nr
      write(900, *) rcoll(i), profile(i)
   enddo
   close(900)
   spew 'Saved file ', trim(io_calc_file_in)//'.rprof'
contains

   !--------------------------------------------
   !> Initialise things
   subroutine init(error)
      implicit none
      integer, intent(inout):: error
      error = 0
      call parseConfig(error)
      if (error.ne.0) then
         call usage()
         call drs_abort(error)
      endif
      call setWhatName()
      ! Read parameters in
      open(900, file=trim(io_calc_file_in)//'.par', status='OLD')
      call drs_read_state_par(900) ! This is needed before initialising the modules
      close(900)

      ! We do not use hyper diffusivity.
      drs_want_hypDiff = .FALSE.

      ! Set the resolution to use in the calculations equal to the resolution on
      ! the files
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
      allocate(profile(Nr))

   end subroutine init

   !--------------------------------------------
   !> Prints the usage on screen
   subroutine usage()
      implicit none
      Write(*,*) " In order to run this program the file getProfile.conf"
      Write(*,*) " needs to be present in the current folder."
      Write(*,*) " Please get the corresponding template from the templates folder."
   end subroutine

   !--------------------------------------------
   !> Parses the configuration file
   !! ~~~~~
   !! state = <state base name>
   !! what = quantity to generate a profile for
   !! ~~~~~
   subroutine parseConfig(error)
      use parser
      implicit none
      integer, intent(out):: error
      character(len=60):: varname
      character(len=256):: line
      integer:: nvar

      error = 0
      open(700,file='getProfile.conf',status="OLD", iostat=error)
      if(error.ne.0) then
         error=ERR_CONFIG_UNOPENABLE
         return 
      endif
      ! Read until the end of the file or until we got all the values we need.
      nvar = 0 
      do while (error.eq.0.and.nvar.lt.2)
         call parse(700, varname, line, error)
         select case(varname)
            case('state')
               call read_val(line, io_calc_file_in)
               nvar = nvar + 1
            case('what')
               call read_val(line, what)
               nvar = nvar + 1
            case default
               cycle
         end select
      enddo

      close(700)
   end subroutine parseConfig

   !--------------------------------------------
   !> Sets a whuman readable name for what is being computed.
   subroutine setWhatName()
      implicit none
      select case(what)
         case (1)
            whatName = 'Temperature anomaly'
         case (2)
            whatName = 'Temperature'
#ifdef COMP
         case (3)
            whatName = 'Composition anomaly'
         case (4)
            whatName = 'Composition'
#endif
         case default
            whatName = 'Temperature anomaly'
      end select
   end subroutine setWhatName
end program getProfile
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
