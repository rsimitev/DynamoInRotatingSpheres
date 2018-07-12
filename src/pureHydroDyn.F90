! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!-- drs: Dynamo in Rotating Sphere.
!-- sequential and parallel version (choose with compiler option -DMPI).
! *****************************************************************************
#include "drsDefs.F90"
program drs
   use drs_time
   use drs_params
   use drs_dims
   use drs_mpi
   use drs_radial
   use drs_io_conf
   use drs_io_state
   use drs_io
   use drs_probes
   use drs_real_space
   use drs_debug
   use drs_lock
   use drs_error_codes
   use drs_momentum_equation
   use drs_flow
   use drs_temp
   implicit none

   double precision:: h1, h2, alpha !< These will cache h/2 times a few factors

   integer:: error = 0
   integer, pointer:: blk_size, blk_start
   integer::   nargs
   character(len=2):: cargs

   nargs=command_argument_count()
   if(nargs.gt.0) then
      call get_command_argument(1,cargs)
      if (cargs=='-v') then
         Write(*,*) 'DRS Dynamo in Rotating Sphere ',VERSION
         Write(*,*) '(Pure HydroDynamics test program.)'
         stop 1
      endif
   endif
   call drs_init(error)
   call evaluate_real_space()
   call drs_real_space_compute_cfl()
   h = 5.0d-2*minval(cfl)
   h_old = h


   if (need_to_save()) then
      call save_cfl()
      call save_stuff(nsample)
   endif

   !-- the rhs-fields rhs_NS_pol, rhs_NS_tor, .., rhs_IE_pol are in (l, m, r) space.
   !  take one Newton step coupled to Crank-Nicholson (h/2 for the laplacian)
   h1 = h
   h2 = 0.0d0
   ! Compute the RHS using the old values
   call NavierStokes(rhs_NS_tor, rhs_NS_pol)
   ! since the lhs whas already computed all we need to do is update
   call update_flow(h, h1, h2, error)

   ! debug_call save_lmr_quantity(flow_pol, 'flowPol')
   ! debug_call save_lmr_quantity(flow_tor, 'flowTor')

   call cpu_time(cpu_time_first_step)
   ! all fields in (lmr) now.
   ! end of first Newton step
   !-------------------------------------------------------------
   ! loop over time steps:
   main_loop: do while (need_to_step())
      call drs_time_update()
      !spew 'Done drs_time_update()'
      call evaluate_real_space()
      !spew 'Done evaluate_real_space()'
      !-- rhs calculates the right hand sides rhs_NS_tor, rhs_TE, rhs_NS_pol, ...:
      ! This needs to be done before any update is made otherwise we
      ! screw up the physics.
      if (need_to_save()) then
         call drs_real_space_compute_cfl()
         !spew 'Done drs_real_space_compute_cfl(cfl)'
         ! Check that the temporal resolution is appropriate
         ! The adaptive time-stepping is unstable so lets only adapt
         ! the time step when we save quantities
         call update_timestep(h, h_old, error)
         call save_cfl()
         !spew 'Done update_timestep(cfl, h, h_old, error)'
         call save_stuff(nsample)
         !spew 'Done save_stuff(nsample)'
      endif

      ! This is a two step Adams-Bashford/Crank-Nicholson integration
      ! with variable time step
      alpha = 0.5d0
      h1 = alpha*3.0d0*h
      h2 = alpha*h

      ! Compute the RHS using the old values
      call NavierStokes(rhs_NS_tor, rhs_NS_pol)
      !spew 'Done NavierStokes(rhs_NS_tor, rhs_NS_pol)'

      ! Update quantities
      ! Flow
      call update_flow(h, h1, h2, error)
      !spew 'Done update_flow(h, h1, h2, error)'
   enddo main_loop
   !-------------------- end of time step loop ---------------------------
   call drs_time_update()
   call cpu_time(cpu_time_now)
   spew 'Average time per step =', (cpu_time_now-cpu_time_first_step)/(steps-stepstart)

   call dump_state()

   if(mpi_rank.eq.0) then
      close (11) !cfl
      close (13) ! u_r @Nr/2
      close (14) ! Ek
      close (15) ! flow coeffs
      close (18) ! flow dissipation
      close (28) ! Angular momentum
      ! Close the lock file
      call rm_lock(error)
      if(error.ne.0) call drs_abort(error)
   endif  !  root

   call mpi_cleanup()

contains

   !> Initialize quantities and modules necessary for the program to run.
   subroutine drs_init(error)
      implicit none
      integer, intent(inout):: error
      character(len=9), parameter:: this(2)=(/'the inner','the outer'/)
      integer, parameter:: unit_lock=987
      integer:: i

      call cpu_time(cpu_time_start)

      call drs_mpi_init(error)
      if(error.ne.0) stop

      spew 'DRS Dynamo in Rotating Sphere ',VERSION
      spew '(Pure HydroDynamics test program.)'
      ! Hyper-difusivity is hardcoded to false for now.
      drs_want_hypDiff = .FALSE.

      !-- in the first part there's nothing to do for slave threads
      if(mpi_rank.eq.0) then
         ! Create the lock file
         call drs_lock_init(unit_lock,'pureHydroDyn.lock')
         call add_lock(error)
         if (error.ne.0) call drs_abort(error)
         !call drs_read_conf(io_calc_file_in, io_calc_file_out, comment, error)
         ! Instead of reading from file we hardcode all of the parameters
         io_calc_file_in  = ''
         io_calc_file_out = 'HydroDynamics'
         noise = 0.0d0
         lform = 1
         drs_calc_type = 3
         flowBC = 0
         eta = 0.4
         Pt = 0.0d0
         Ta = 1.0d2
         Ra_t = 0.0d0
#ifdef COMP
         Ra_c = 0.0d0
#endif
         Pm = 0.0d0
         Nr = 31
         Nt = 64
         Np = 129
         Nr_s = 31
         Nt_s = 64
         Np_s = 129
         lsymm = 0
         m0 = 1
         h = -2.0d-4
         stepmax = 2000
         cpu_max_time = 2.0d0
         transient = 0
         sample_rate_steps = 1
         time=0.0d0
         steps = 0
         comment = 'HydroDynamic test.'

         ! Bypass the calc type selection
         !call select_calc_type(error)
         flow_present  = .true.
         flow_evolves  = .true.
         field_present = .false.
         field_evolves = .false.
         temp_present  = .false.
         temp_evolves  = .false.
#ifdef COMP
         comp_present  = .false.
         comp_evolves  = .false.
#endif
         call check_dims(error)
         if (error.ne.0) call drs_abort(error)
         Write(*,*) 'Computations will be carried with the following parameters:'
         write(*,'(A11,8A7)') 'Dimensions:','Nr','Nt','Np','Nr_s','Nt_s','Np_s','lsymm','m0'
         write(*,'(A11,8I7)') ' ', Nr,  Nt,  Np,  Nr_s,  Nt_s,  Np_s,  lsymm,  m0
         write(*,'(6A11)') 'Parameters:','eta','Pt','Tau','Ra_t','Pm'
         write(*,'(A11,5D11.4)') ' ',eta,  Pt,  Ta,  Ra_t,  Pm

         boundaries: do i=1,2 ! loop over boundaries 1=inner, 2=outer
            spew 'Conditions for the '//this(i)//' boundary:'
            ! Boundary conditions for the flow
            if(flowBC(i)==FreeSlip) then
               call logFeature('Velocity',flowBC(i),'Free slip at '//this(i)//' boundary.')
            else
               Write(*,*) 'Wrong boundary type'
               call drs_abort(1)
            endif
         enddo boundaries

         ! Now that we have all the configuration parameters let us decide whether
         ! we want a variable tim-step or not
         if (h .le. 0.0d0) then
            h = -h
            variable_h = .TRUE.
            spew 'Using a variable time-step with initial value (h) =', h
         else
            variable_h = .FALSE.
            spew 'Using a constant time-step with value         (h) =', h
         endif

         if(cpu_max_time.gt.0) then
            spew 'Will stop after ',cpu_max_time,'cpu hours.'
         else
            spew 'No limit on cpu time.'
         endif
         if(stepmax.gt.0) then
            spew 'Will stop after ',stepmax, 'steps.'
         else
            spew 'No limit on iteration number.'
         endif
      endif
      if (error.ne.0) call drs_abort(error)

      ! Start the initializations
      call drs_time_init()
      call drs_params_init()
      call drs_dims_init(error)
      if (error.ne.0) call drs_abort(error)

      ! broadcast input parameters:
      call drs_bcast(dconsts, size(dconsts))
      call drs_bcast(variable_h)
      call drs_bcast(dtimestep, size(dtimestep))
      call drs_bcast(usr_dims, size(usr_dims))
      call drs_bcast(imeasure, size(imeasure))
      call drs_bcast(models, size(models))
      call drs_bcast(flow_evolves)
      call drs_bcast(flow_present)

      eta   = dconsts(1)
      Ra_t  = dconsts(2)
      Pt    = dconsts(3)
      Ta    = dconsts(4)
      Pm    = dconsts(5)
      pi    = dconsts(6)
#ifdef COMP
      Ra_c  = dconsts(7)
      Pc    = dconsts(8)
#endif

      h     = dtimestep(1)
      time  = dtimestep(2)
      drift = dtimestep(3)
      h_old = h

      Nr    = usr_dims(1)
      Nt    = usr_dims(2)
      Np    = usr_dims(3)
      Nr_s  = usr_dims(4)
      Nt_s  = usr_dims(5)
      Np_s  = usr_dims(6)
      lsymm = usr_dims(7)
      m0    = usr_dims(8)

      transient     = imeasure(1)
      sample_rate_steps = imeasure(2)
      stepmax       = imeasure(3)
      nsample       = imeasure(4)
      steps         = imeasure(5)

      lform         = models(1)
      drs_calc_type = models(2)

      ! Now that we know how big the problem is, we can split it across the
      ! CPU's.
      call mpi_dims_init(Nt, Np_s, m0, error)
      blk_size  => blk_ps_size(mpi_rank)
      blk_start => blk_ps_start(mpi_rank)
      ! From now on we can use the per CPU sizes

      ! Initialize the dft routines.
      call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)
      call drs_legendre_allocation()
      call drs_legendre_init()
      call drs_radial_init(eta)

!      call CrankNicholson_init()
      if (drs_want_hypDiff) then
        spew '(Using Hyper-Diffusivity)'
      endif
      call drs_hypDiff_init(Nt)

      ! Allocate space for the input state.
      call drs_flow_allocation()
      call drs_temp_allocation()
      call drs_probes_allocation()
      ! Load the input state
      call load_synt_state(error)
      if (error>0) call drs_abort(error)

      ! We now know a bunch of quantities. Let's share them with the world.
      call drs_bcast(steps)
      call drs_bcast(stepstart)
      call drs_bcast(stepmax)
      call drs_bcast(time)
      call drs_bcast(drift)

      ! The probes require knowledge of the time
      call drs_probes_init(time)
      ! Initialise the derivatives
      ! and equation modules.
      call drs_flow_init(error)
      call drs_momentum_equation_init(error)
      !-- all fields are in (l, m, r) space now.

      ! Open output files
      call drs_open_output()
   end subroutine drs_init

   !----------------------------------------------------------------------
   !> This function will determine whether we need to take another time-step 
   !! or not. There are several stoping conditions:
   logical function need_to_step()
      implicit none
      real:: cpu_time_now
      if (mpi_rank==0) then
         need_to_step = .TRUE.
         !> - a stopping condition on number of steps;
         if(stepmax.gt.0) need_to_step = need_to_step .and. (steps <= stepstart + stepmax-1)
         !> - a stopping condition on cpu time;
         if (cpu_max_time .gt. 0.0e0) then
            call cpu_time(cpu_time_now)
            need_to_step = need_to_step .and. ((cpu_time_now-cpu_time_start)/3600 < cpu_max_time)
         endif
         !> - a stoping condition on the existence of a lock file;
         need_to_step = need_to_step .and. lockExists()
         !> - a stoping condition on simulation time;
         need_to_step = need_to_step .and. (time.lt.max_time)
      endif
      call drs_bcast(need_to_step)
      call wait_for_everyone()
   end function need_to_step

   !--------------------------------------------------------------------------
   !> Returns .TRUE. if it is time to save the probed quantities.
   pure logical function need_to_save()
      implicit none
      need_to_save = steps.ge.transient
      need_to_save = need_to_save .and. mod((steps-transient),sample_rate_steps).eq.0
   end function need_to_save

   !---------------------------------------------------------------
   !> Reads a state performing interpolation as needed.
   !! The state is stored in the files with name
   !! given by \a io_calc_file_in and are described by the file with extension
   !! .par.
   subroutine load_synt_state(error)
      implicit none
      integer, intent(out):: error
      double precision, allocatable:: cache(:,:,:)
      integer:: i,j,l, m
      integer:: blk_size,blk_start

      error = 0
      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)
      ! This is the model time of the input state
      time_start = 0.0d0
      if(mpi_rank.eq.0) then
         allocate(cache(Nr, 0:Nt_s, Np_s))
      else
         allocate(cache(Nr, 0:Nt_s, blk_ps_size(mpi_rank)))
      endif
      cache = 0.0d0
      if(mpi_rank.eq.0) then
         !-- Npi_s is dim of input file, blk_ps_size(mpi_rank) is dim for this calculation!
         ! If the input state contains unnormalized coefficients
         ! apply a norm transformation
         do j=1, Np_s
            m = m0i*(j/2)
            do l = m, Nt_s
               do i=1, Nr
                  if(m.eq.0 .and. l==1) then
                     cache(i,l,j) = rcoll(i)
                  else
                     cache(i,l,j) = 0.0d0
                  endif
               enddo
            enddo
         enddo
      endif
      call distribute_in_m(cache, Nt_s, Nr)

      do j=1, blk_ps_size(mpi_rank)
         do l = m0*(j/2), Nt_s
            flow_tor(l,j,1:Nr) = cache(1:Nr,l,j)
         enddo
      enddo
      deallocate(cache)
      flow_pol=0.0d0
      temp=0.0d0

   end subroutine load_synt_state
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
