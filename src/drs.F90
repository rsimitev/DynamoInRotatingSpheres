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
   use drs_io_conf
   use drs_io_state
   use drs_io
   use drs_io_units
   use drs_probes
   use drs_real_space
   use drs_debug
   use drs_lock
   use drs_error_codes
   use drs_momentum_equation
   use drs_heat_equation
   use drs_induction_equation
#ifdef COMP
   use drs_composition_equation
#endif
   implicit none

   double precision:: h1, h2, alpha !< These will cache h/2 times a few factors

   integer:: error = 0
#ifdef STAR
   double precision:: Ra_0, Ra_tau
#endif
   integer::   nargs
   character(len=2):: cargs

   nargs=command_argument_count()
   if(nargs.gt.0) then
      call get_command_argument(1,cargs)
      if (cargs=='-v') then
         Write(*,*) 'DRS Dynamo in Rotating Sphere ',VERSION
#ifdef STAR
         Write(*,*) '(using time dependent Rayleigh number: Rt_tau = ', Ra_tau, ')'
#endif
#ifdef COMP
         Write(*,*) '(Compiled with support for compositional convection.)'
#endif
         stop 1
      endif
   endif
   call drs_init(error)
   call evaluate_real_space()

   if (need_to_save()) then
      call save_stuff(nsample)
   endif

   !-- the rhs-fields rhs_NS_pol, rhs_NS_tor, .., rhs_IE_pol are in (l, m, r) space.
   !  take one Newton step coupled to Crank-Nicholson (h/2 for the laplacian)
   h1 = h
   h2 = 0.0d0
   ! Compute the RHS using the old values
   if(flow_evolves)  call NavierStokes(rhs_NS_tor, rhs_NS_pol)
   if(temp_evolves)  call TemperatureEquation(rhs_TE, error)
   if(field_evolves) call InductionEquation(rhs_IE_tor, rhs_IE_pol)
#ifdef COMP
   if(comp_evolves)  call CompositionEquation(rhs_CE, error)
#endif
   ! since the lhs whas already computed all we need to do is update
   if(flow_evolves) then
      call update_flow(h, h1, h2, error)
   else ! Just apply the boundary conditions
      call just_apply_flow_BC()
   endif ! flow_evolves

   ! Temperature
   if(temp_evolves) then
      call update_temp(h, h1, h2, Pt, error)
   elseif(temp_present) then
      call just_apply_temp_BC()
   endif

   ! Field
   if(field_evolves) then
      call update_field(h, h1, h2, Pm, error)
   elseif(field_present) then
      call just_apply_mag_BC()
   endif ! field_evolves

#ifdef COMP
   ! Composition
   if(comp_evolves) then
      call update_comp(h, h1, h2, Pc, error)
   elseif(comp_present) then
      call just_apply_comp_BC()
   endif
#endif

   !debug_call save_lmr_quantity(flow_pol, 'flowPol')
   !debug_call save_lmr_quantity(flow_tor, 'flowTor')

   call cpu_time(cpu_time_first_step)
   ! all fields in (lmr) now.
   ! end of first Newton step
   !-------------------------------------------------------------
   ! loop over time steps:
   main_loop: do while (need_to_step())
      call drs_time_update()
#ifdef STAR
      Ra_t = Ra_0 / exp(time/Ra_tau)
#endif
      call evaluate_real_space()
      ! This needs to be done before any update is made otherwise we
      ! screw up the physics.
      if (need_to_save()) then
         call drs_real_space_compute_cfl()
         call save_cfl()
         ! Check that the temporal resolution is appropriate
         ! The adaptive time-stepping is unstable so lets only adapt
         ! the time step when we save quantities
         call update_timestep(h, h_old, error)
         if (error.lt.0) Write(*,*) "Warning: Timestep is bigger than suggested by the CFL conditions!"
         call save_stuff(nsample)
      endif

      if (variable_h .and. (abs(h-h_old).gt.tiny(0.0d0))) then
         ! It generates matrices for the field with factor Pm.
         if(field_evolves) then
            call cnBp_init(inv_lhs_IE_pol, h, Pm, error)
            if(error.gt.0) call drs_abort(error)
            call cnBt_init(inv_lhs_IE_tor, h, Pm, error)
            if(error.gt.0) call drs_abort(error)
         endif

         ! It generates matrices for the flow with factor 1.0.
         if(flow_evolves) then
            call cnp_init( inv_lhs_NS_pol, h, 1.0d0, error)
            if(error.gt.0) call drs_abort(error)
            call cnt_init( inv_lhs_NS_tor, h, 1.0d0, error)
            if(error.gt.0) call drs_abort(error)
         endif

         ! It generates matrices for the temperature with factor Pt.
         if(temp_evolves) call cntemp_init(inv_lhs_TE, h, Pt,error)
         if(error.gt.0) call drs_abort(error)
#ifdef COMP
         ! It generates matrices for the composition with factor Pc.
         if(comp_evolves) call cncomp_init(inv_lhs_CE, h, Pc, error)
         if(error.gt.0) call drs_abort(error)
#endif
         ! Update alpha
         ! This is a two step Adams-Bashford/Crank-Nicholson integration
         ! with variable time step
         alpha = 0.5d0*h/h_old
         h1 = alpha*(2*h_old+h)
         h2 = alpha*h
      else
         ! This is a two step Adams-Bashford/Crank-Nicholson integration
         ! with constant time step
         h=h_old
         alpha = 0.5d0
         h1 = alpha*3.0d0*h
         h2 = alpha*h
      endif
      if (drs_calc_type.eq.KinematicDynamo .and. mod(steps, sample_rate_steps).eq.0) call kd_grothrate()
      ! Check that the spatial resolution is appropriate
      if(field_present) call check_resolution_Hartman(Rm, error)
      if(error>0) call drs_abort(error)

      ! Compute the RHS using the old values
      if(flow_evolves)  call NavierStokes(rhs_NS_tor, rhs_NS_pol)
      if(temp_evolves)  call TemperatureEquation(rhs_TE, error)
      if(field_evolves) call InductionEquation(rhs_IE_tor, rhs_IE_pol)
#ifdef COMP
      if(comp_evolves)  call CompositionEquation(rhs_CE, error)
#endif

      ! Update quantities
      ! Flow
      if (flow_evolves) call update_flow(h, h1, h2, error)
      ! Temperature
      if (temp_evolves) call update_temp(h, h1, h2, Pt, error)
      ! Magnetic Field
      if (field_evolves) call update_field(h, h1, h2, Pm, error)
#ifdef COMP
      ! Composition
      if (comp_evolves) call update_comp(h, h1, h2, Pc, error)
#endif
   enddo main_loop
   !-------------------- end of time step loop ---------------------------
   call drs_time_update()
   call cpu_time(cpu_time_now)
   spew 'Average time per step =', (cpu_time_now-cpu_time_first_step)/(steps-stepstart)

   call dump_state()

   call drs_close_output()
   if(mpi_rank.eq.0) then
      ! Close the lock file
      call rm_lock(error)
   endif
   if(error.ne.0) call drs_abort(error)
   call mpi_cleanup()

contains

   !------------------------------------------------------------------
   !> Initialize quantities and modules necessary for the program to run.
   subroutine drs_init(error)
      implicit none
      integer, intent(inout):: error
      integer, parameter:: unit_lock=987

      call cpu_time(cpu_time_start)

      call drs_mpi_init(error)
      if(error.ne.0) stop

      spew 'DRS Dynamo in Rotating Sphere ',VERSION
#ifdef STAR
      spew '(using time dependent Rayleigh number: tau = ', Ra_tau, ')'
#endif
#ifdef COMP
      spew '(Compiled with support for compositional convection.)'
#endif
      ! Hyper-difusivity is hardcoded to false for now.
      drs_want_hypDiff = .FALSE.

      !-- in the first part there's nothing to do for slave threads
      if(mpi_rank.eq.0) then
         ! Create the lock file
         call drs_lock_init(unit_lock,'drs.lock')
         call add_lock(error)
         if (error.ne.0) call drs_abort(error)
         call drs_read_conf(error)
         if (error.ne.0) call drs_abort(error)
         call select_calc_type(error)
         if (error.ne.0) call drs_abort(error)
         if(drs_calc_type.lt.1  .or. drs_calc_type.gt.CALCTYPEMAX+Compositional) then
            spew 'Cannot perform the specified calculation LCALC: ', drs_calc_type
            error = error + 1
         endif
         call check_dims(error)
         if (error.ne.0) call drs_abort(error)
         call advert_computation_parameters()
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
      call drs_bcast(temp_evolves)
      call drs_bcast(temp_present)
#ifdef COMP
      call drs_bcast(comp_evolves)
      call drs_bcast(comp_present)
#endif
      call drs_bcast(field_evolves)
      call drs_bcast(field_present)
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
      if (flow_present)  call drs_flow_allocation()
      if (field_present) call drs_field_allocation()
      if (temp_present)  call drs_temp_allocation()
#ifdef COMP
      if(comp_present)  call drs_comp_allocation()
#endif
      call drs_probes_allocation()
      ! Load the input state
      call drs_load_state(error)
      if (error>0) call drs_abort(error)

      ! We now know a bunch of quantities. Let's share them with the world.
      call drs_bcast(steps)
      call drs_bcast(stepstart)
      call drs_bcast(stepmax)
      call drs_bcast(time)
      call drs_bcast(drift)

      ! The probes require knowledge of the time
      call drs_probes_init(time)
      ! Use one tenth of the time-step from the previous computation for our
      ! first step
      if(h.lt.tiny(0.0d0)) then
         h = hi/10
         call drs_bcast(h)
         spew 'Updating initial time-step to',h
      endif
      ! Initialise the derivatives
      ! and equation modules.
      if (flow_present)  call drs_flow_init(error)
      if (flow_evolves)  call drs_momentum_equation_init(error)
      if (field_present) call drs_field_init(error)
      if (field_evolves) call drs_induction_equation_init(error)
      if (temp_present)  call drs_temp_init(error)
      if (temp_evolves)  call drs_heat_equation_init(error)
#ifdef COMP
      if (comp_present) call drs_comp_init(error)
      if (comp_evolves) call drs_composition_equation_init(error)
#endif
      !-- all fields are in (l, m, r) space now.

#ifdef STAR
      Ra_tau = 1.0
      Ra_t = Ra_0*exp(-time/Ra_tau)
      Ra_0 = Ra_t !-- Set initial value of Rayleigh number
#endif
      ! Open output files
      call drs_open_output()
   end subroutine drs_init

   !------------------------------------------------------------------
   !> Adverts the parameters taht are to be used in this computation.
   subroutine advert_computation_parameters()
      implicit none
      character(len=9), parameter:: this(2)=(/'the inner','the outer'/)
      integer:: i
      Write(*,*) 'Computations will be carried with the following parameters:'
      Write(*,'(A11,8A7)') 'Dimensions:','Nr','Nt','Np','Nr_s','Nt_s','Np_s','lsymm','m0'
      Write(*,'(A11,8I7)') ' ', Nr,  Nt,  Np,  Nr_s,  Nt_s,  Np_s,  lsymm,  m0
#ifdef COMP
      Write(*,'(8A11)') 'Parameters:','eta','Pt','Pc','Tau','Ra_t','Ra_c','Pm'
      Write(*,'(A11,7D11.4)') ' ',eta,  Pt,  Pc,  Ta,  Ra_t,  Ra_c,  Pm
      if (comp_present) then
         spew 'Compositional profile chosen was:', compProf, ' - ', compProfName()
      endif
#else
      Write(*,'(6A11)') 'Parameters:','eta','Pt','Tau','Ra_t','Pm'
      Write(*,'(A11,5D11.4)') ' ',eta,  Pt,  Ta,  Ra_t,  Pm
#endif
      if (temp_present) then
         spew 'Thermal profile chosen was:      ', tempProf, ' - ', tempProfName()
      endif
      boundaries: do i=1,2 ! loop over boundaries 1=inner, 2=outer
         Write(*,*) 'Conditions for the '//this(i)//' boundary:'
         ! Boundary conditions for the temperature
         select case(tempBC(i))
            case(FixedTemperature)
               call logFeature('Heating',tempBC(i),'Fixed temperature at '//this(i)//' boundary.')
            case(FixedHeatFlux)
               call logFeature('Heating',tempBC(i),'Fixed Heat flux at '//this(i)//' boundary.')
            case default
               Write(*,*) 'Unsupported temperature boundary conditions for '//this(i)//' boundary.',tempBC(i)
               error = ERR_UNKNOWN_TEMP_BC
         end select
#ifdef COMP
         ! Boundary conditions for the composition
         if(comp_present) then
            select case(compBC(i))
               case(FixedComposition)
                  call logFeature('Composition',compBC(i),'Fixed composition at '//this(i)//' boundary.')
               case(FixedChemFlux)
                  call logFeature('Composition',compBC(i),'Compositional flux release at '//this(i)//' boundary.')
               case default
                  Write(*,*) 'Unsupported composition boundary conditions for '//this(i)//' boundary: ',compBC(i)
                  error = ERR_UNKNOWN_COMP_BC
            end select
         endif
#endif
         ! Boundary conditions for the flow
         select case(flowBC(i))
            case(FreeSlip)
               call logFeature('Velocity',flowBC(i),'Free slip at '//this(i)//' boundary.')
            case(NoSlip)
               call logFeature('Velocity',flowBC(i),'No slip at '//this(i)//' boundary')
            case default
               Write(*,*) 'Unsuported velocity boundary conditions for '//this(i)//' boundary: ',flowBC(i)
               error = ERR_UNKNOWN_FLOW_BC
         end select

         ! Boundary conditions for the magnetic field
         select case(magBC(i))
            case(Vacuum)
               call logFeature('Magnetic',magBC(i),'vacuum at '//this(i)//' boundary.')
            case(PseudoVacuum)
               call logFeature('Magnetic',magBC(i),'Pseudo-vacuum at '//this(i)//' boundary.')
            case default
               Write(*,*) 'Unsuported magnetic boundary conditions for '//this(i)//' boundary: ',magBC(i)
               error = ERR_UNKNOWN_MAG_BC
         end select
      enddo boundaries
      ! Now that we have all the configuration parameters let us decide whether
      ! we want a variable tim-step or not
      if (h .le. 0.0d0) then
         h = -h
         variable_h = .TRUE.
         Write(*,*) 'Using a variable time-step with initial value (h) =', h
      else
         variable_h = .FALSE.
         Write(*,*) 'Using a constant time-step with value (h) =', h
      endif

      if(cpu_max_time.gt.0) then
         Write(*,*) 'Will stop after ',cpu_max_time,'cpu hours.'
      else
         Write(*,*) 'No limit on cpu time.'
      endif
      if(stepmax.gt.0) then
         Write(*,*) 'Will stop after ',stepmax, 'steps.'
      else
         Write(*,*) 'No limit on iteration number.'
      endif
   end subroutine advert_computation_parameters


   !----------------------------------------------------------------------
   !> This function will determine whether we need to take another time-step or not. There are several stoping conditions:
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
         if(need_to_save()) need_to_step = need_to_step .and. lockExists()
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

   !----------------------------------------------------------------------
   !> Computes the grothrate of the kinamatic dynamo.
   subroutine kd_grothrate()
      implicit none

      double precision:: norm, evr, evi, evr_avg, evi_avg
      double precision:: ar, bre, ai, bi, evi_min, evi_max, evr_min, evr_max
      integer:: i1, l1, m1, i2, l2, m2, i3, l3, m3, i4, l4, m4
      integer:: i, j, jg, m, l
      integer, pointer:: blk_size, blk_start

      blk_size  => blk_ps_size(mpi_rank)
      blk_start => blk_ps_start(mpi_rank)

      !-- estimate the eigenvalue from the poloidal field
      evi_min = 900000.d0
      evr_min = 900000.d0
      evi_max = -900000.d0
      evr_max = -900000.d0
      evr_avg = 0.d0
      evi_avg = 0.d0
      norm    = 0.d0
      do i=2, Nr-1
         do m=m0, m0*(blk_size-1)/2, m0
            jg = 2*(m/m0)
            j  = jg - blk_start+1
            do l=m, Nt_s
               bre = field_pol(l, j, i)
               bi  = field_pol(l, j+1, i)
               ar  = field_pol_ddr(l, j, i)/Pm   + rhs_IE_pol(l, j, i)
               ai  = field_pol_ddr(l, j+1, i)/Pm + rhs_IE_pol(l, j+1, i)
               evr = ar*bre + ai*bi
               evi = ai*bre - ar*bi
               if ((bre**2+bi**2).gt.tiny(0.0d0)) then
                  evr = evr/(bre**2+bi**2)
                  evi = abs(evi)/(bre*bre+bi*bi)
               endif
               evr = evr - llp1(l)/rcoll(i)**2/Pm
               if (evi.lt.evi_min) then
                  evi_min = evi
                  i1 = i
                  l1 = l
                  m1 = m
               endif
               if (evi.gt.evi_max) then
                  evi_max = evi
                  i2 = i
                  l2 = l
                  m2 = m
               endif
               if (evr.lt.evr_min) then
                  evr_min = evr
                  i3 = i
                  l3 = l
                  m3 = m
               endif
               if (evr.gt.evr_max) then
                  evr_max = evr
                  i4 = i
                  l4 = l
                  m4 = m
               endif
               evr_avg = evr_avg + evr*sqrt(field_pol(l, j, i)**2 + field_pol(l, j+1, i)**2)
               evi_avg = evi_avg + evi*sqrt(field_pol(l, j, i)**2 + field_pol(l, j+1, i)**2)
               norm    = norm    +     sqrt(field_pol(l, j, i)**2 + field_pol(l, j+1, i)**2)
            enddo
         enddo
      enddo

      call drs_minimize(evi_min)
      call drs_minimize(evr_min)
      call drs_maximize(evi_max)
      call drs_maximize(evr_max)
      call sum_over_all_cpus(evr_avg)
      call sum_over_all_cpus(evi_avg)
      call sum_over_all_cpus(norm)

      evr_avg = evr_avg/norm
      evi_avg = evi_avg/norm
      if(mpi_rank.eq.0) then
         spew 'average pol. ev(r, i) = ', evr_avg, evi_avg
         !      write(*, *) 'minimum evi= ', evi_min, ' at ', i1, l1, m1
         !      write(*, *) 'maximum evi= ', evi_max, ' at ', i2, l2, m2
         !      write(*, *) 'minimum evr= ', evr_min, ' at ', i3, l3, m3
         !      write(*, *) 'maximum evr= ', evr_max, ' at ', i4, l4, m4
         !      write(*, *) ' '
         write (1, *) steps, evr_avg, evi_avg
      endif

      !-- estimate eigenvalues from toroidal field
      evi_min=900000.d0
      evr_min=900000.d0
      evi_max=-900000.d0
      evr_max=-900000.d0
      evr_avg=0.d0
      evi_avg=0.d0
      norm=0.d0
      do i=2, Nr-1
         do m=m0, m0*(Np_s-1)/2, m0
            jg=2*(m/m0)
            j = jg-blk_start+1
            do l=m, Nt_s
               bre=field_tor(l, j, i)
               bi=field_tor(l, j+1, i)
               ar=field_tor_ddr(l, j, i)/pm+rhs_IE_tor(l, j, i)
               ai=field_tor_ddr(l, j+1, i)/pm+rhs_IE_tor(l, j+1, i)
               evr=ar*bre+ai*bi
               evi=ai*bre-ar*bi
               if ((bre*bre+bi*bi).gt.tiny(0.0d0)) then
                  evr=evr/(bre*bre+bi*bi)
                  evi=abs(evi)/(bre*bre+bi*bi)
               endif
               evr=evr-llp1(l)/rcoll(i)**2/pm
               if (evi.lt.evi_min) then
                  evi_min=evi
                  i1=i
                  l1=l
                  m1=m
               endif
               if (evi.gt.evi_max) then
                  evi_max=evi
                  i2=i
                  l2=l
                  m2=m
               endif
               if (evr.lt.evr_min) then
                  evr_min=evr
                  i3=i
                  l3=l
                  m3=m
               endif
               if (evr.gt.evr_max) then
                  evr_max=evr
                  i4=i
                  l4=l
                  m4=m
               endif
               evr_avg = evr_avg + evr*sqrt(field_tor(l, j, i)**2+ field_tor(l, j+1, i)**2)
               evi_avg = evi_avg + evi*sqrt(field_tor(l, j, i)**2+ field_tor(l, j+1, i)**2)
               norm    = norm    +     sqrt(field_tor(l, j, i)**2+ field_tor(l, j+1, i)**2)
            enddo
         enddo
      enddo

      call drs_minimize(evi_min)
      call drs_minimize(evr_min)
      call drs_maximize(evi_max)
      call drs_maximize(evr_max)
      call sum_over_all_cpus(evr_avg)
      call sum_over_all_cpus(evi_avg)
      call sum_over_all_cpus(norm)

      evr_avg = evr_avg/norm
      evi_avg = evi_avg/norm
      if(mpi_rank.eq.0) then
      spew 'average tor. ev(r, i) = ', evr_avg, evi_avg
      !      print *, 'minimum evi= ', evi_min, ' at ', i1, l1, m1
      !      print *, 'maximum evi= ', evi_max, ' at ', i2, l2, m2
      !      print *, 'minimum evr= ', evr_min, ' at ', i3, l3, m3
      !      print *, 'maximum evr= ', evr_max, ' at ', i4, l4, m4
      !      print *, ' '
      write (2, *) steps, evr_avg, evi_avg
      endif
   end subroutine kd_grothrate

end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
