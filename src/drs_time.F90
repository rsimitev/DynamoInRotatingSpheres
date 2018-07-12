! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Module to deal with time. It deals with both wall time and simulation time.
!! It also deals wit the time-stepping.
module drs_time
   use drs_error_codes
   implicit none
   save
   ! Simulation time
   double precision:: time_start !< The simulation time before the first step.
   double precision:: time !< The simulation time.
   double precision:: max_time !< The maximum simulation time.
   double precision:: h !< The simulation time-step.
   double precision:: h_old !< The simulation time step at the previous step.
   double precision:: drift !< The drift rate in radians per simulation time unit.
   double precision:: time_last_sample !< Time we last saved the probes values.
   double precision:: time_since_last_sample !< Time since we last saved the probes values.
   double precision:: time_last_save !< Time we last saved the probes values.
   double precision:: time_since_last_save !< Time since we last saved the probes values.
   double precision:: sample_rate_sim_time !< The sampling rate in simulation time intervals.
   double precision:: state_save_rate_sim_time !< The sampling rate in simulation time intervals.
   logical:: variable_h !< Do we have a variable time step?
   !> By default we probe at every \a sample_rate_steps steps.
   !! Setting this flag to .true. allows us to probe on \a sample_rate_sim_time
   !! time intervals.
   logical:: sample_on_sim_time=.false.
   !> By default we save intermediate states at every \a state_save_rate_steps steps.
   !! Setting this flag to .true. allows us to probe on \a state_save_rate_sim_time
   !! time intervals.
   logical:: state_save_on_sim_time=.false.
   ! Wall time
   real:: cpu_time_start!< Wall time when we started the run in seconds.
   real:: cpu_time_now!< Wall time now, in seconds.
   real:: cpu_time_first_step!< Wall time after we finished the first Newton step, in seconds.
   real:: cpu_max_time !< The maximum wall time in hours.
   ! Steps
   integer:: transient !< How many steps to consider as a transient and discard?
   integer:: sample_rate_steps !< The sampling rate in number of steps.
   integer:: state_save_rate_steps !< The intermediate state save rate in number of steps.
   integer:: stepmax !< The maximum number of steps to take
   integer:: nsample !< How many samples we took so far.
   integer:: steps !< Steps taken so far in total (includes stepstart).
   integer:: stepstart !< The step count we start the run with.
   integer:: steps_last_save !< Step index we last saved an intermediary state.
   integer:: steps_since_last_save !< Number of steps since we last saved an intermediary state.
   double precision, target:: dtimestep(3)
   integer, target:: imeasure(5)
   integer, parameter::ncfl=4
   double precision:: cfl(ncfl) !< 

contains

   !-----------------------------------------------------------------
   subroutine drs_time_init()
      implicit none
      nsample = 0
      dtimestep(1) = h
      dtimestep(2) = time
      dtimestep(3) = drift
      imeasure(1)  = transient
      imeasure(2)  = sample_rate_steps
      imeasure(3)  = stepmax
      imeasure(4)  = nsample
      imeasure(5)  = steps
      max_time = 9.0d99
   end subroutine

   !-----------------------------------------------------------------
   !> Updates the current simulation time and step index.
   subroutine drs_time_update()
      implicit none
      time  = time  + h
      steps = steps + 1
   end subroutine

   !-----------------------------------------------------------------
   !> Updates the time-step according to the cfl numbers.
   subroutine update_timestep(h, h_old, stat)
      implicit none
      double precision, intent(inout)::h
      double precision, intent(out)::h_old
      integer, intent(out):: stat
      double precision:: aux

      stat = 0
      ! Make the time-step 100 times smaller than the
      ! smallest cfl, just in case.
      !aux = 0.005D0*minval(cfl)
      aux = 0.5D0*minval(cfl)
      if (variable_h) then
         h_old = h
         if (aux.gt.1.01d0*h_old) then
            h = 1.01d0*h_old ! increase the time step more gradualy
         elseif (aux.lt.0.09d0*h_old) then
            h = 0.09d0*h_old ! decrease the time step more gradualy
         else
            h = aux
         endif
      else
         h_old = h
         if(aux.lt.h) stat = WARN_TIME_STEP_TOO_BIG
      endif
   end subroutine update_timestep

   !-----------------------------------------------------------------
   !> Updates the time we last took a sample of some quantities
   subroutine update_time_last_sample(tnew)
      implicit none
      double precision:: tnew
      time_last_sample = tnew
   end subroutine

   !-----------------------------------------------------------------
   !> Updates the time we last saved an intermediary snapshot
   subroutine update_time_last_save(tnew)
      implicit none
      double precision:: tnew
      time_last_save = tnew
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
