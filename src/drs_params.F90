! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Provides the adimensional parameters for the simulation.
!! Also provides the internal flags that allow to perform different types of simulations.
#include "drsDefs.F90"
module drs_params
   use drs_error_codes
   implicit none
   save
   double precision:: eta  !< The aspect ratio
   double precision:: Pt   !< The thermal Prandtl number
   double precision:: Ra_t !< The thermal Rayleigh number
#ifdef COMP
   double precision:: Pc   !< The compositional Prandtl number
   double precision:: Ra_c !< The compositional Rayleigh number
#endif
   double precision:: Ta !< The taylor number
   double precision:: Pm !< The magnetic Prandtl number
   double precision:: pi = 3.141592653589793d0

#ifdef COMP
   double precision, target:: dconsts(8)
#else
   double precision, target:: dconsts(6)
#endif
   !< Format of the input file
   integer:: lform
   !< The type of calculation.
   !!~~~~~
   !! LinearThermalOnset  1
   !! KinematicDynamo     2
   !! NonlinearConvection 3
   !! NonlinearDynamo     4
   !! MagnetoConvection   5
   !!~~~~~
   integer:: drs_calc_type
   !< Convection type
   !!~~~~~
   !! Thermal only        1
   !! Composition only    2
   !! Doubly diffusive    3
   !!~~~~~
   integer:: conv_type
   !< The type of boundary conditions to apply to the temperature field.
   !!~~~~~
   !! FixedTemperature   0
   !! FixedHeatFlux      1
   !!~~~~~
   integer:: tempBC(2)
   !< Temperature profiles:
   !!~~~~~
   !! Conduction       0
   !! InternalHeating  1
   !!~~~~~
   integer:: tempProf
   !< Composition profiles:
   !!~~~~~
   !! WellMixed       0
   !! InternalSources 1
   !! Diffusive       2
   !!~~~~~
   integer:: compProf
#ifdef COMP
   !< The type of boundary conditions to apply to the composition field.
   !!~~~~~
   !! BottomRelease   0
   !!~~~~~
   integer:: compBC(2)
#endif
   !< The type of boundary conditions to apply to the flow.
   !!~~~~~
   !! FreeSlip   0
   !! NoSlip     1
   !!~~~~~
   integer:: flowBC(2)
   !< The type of boundary conditions to apply to the magnetic field.
   !!~~~~~
   !! Vacuum          0
   !! PseudoVacuum    1
   !!~~~~~
   integer:: magBC(2)
   !> Utility array used to pack the model options befor distribution through mpi.
   integer, target:: models(5)
   !> Does the field evolve?
   logical:: field_evolves = .TRUE.
   !> Does the flow evolve?
   logical:: flow_evolves = .TRUE.
   !> Does the temperature evolve?
   logical:: temp_evolves = .TRUE.
#ifdef COMP
   !> Does the composition evolve?
   logical:: comp_evolves = .FALSE.
#endif
   !> Is the field present?
   logical:: field_present = .TRUE.
   !> Is the flow present?
   logical:: flow_present = .TRUE.
   !> Is the temperature present?
   logical:: temp_present = .TRUE.
#ifdef COMP
   !> Is the composition present?
   logical:: comp_present = .FALSE.
#endif
   !> Noise
   double precision:: noise
   !> The input and output file names
   character(len=256):: io_calc_file_in, io_calc_file_out
   !> A comment that will identify the simulation.
   character(len=256):: comment


contains
   !------------------------------------------------------------------
   !> Packs constants and model parameters for mpi.
   subroutine drs_params_init()
      implicit none
      dconsts(1) = eta
      dconsts(2) = Ra_t
      dconsts(3) = Pt
      dconsts(4) = Ta
      dconsts(5) = Pm
      pi  = 3.141592653589793d0
      dconsts(6) = pi
#ifdef COMP
      dconsts(7) = Ra_c
      dconsts(8) = Pc
#endif

      models(1) = lform
      models(2) = drs_calc_type
   end subroutine drs_params_init

   !------------------------------------------------------------------
   !> Utility subroutine to pretty-print features
   !! @par feature is the feature to be logged.
   !! @par t is the value for the feature.
   !! @par desc is a description of the feature and the effect of the value.
   subroutine logFeature(feature,t,desc)
      implicit none
      integer, intent(in):: t
      character(len=*), intent(in)::feature, desc
      character(len=20):: fname

      fname = adjustl(feature//':')
      if (fname(20:20).ne.' ') fname(17:20)='...:'
      Write(*,'(A20,I5,A,A)') fname, t, ' - ', trim(desc)
   end subroutine

   !------------------------------------------------------------------
   !> Sets the *_present and *_evolves flags based on the value of drs_calc_type
   subroutine select_calc_type(error)
      implicit none
      integer, intent(out):: error

      ! What are we going to calculate?
      !TODO Figure out a better way to do this. Maybe something like the file
      !     check algorithm in drs_load_init_state()?
      select case(drs_calc_type)
         case(LinearThermalOnset)
            call logFeature('Calculation',drs_calc_type, 'Linear thermal onset')
            field_present = .false.
            field_evolves = .false.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
#ifdef COMP
            comp_present  = .false.
            comp_evolves  = .false.
         case(LinearCompositionalOnset)
            call logFeature('Calculation',drs_calc_type, 'Linear compositional onset')
            field_present = .false.
            field_evolves = .false.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .false.
            temp_evolves  = .false.
            comp_present  = .true.
            comp_evolves  = .true.
#endif
         case(KinematicDynamo)
            call logFeature('Calculation',drs_calc_type, 'Kinematic dynamo')
            field_present = .true.
            field_evolves = .true.
            flow_present  = .true.
            flow_evolves  = .false.
            temp_present  = .false.
            temp_evolves  = .false.
#ifdef COMP
            comp_present  = .false.
            comp_evolves  = .false.
#endif
         case(NonlinearConvection)
            call logFeature('Calculation',drs_calc_type, 'Nonlinear thermal convection')
            field_present = .false.
            field_evolves = .false.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
#ifdef COMP
            comp_present  = .false.
            comp_evolves  = .false.
#endif
         case(NonlinearDynamo)
            call logFeature('Calculation',drs_calc_type, 'Nonlinear dynamo (thermal)')
            field_present = .true.
            field_evolves = .true.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
#ifdef COMP
            comp_present  = .false.
            comp_evolves  = .false.
#endif
         case(MagnetoConvection)
            call logFeature('Calculation',drs_calc_type, 'Magnetoconvection (thermal)')
            field_present = .true.
            field_evolves = .true.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
#ifdef COMP
            comp_present  = .false.
            comp_evolves  = .false.
         case(LinearThermalOnset+Compositional)
            call logFeature('Calculation',drs_calc_type, 'Linear onset (thermal+compositional)')
            field_present = .false.
            field_evolves = .false.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
            comp_present  = .true.
            comp_evolves  = .true.
         case(NonlinearConvection+Compositional)
            call logFeature('Calculation: ',drs_calc_type, 'Nonlinear convection (thermal+compositional)')
            field_present = .false.
            field_evolves = .false.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
            comp_present  = .true.
            comp_evolves  = .true.
         case(NonlinearDynamo+Compositional)
            call logFeature('Calculation: ',drs_calc_type, 'Nonlinear dynamo (thermal+compositional)')
            field_present = .true.
            field_evolves = .true.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
            comp_present  = .true.
            comp_evolves  = .true.
         case(MagnetoConvection+Compositional)
            call logFeature('Calculation: ',drs_calc_type, 'Magnetoconvection (thermal+compositional)')
            field_present = .true.
            field_evolves = .true.
            flow_present  = .true.
            flow_evolves  = .true.
            temp_present  = .true.
            temp_evolves  = .true.
            comp_present  = .true.
            comp_evolves  = .true.
#endif
         case default
#ifndef COMP
            if(drs_calc_type.gt.10) then
               Write(*,*) 'This code was compiled without composition suport'
               error = ERR_NO_COMP_SUPPORT
               return
            endif
#endif
            Write(*,*) 'Unknown computation type: ', drs_calc_type
            error = ERR_UNKNOWN_CALC_TYPE
            return
      end select
      if((.not.field_evolves).or.(.not.field_present)) Pm=1.0d0
   end subroutine select_calc_type
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
