! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Module that takes charge of reading and writing .conf and .log files.
module drs_io_conf
#include "drsDefs.F90"
   use drs_params
   use drs_dims
   use drs_time
   use parser
   use drs_error_codes
   implicit none
   private
   save

   integer, parameter:: nConfSections=5
   integer, parameter:: nConfEntries=39

   type(parser_var):: conf_entries(nConfEntries)
   type(parser_section):: conf_sections(nConfSections)

   public:: drs_read_conf, drs_write_conf

contains

   !-----------------------------------------------------------------
   !> Describes the known config key and section entries
   subroutine drs_conf_known_vars()
      implicit none
      integer:: i
      ! Input and output file names.
      conf_sections(1)%sectionName='Filenames'
      do i=1, 2
         conf_entries(i)%section=conf_sections(1)%sectionName
      enddo
      conf_entries(1)%varName='io_calc_file_in'
      conf_entries(2)%varName='io_calc_file_out'
      ! Parameters describing the geometry.
      conf_sections(2)%sectionName='Geometry'
      do i=3,11
         conf_entries(i)%section=conf_sections(2)%sectionName
      enddo
      conf_entries(3)%varName='eta' ! Aspect ratio
      conf_entries(4)%varName='lsymm' ! Equatorial symmetry
      conf_entries(5)%varName='m0' ! Azimuthal symmetry
      conf_entries(6)%varName='Nr' ! Radial resolution
      conf_entries(7)%varName='Nt' ! Meridional resolution
      conf_entries(8)%varName='Np' ! Azimuthal resolution
      conf_entries(9)%varName='Nr_s' !
      conf_entries(10)%varName='Nt_s'
      conf_entries(11)%varName='Np_s'
      ! Parameters describing the boundary conditions.
      conf_sections(3)%sectionName='Boundaries'
      do i=12, 21
         conf_entries(i)%section=conf_sections(3)%sectionName
         ! These will be checked later
         conf_entries(i)%isRequired=.false.
      enddo
      conf_entries(12)%varName='tempBC_i'
      conf_entries(13)%varName='tempBC_o'
      conf_entries(14)%varName='tempProf'
#ifdef COMP
      conf_entries(15)%varName='compBC_i'
      conf_entries(16)%varName='compBC_o'
      conf_entries(17)%varName='compProf'
#endif
      conf_entries(18)%varName='flowBC_i'
      conf_entries(19)%varName='flowBC_o'
      conf_entries(20)%varName='magBC_i'
      conf_entries(21)%varName='magBC_o'
      ! adimensional parameters to be used.
      conf_sections(4)%sectionName='Adimensional'
      do i=22, 27
         conf_entries(i)%section=conf_sections(4)%sectionName
         ! These will be checked later
         conf_entries(i)%isRequired=.false.
      enddo
      conf_entries(22)%varName='Ta' ! Taylor
      conf_entries(22)%isRequired=.true. ! Taylor number is always required
      conf_entries(23)%varName='Pt' ! Thermal Prandtl
      conf_entries(26)%varName='Ra_t' ! Rayleigh thermal
#ifdef COMP
      conf_entries(24)%varName='Pc' ! Compositional Prandtl
      conf_entries(27)%varName='Ra_c' ! Rayleigh compositional
#endif
      conf_entries(25)%varName='Pm' ! Magnetic Prandtl
      ! Parameters that controll the runtime aspects of the simulation.
      conf_sections(5)%sectionName='Runtime'
      do i=28, nConfEntries
         conf_entries(i)%section=conf_sections(5)%sectionName
      enddo
      conf_entries(28)%varName='lform' ! Input/output data format.
      conf_entries(29)%varName='drs_calc_type' ! Calculation to be performed.
      conf_entries(30)%varName='h' ! Simulation time step.
      conf_entries(31)%varName='stepmax' ! Maximum number of time steps.
      conf_entries(32)%varName='cpu_max_time' ! Maximum cpu time.
      conf_entries(33)%varName='transient' ! How many steps constitute a transient.
      conf_entries(34)%varName='sample_on_sim_time' ! True if the sampling is to be made on simulation time intervals
      conf_entries(35)%varName='sample_rate' ! Simulation time or steps between probe samples.
      conf_entries(36)%varName='state_save_rate' ! Number of steps between intermediate state saves.
      conf_entries(37)%varName='state_save_on_sim_time' ! True if the intermediate saving is to be made on simulation time intervals
      conf_entries(38)%varName='comment' ! A comment describing the calculation.
      conf_entries(39)%varName='noise' ! Some noise in the initial temperature.
   end subroutine

   !-----------------------------------------------------------------
   !> Read the config file for the main program.
   subroutine drs_read_conf(error)
      implicit none
      integer, intent(out):: error
      logical:: isValidSection, isValidEntry
      character(len=256):: line
      character(len=60):: varname, currentSection

      ! 
      call drs_conf_known_vars()
      open(unit=444, file='drs.conf', status='old', iostat=error)
      if (error.ne.0) then
         ! The errors from the open statement are not standard accross
         ! compilers. Set our own.
         error = ERR_CONFIG_UNOPENABLE
         return
      endif
      ! Read the config file line by line and set variables accordingly
      conf_line:do
         call parse(444, varname, line, error)
         if (error.gt.0) return
         if (error.lt.0) then ! We reached the end of the file
            error = 0
            exit
         endif
         call checkEntryValidity(varname, line, conf_sections, conf_entries, currentSection, isValidEntry, isValidSection, error)
         if (error.ne.0) return
         ! Ignore section headers and read in values from valid entries.
         if (isValidEntry) call setVariableValueFromConf(line, varname, error)
         if (error.ne.0) return
      enddo conf_line
      ! Double check that we have all the required parameters for the given
      ! computation. This avoids having to set a magic number in the config
      ! file.
      call doubleCheckRequiredPars(error)
      close(444)
   end subroutine

   !------------------------------------------------------------------
   !> Double check that for a given computation, the parameters that need 
   !! to be set are indeed set.
   subroutine doubleCheckRequiredPars(error)
      implicit none
      integer, intent(out):: error
      integer:: j
      call select_calc_type(error)
      if (error.ne.0) return
      if(field_present) then
         j=entryIndexForVarName(conf_entries,'magBC_i')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'magBC_o')
         conf_entries(j)%isRequired=.true.
         if (field_evolves) then
            j=entryIndexForVarName(conf_entries,'Pm')
            conf_entries(j)%isRequired=.true.
         endif
      endif
      if(flow_present) then
         j=entryIndexForVarName(conf_entries,'flowBC_i')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'flowBC_o')
         conf_entries(j)%isRequired=.true.
         if (flow_evolves) then
            j=entryIndexForVarName(conf_entries,'Ta')
            conf_entries(j)%isRequired=.true.
         endif
      endif
      if(temp_present) then
         j=entryIndexForVarName(conf_entries,'tempBC_i')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'tempBC_o')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'tempProf')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'Ra_t')
         conf_entries(j)%isRequired=.true.
         if (temp_evolves) then
            j=entryIndexForVarName(conf_entries,'Pt')
            conf_entries(j)%isRequired=.true.
         endif
      endif
#ifdef COMP
      if(comp_present) then
         j=entryIndexForVarName(conf_entries,'compBC_i')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'compBC_o')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'compProf')
         conf_entries(j)%isRequired=.true.
         j=entryIndexForVarName(conf_entries,'Ra_c')
         conf_entries(j)%isRequired=.true.
         if (comp_evolves) then
            j=entryIndexForVarName(conf_entries,'Pc')
            conf_entries(j)%isRequired=.true.
         endif
      endif
#endif
      do j=1, nConfEntries
         if (conf_entries(j)%isRequired .and. (.not.conf_entries(j)%isSet)) then
            Write(*,*) trim(conf_entries(j)%varName), ': entry is required but could not be found.', j
            Write(*,*) trim(conf_entries(j)%varName), ': please add it under ', trim(conf_entries(j)%section)
         endif
      enddo
      ! Validate the periodic probing and saving.
      if(sample_on_sim_time) then
         sample_rate_steps = -1
      else
         sample_rate_steps = nint(sample_rate_sim_time)
         sample_rate_sim_time=-1.0d0
      endif
      if(state_save_on_sim_time) then
         state_save_rate_steps=-1
      else
         state_save_rate_steps=nint(state_save_rate_sim_time)
         state_save_rate_sim_time=-1.0d0
      endif

   end subroutine

   !------------------------------------------------------------------
   !> Actually reads the configuration values into the correct variable.
   subroutine setVariableValueFromConf(line, varname, error)
      implicit none
      character(len=*), intent(in):: line
      character(len=60), intent(in):: varname !< Key name in the config file.
      integer, intent(out):: error !< 1 if the variable name was not found; 0 otherwise.
      error = 0
      select case(varname)
         case('io_calc_file_in')
            call read_val(line, io_calc_file_in)
         case('io_calc_file_out')
            call read_val(line, io_calc_file_out)
         case('comment')
            call read_val(line, comment)
         case('lform')
            call read_val(line, lform)
         case('drs_calc_type')
            call read_val(line, drs_calc_type)
         case('tempBC_i')
            call read_val(line, tempBC(1))
         case('tempBC_o')
            call read_val(line, tempBC(2))
         case('tempProf')
            call read_val(line, tempProf)
#ifdef COMP
         case('compBC_i')
            call read_val(line, compBC(1))
         case('compBC_o')
            call read_val(line, compBC(2))
         case('compProf')
            call read_val(line, compProf)
#endif
         case('flowBC_i')
            call read_val(line, flowBC(1))
         case('flowBC_o')
            call read_val(line, flowBC(2))
         case('magBC_i')
            call read_val(line, magBC(1))
         case('magBC_o')
            call read_val(line, magBC(2))
         case('eta')
            call read_val(line, eta)
         case('Pt')
            call read_val(line, Pt)
#ifdef COMP
         case('Pc')
            call read_val(line, Pt)
         case('Ra_c')
            call read_val(line, Ra_c)
#endif
         case('Ta')
            call read_val(line, Ta)
            Ta = sqrt(Ta)
         case('Ra_t')
            call read_val(line, Ra_t)
         case('Pm')
            call read_val(line, Pm)
         case('Nr')
            call read_val(line, Nr)
         case('Nt')
            call read_val(line, Nt)
         case('Np')
            call read_val(line, Np)
         case('Nr_s')
            call read_val(line, Nr_s)
         case('Nt_s')
            call read_val(line, Nt_s)
         case('Np_s')
            call read_val(line, Np_s)
         case('lsymm')
            call read_val(line, lsymm)
         case('m0')
            call read_val(line, m0)
         case('h')
            call read_val(line, h)
         case('stepmax')
            call read_val(line, stepmax)
         case('cpu_max_time')
            call read_val(line, cpu_max_time)
         case('transient')
            call read_val(line, transient)
         case('noise')
            call read_val(line, noise)
         case('sample_rate')
            call read_val(line, sample_rate_sim_time)
         case('sample_on_sim_time')
            call read_val(line, sample_on_sim_time)
         case('state_save_rate')
            call read_val(line, state_save_rate_sim_time)
         case('state_save_on_sim_time')
            call read_val(line, state_save_on_sim_time)
         case default
            error = ERR_UNKNOWN_VARIABLE
      end select
   end subroutine

   !------------------------------------------------------------------
   !> Write a config file.
   subroutine drs_write_conf(error)
      implicit none
      integer, intent(out):: error
      error = 0
   end subroutine
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
