! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> This module provides a locking mechanism for the dynamo code.
!! @since 1.6.1
module drs_lock
   use drs_error_codes
   implicit none
   private
   character(len=128):: lockFileName='Unknown'
   integer:: lockFileUnit=-1
   public:: drs_lock_init, add_lock, rm_lock, lockExists
contains

   !------------------------------------------------------------------
   !> Sets the lock file name to f and manages it on unit u.
   subroutine drs_lock_init(u, f)
      implicit none
      character(len=*):: f !< The name of the lock file.
      integer:: u !< The unit it is going to be openned on.
      lockFileUnit = u
      lockFileName = trim(f)
   end subroutine

   !------------------------------------------------------------------
   !> Creates the lock file
   subroutine add_lock(error)
      implicit none
      integer, intent(inout):: error
      integer:: err
      open(unit=lockFileUnit, file=trim(lockFileName), status='NEW', iostat=err)
      if (err.ne.0) then
         Write(*,*) 'Could not create lock file or lock file already exists.'
         Write(*,*) 'Probably we are already running here!'
         error = ERR_CREATING_LOCK
      else
         close(lockFileUnit)
      endif
   end subroutine

   !------------------------------------------------------------------
   !> Removes the lock file
   subroutine rm_lock(error)
      implicit none
      integer, intent(inout):: error
      integer:: err
      open(unit=lockFileUnit, file=trim(lockFileName), status='UNKNOWN', iostat=err)
      if (err.ne.0) then
         Write(*,*) 'Could not find lock file.'
         Write(*,*) 'Probably we are no longer running here!'
         error = ERR_DELETING_LOCK
      else
         close(unit=lockFileUnit, status='delete')
      endif
   end subroutine

   !------------------------------------------------------------------
   !> Checks whether the lock file exists.
   logical function lockExists()
      inquire(file=trim(lockFileName),exist=lockExists)
   end function

end module drs_lock
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
