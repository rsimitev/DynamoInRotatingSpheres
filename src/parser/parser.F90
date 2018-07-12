!     Copyright (C) 2005-2015 Luis Silva <lacsilva@gmail.com>
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!>   This module provides a simple parser for input files of the type
!!   'key = val'
!!   eventually separated by '[Sections]'
!! 
!!   Typical use of this module goes like:
!! ~~~~~~~ 
!! open(unit=444, file='myconfig.conf', status='old', iostat=error)
!! if (error.ne.0) return
!! do
!!    call parse(444, varname, line, error)
!!    select case(varname)
!!       case('var1')
!!          call read_val(line, var1)
!!       case('var2')
!!          call read_val(line, var2)
!!       case default
!!          cycle
!!    end select
!! enddo
!! close(444)
!! ~~~~~~
!
#include "parser_error_codes.h"
module parser
   use iso_fortran_env
   implicit none
   private
   !> A description of the keys we will find and their properties.
   type parser_var
      character(len=60):: varName='' !< The key name
      character(len=60):: section='' !< The section the key belongs to
      logical:: isSet=.false. !< Did we find and parsed a key with this name already?
      logical:: isRequired=.true. !< Is this key required?
   end type
   !> A description of the sections we will find and their properties.
   !! @since 2.6.3
   type parser_section
      character(len=60):: sectionName='' !< The section name
      logical:: isRequired=.true. !< Is this section required?
      logical:: isParsed=.false. !< Was a section with this name parsed?
   end type
   !> Reads a value
   interface read_val
      module procedure read_val_int
      module procedure read_val_real
      module procedure read_val_double
      module procedure read_val_logical
      module procedure read_val_string
   end interface

   public:: read_val, parse, parser_var, parser_section
   public:: checkEntryValidity, entryIndexForVarName

contains

   !********************************************************************
   !> Reads one line from the specified unit and outputs the key 
   !! name and possible values as a string.\n
   !! @par key is the key name.
   !! @par line is the rhs of the attribution.
   !! @par unit_in is the unit to read from
   !! @par error is an error code.
   !! \todo parse comments
   subroutine parse(unit_in, key, line, error)
      implicit none
      integer, intent(in):: unit_in
      character(len=*), intent(out):: key
      character(len=256), intent(out):: line
      integer, intent(out):: error
      character(len=256):: line2
      character:: first
      integer:: eq
      key   = ""
      error = 0
      do
         read(unit_in,'(A)',iostat=error) line
         if (error.ne.0) then
            if (error == iostat_end) then
               error = WARN_END_OF_FILE
            else
               error = ERR_CONFIG_UNREADABLE
            endif
            return
         endif
         line2 = trim(adjustl(line))
         first = line2(1:1)
         select case(first)
            case('#','*','')
               cycle
            case('[')
               eq = scan(line2,"]")
               key = 'SECTION'
               line = trim(adjustl(line2(2:eq-1)))
               exit
            case default
               eq = scan(line2,"=")
               key = trim(adjustl(line2(1:eq-1)))
               line = adjustl(line2(eq+1:))
               exit
         endselect
      enddo
   end subroutine parse

   !*********************************************************************
   !> Checks whether a key is valid according to a
   !! list of keys and a list of sections provided.\n
   !! Secion headers are always considered invalid keys.\n
   !! When a key is valid, the \a isSet attribute is set to .true.
   subroutine checkEntryValidity(key, line, sections, keys,  currentSection, isValidKey, isValidSection, error)
      implicit none
      integer, intent(out):: error
      logical, intent(out):: isValidKey
      logical, intent(inout):: isValidSection
      character(len=*), intent(inout):: currentSection
      character(len=*), intent(in):: key, line
      type(parser_var), intent(inout):: keys(:)
      type(parser_section), intent(in):: sections(:)
      integer:: i, nSections, nKeys
      error        = 0
      isValidKey = .false.
      nSections    = size(sections)
      nKeys        = size(keys) 

      ! We have a section header. Get its name and check it is valid.
      if (trim(key)=='SECTION') then
         call read_val(line, currentSection)
         isValidSection=.false.
         ! Check that the section is known.
         do i=1, nSections
            if (trim(currentSection) == trim(sections(i)%sectionName)) then
                  isValidSection=.true.
               exit
            endif
         enddo
         if (.not.isValidSection) then
            Write(*,*) 'Found an unknown section: ', trim(currentSection)
            error = ERR_INVALID_SECTION
            return
         endif
      else ! We have a key entry. Check that it is valid.
         entries:do i=1, nKeys
            if ( trim(keys(i)%varName) == trim(key)) then
               if ( trim(keys(i)%section) == trim(currentSection)) then
                  isValidKey  = .true.
                  keys(i)%isSet = .true.
                  exit
               else
                  Write(*,*) 'Entry ',trim(key), ' filed under the wrong section: ', trim(currentSection)
                  Write(*,*) 'Entry ',trim(key), ' should have been under: ', trim(keys(i)%section)
                  error = ERR_WRONG_ENTRY_IN_SECTION
                  return
               endif
            endif
         enddo entries
         if (.not.isValidKey) then
            Write(*,*) 'Found an invalid entry: ', trim(key)
            error = ERR_INVALID_ENTRY
            return
         endif
      endif
   end subroutine

   !*************************************************************
   !> Finds the key entry index from its name in a given list.
   pure integer function entryIndexForVarName(keyList,key)
      implicit none
      character(len=*), intent(in):: key
      type(parser_var), intent(in):: keyList(:)
      integer:: i, nKeys
      nKeys = size(keyList)
      do i=1, nKeys
         if (trim(keyList(i)%varName) == trim(key)) then
            entryIndexForVarName = i
            exit
         endif
      enddo
   end function

   !***************************************************
   pure subroutine read_val_int(line, value)
      implicit none
      character(len=256), intent(in):: line
      integer, intent(out):: value
      read(line,*) value
   end subroutine read_val_int
   
   !***************************************************
   pure subroutine read_val_real(line, value)
      implicit none
      character(len=256), intent(in):: line
      real, intent(out):: value
      read(line,*) value
   end subroutine read_val_real
   
   !***************************************************
   pure subroutine read_val_double(line, value)
      implicit none
      character(len=256), intent(in):: line
      double precision, intent(out):: value
      read(line,*) value
   end subroutine read_val_double
   
   !***************************************************
   pure subroutine read_val_logical(line, value)
      implicit none
      character(len=256), intent(in):: line
      logical, intent(out):: value
      read(line,*) value
   end subroutine read_val_logical
   
   !***************************************************
   pure subroutine read_val_string(line, value)
      implicit none
      character(len=256), intent(in):: line
      character(len=*), intent(out):: value
      integer:: spc
      spc = scan(line," ")
      value = line(:spc-1)
   end subroutine read_val_string

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
