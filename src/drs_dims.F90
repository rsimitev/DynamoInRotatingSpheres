! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Provides variables to store the real space and spectral space dimensions of the problem.
module drs_dims
   use drs_error_codes
   implicit none
   save

   integer, target::  usr_dims(8)
   integer:: Nr !< Number of radial points.
   integer:: Nt !< Number of meridional points.
   integer:: Np !< Number of azimuthal points.
   integer:: Nr_s !< Highest index for the polynomials in the radial direction.
   integer:: Nt_s !< Number of spherical harmonic degrees to use, including 0.
   integer:: Np_s !< Number of spherical harmonic orders (positive, negative and zero) to use.
   integer:: lsymm !< Equatorial symmetry.
   integer:: m0 !< Axial symmetry to use.

contains
   !< Packs all the dimensions into a neat array for distribution through mpi.
   subroutine drs_dims_init(error)
      implicit none
      integer, intent(out):: error
      
      error = 0 
      usr_dims(1:8) = (/Nr, Nt, Np, Nr_s, Nt_s, Np_s, lsymm, m0/)
   end subroutine drs_dims_init
   
   !****************************************************************
   !>     Checks consistency of input parameters.
   !****************************************************************
   subroutine check_dims(error)
      implicit none
      integer, intent(out):: error

      if(mod(Np_s,2)==0) then
         Write(*,*) 'Np_s should be an odd number. Correcting this.'
         Np_s = Np_s + 1
      endif
      if(mod(Np,2)==0) then
         Write(*,*) 'Np should be an odd number. Correcting this.'
         Np = Np + 1
      endif
      if(Nr_s.gt.Nr) then
         Write(*,*) 'Radial number of points should be at least equal to the number of modes.'
         Write(*,*) 'Adjusting Nr to ', Nr_s,'.'
         Nr = Nr_s
      endif
      if(Nt_s.gt.Nt) then
         Write(*,*) 'Meridional number of points should be at least equal to the number of modes.'
         Write(*,*) 'Adjusting Nt to ', Nt_s,'.'
         Nt = Nt_s
      endif
      if(Np_s.gt.Np) then
         Write(*,*) 'Azimuthal number of points should be at least equal to the number of modes.'
         Write(*,*) 'Adjusting Np to ', Np_s,'.'
         Np = Np_s
      endif

      if(m0*(Np_s/2).gt.Nt_s) Write(*,*) 'Warning: m0*Np_s/2>Nt_s: Some m''s with no l''s!'
      
      if( mod(Nr-1,2).ne.0 .and. mod(Nr-1,3).ne.0 ) then
         Write(*,*) 'Nr-1 should be power of 2 or 3.'
         error = ERR_NR1_NOT_POW2
      endif
      if( Nt.lt.Np/2 .or. Nt_s.lt.Np_s/2) then
         Write(*,*) 'Nt must be at least Np/2.'
         error = ERR_NT_TOO_SMALL
      endif
    end subroutine check_dims


end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
