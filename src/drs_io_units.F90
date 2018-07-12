! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2014
!> Manages the I/O units of DRS
module drs_io_units
   implicit none
   save
   ! Flow
   integer, parameter:: unit_ek    = 11
   integer, parameter:: unit_ur    = 12
   integer, parameter:: unit_uzon  = 13
   integer, parameter:: unit_koeu  = 14
   integer, parameter:: unit_uaz   = 15
   integer, parameter:: unit_u_mid = 16
   integer, parameter:: unit_am    = 17
   ! Heat
   integer, parameter:: unit_nu    = 21
   integer, parameter:: unit_adv   = 22
   integer, parameter:: unit_t     = 23
   ! Magnetic Field
   integer, parameter:: unit_eb    = 31
   integer, parameter:: unit_koeb  = 32
   integer, parameter:: unit_dissu = 33
   integer, parameter:: unit_dissB = 34
   ! Spectra
   integer, parameter:: unit_mspec = 41
   integer, parameter:: unit_lspec = 42
   integer, parameter:: unit_nspec = 43
   ! Kinematic dynamos
   integer, parameter:: unit_evp   = 51 
   integer, parameter:: unit_evt   = 52 
   ! Time
   integer, parameter:: unit_cfl   = 88

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
