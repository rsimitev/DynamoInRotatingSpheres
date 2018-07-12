program drsVersion
#include "drsDefs.F90"
   implicit none

   Write(*,*) 'DRS Dynamo in Rotating Sphere ',VERSION
#ifdef STAR
   Write(*,*) '(using time dependent Rayleigh number: tau = ', Ra_tau, ')'
#endif
#ifdef COMP
   Write(*,*) '(Compiled with support for compositional convection.)'
#endif
end program


! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
