! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
program test_drs_read_conf_v3
   use drs_io_par
   implicit none
   integer:: error

   call drs_conf_known_vars()
   call drs_read_conf_v3(error)
   
   if (error.ne.0) then
      Write(*,*) 'Fail, ', error
   else
      Write(*,*) 'Pass.'
   endif
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
