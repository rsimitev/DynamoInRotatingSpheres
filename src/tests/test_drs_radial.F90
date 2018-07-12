! Copyright  L. Silva (lacsilva@gmail.com), 2015
program test_drs_radial
   use drs_dims
   use drs_params
   use drs_radial
   use drs_fftw3
   implicit none
   logical:: pass
   integer:: error, i
   character(len=7):: test

   eta = 0.4
   Nr  = 21
   Nr_s  = 21
   Nt  = 3
   Np  = 7
   pass=.true.

   call get_command_argument(1,test)

   call drs_params_init()
   call drs_dims_init(error)
   call drs_fftw3_init(Nr, Nt, Np)
   call drs_radial_init(eta)

   select case(trim(test))
   case('col')
      pass = test_radial_colocation_points()
   case('dr1')
      ! Test if the derivative of r**2 is 2*r.
      pass = test_radial_derivative_r2r()
   case('1D_r2r')
      pass = test_radial_dr_ddr_1D_r2r()
   case default
      pass = .false.
   end select

   if(pass) Write(*,*) 'Pass.'
contains
   !-----------------------------------------------------------------------------
   !> Tests the correctnes of the collocation points
   logical function test_radial_colocation_points()
      implicit none
      double precision:: my_ri, my_ro, my_r
      my_ri = eta/(1.0d0-eta)
      my_ro = 1.0d0/(1.0d0-eta)
      test_radial_colocation_points=.true.

      ! Test the end points
      if (dabs(rcoll(1)-my_ro).gt.1.0d-10) then
         Write(*,*) 'rcoll(1)=', rcoll(1),' should be ', my_ro
         test_radial_colocation_points=.false.
      endif
      if (dabs(rcoll(Nr)-my_ri).gt.1.0d-10) then
         Write(*,*) 'rcoll(1)=', rcoll(Nr),' should be ', my_ri
         test_radial_colocation_points=.false.
      endif

      ! Test two points in the middle, i=7, i=15
      i = 7
      my_r = my_ri + 0.5d0*(cos(pi*(i-1)/(dble(Nr)-1.d0))+1)
      if (dabs(rcoll(i)-my_r).gt.1.0d-10) then
         Write(*,*) 'rcoll(7)=', rcoll(i),' should be ', my_r
         test_radial_colocation_points=.false.
      endif
      i = 15
      my_r = my_ri + 0.5d0*(cos(pi*(i-1)/(dble(Nr)-1.d0))+1)
      if (dabs(rcoll(i)-my_r).gt.1.0d-10) then
         Write(*,*) 'rcoll(15)=', rcoll(i),' should be ', my_r
         test_radial_colocation_points=.false.
      endif
   end function
   
   !-----------------------------------------------------------------------------
   !> Tests if the derivative of r**2 is 2*r.
   logical function test_radial_derivative_r2r()
      implicit none
      double precision:: aux2(Nr)
      test_radial_derivative_r2r=.true.
      aux2 = radial_derivative_r2r(rcoll2)
      do i=1, Nr
         if(dabs(aux2(i)-2.0d0*rcoll(i)).gt.1.0d-10) then
            Write(*,*) rcoll(i), rcoll2(i), aux2(i)
            test_radial_derivative_r2r=.false.
         endif
      enddo
   end function

   !-----------------------------------------------------------------------------
   !> Tests if the first derivative of r**2 is 2*r and the second is 2.
   logical function test_radial_dr_ddr_1D_r2r()
      implicit none
      double precision:: aux1(Nr), aux2(Nr)
      test_radial_dr_ddr_1D_r2r=.true.
      call radial_dr_ddr_1D_r2r(rcoll2, aux1, aux2)
      do i=1, Nr
         if(dabs(aux1(i)-2.0d0*rcoll(i)).gt.1.0d-10) then
            Write(*,*) rcoll(i), rcoll2(i), aux1(i)
            test_radial_dr_ddr_1D_r2r=.false.
         endif
         if(dabs(aux2(i)-2.0d0).gt.1.0d-10) then
            Write(*,*) rcoll(i), aux1(i)
            test_radial_dr_ddr_1D_r2r=.false.
         endif
      enddo
   end function

end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
