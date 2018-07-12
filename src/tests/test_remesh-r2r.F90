! Copyright  L. Silva (lacsilva@gmail.com), 2015
program test_drs_fftw
   use drs_fftw3
   implicit none
   integer, parameter::N=120, P=200
   double precision, dimension(N):: x1
   double precision, dimension(P):: z, x2, y2
   double precision, dimension(N):: y
   double precision:: dx1, dx2
   integer:: i

   dx1 = 3.14159265359d0/(N-1)
   dx2 = 3.14159265359d0/(P-1)
   do i=1, N
      x1(i) = dble(i-1)*dx1
   enddo
   do i=1, P
      x2(i) = dble(i-1)*dx2
   enddo
   !Construct the test function on both grids
   y  = test_func(x1)
   y2 = test_func(x2)

   ! Perform the remeshing
   call remesh(N, y, P, z)


   Write(*,*) ' X                        Expected                Obtained'
   do i=1, P
      Write(*,*) x2(i), y2(i),  z(i)
      if(abs(y2(i)).lt.tiny(0.0d0)) then
         if( (y2(i)-z(i))/y2(i) .gt. 1.0d-10 ) then
            Write(*,*) 'Discrepancy bigger than 1d-8% at i = ', i
            stop
         endif
      endif
   enddo
   ! If we got here, no problems were found
   Write(*,*) 'Pass.'
contains
   elemental double precision function test_func(x)
      implicit none
      double precision, intent(in):: x
      test_func = -0.1d0*0.5 + &
              4.0d0*dcos(x) + &
              1.0d0*dcos(2*x) - &
              0.5d0*dcos(3*x) + &
              2.2d0*dcos(7*x)
   end function
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
