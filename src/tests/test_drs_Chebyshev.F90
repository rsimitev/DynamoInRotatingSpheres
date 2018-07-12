! Copyright  L. Silva (lacsilva@gmail.com), 2015
program test_drs_radial
   use drs_Chebyshev
   implicit none
   logical:: pass
   integer:: i, j, Nx
   double precision:: val(6)
   double precision, allocatable:: f(:)
   character(len=3):: test

   Nx=24
   pass=.true.

   call Chebyshev_init(Nx, Nx)
   call get_command_argument(1,test)

   if(test=='') pass=.false.
   Write(*,*) 'Test name', test
   if (test == 'che') then
      open(unit=50, file='Chebyshev.dat', status='OLD')
      ! Check the Chebyshev polynomials
      do i=1, Nx
         read(50, *) val
         do j=1, 5
            if (dabs(Chebyshev(i,j)-val(j+1)).gt.1.0d-10) then
               Write(*,*) 'Chebyshev(',i,',',j,')', Chebyshev(i,j),' should be ', val(j+1)
               pass=.false.
            endif
         enddo
      enddo
   elseif (test=='cdx') then
      ! Check the derivatives of the Chebyshev polynomials
      open(unit=50, file='Chebyshev_dx.dat', status='OLD')
      do i=1, Nx
         read(50, *) val
         do j=1, 5
            if (dabs(Chebyshev_dx(i,j)-val(j+1)).gt.1.0d-10) then
               Write(*,*) 'Chebyshev_dx(',i,',',j,')', Chebyshev_dx(i,j),' should be ', val(j+1)
               pass=.false.
            endif
         enddo
      enddo
   elseif (test=='rdt') then
      allocate(f(Nx))
      f=0.0d0
      f(3)=1.0d0
      call Chebyshev_n2x(f)
      call Chebyshev_x2n(f)
      do i=1, Nx
         Write(*,*) f(i)
      enddo
      if (dabs(f(3)-1.0d0).gt.1.0d-10) then
         pass=.false.
      endif
   endif

   if(pass) Write(*,*) 'Pass.'
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
