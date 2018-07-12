! Copyright  L. Silva (lacsilva@gmail.com), 2015
program test_drs_fftw
   use drs_fftw3
   implicit none
   integer, parameter::N=21, M=1
   double precision, dimension(N):: x
   double precision, dimension(N):: y, y1, z
   double precision:: dx
   integer:: i
   logical:: pass

   pass = .true.
   dx = 3.14159265359d0/(N-1)
   forall(i=1:N)
      x(i) = (i-1)*dx
   endforall
   !Construct a few cosinus
   forall(i=1:N)
      y(i) = -0.1d0*1.0d0/2.0d0 + &
              4.0d0*dcos(x(i)) + &
              1.0d0*dcos(2*x(i)) - &
              0.5d0*dcos(3*x(i)) + &
              2.2d0*dcos((N-1)*x(i))
   endforall
   Write(*,*) '---original series---'
   do i=1, N
      Write(*,*) x(i), y(i)
   enddo
   z = y
   ! Initialise the dft's
   call drs_fftw3_init(N, M+1, N)
   call cos_r2r_1_r2n(z)

   Write(*,*) '---transformed series---'
   do i=1, N
      Write(*,*) i-1, z(i)
   enddo
   if (abs(z(1)+0.05d0).gt.1.0d-10) then
      Write(*,*) 'z(0) = ', z(1),'ref: -0.05'
      pass = .false.
   endif
   if (abs(z(2)-4.0d0).gt.1.0d-10) then
      Write(*,*) 'z(1) = ', z(2),'ref: 4.0d0'
      pass = .false.
   endif
   if (abs(z(3)-1.0d0).gt.1.0d-10) then
      Write(*,*) 'z(2) = ', z(3),'ref: 1.0d0'
      pass = .false.
   endif
   if (abs(z(4)+0.5d0).gt.1.0d-10) then
      Write(*,*) 'z(3) = ', z(4),'ref: -0.5d0'
      pass = .false.
   endif
   if (abs(z(N)-2.2d0).gt.1.0d-10) then
      Write(*,*) 'z(N-1) = ', z(N),'ref: 2.2d0'
      pass = .false.
   endif

   y1 = z
   call cos_r2r_1_n2r(y1)

   Write(*,*) '---roundtrip---'
   dx = 0.0d0
   do i=1, N
      Write(*,*) x(i), y1(i), y1(i)-y(i)
      dx = dx + (y1(i)-y(i))**2
   enddo
   if (dsqrt(dx).gt.1.0d-10) then
      pass = .false.
      Write(*,*) 'r.m.s. = ', dsqrt(dx)
   endif
   if(pass) Write(*,*) 'Pass.'
end program
