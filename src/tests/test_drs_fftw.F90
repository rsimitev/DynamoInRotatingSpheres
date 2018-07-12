! Copyright  L. Silva (lacsilva@gmail.com), 2015
program test_drs_fftw
   use drs_fftw3
   implicit none
   integer, parameter::N=192, M=3
   double precision, dimension(N):: x
   double precision, dimension(0:M,N):: z
   double precision, dimension(0:M,N):: y
   double precision, dimension(0:M,N):: y1
   double precision:: dx
   integer:: i,j
   character(len=8):: f
   z=0.0d0
   y1=0.0d0

   Write(f,'("(",I1,"D15.5)")') M+1

   dx = 4.0d0*3.14159265359d0/N
   forall(i=1:N)
      x(i) = (i-1)*dx
   endforall
   !Construct a few cosinus
   forall(i=1:N, j=0:M)
      y(j,i) = dcos((j+1)*x(i)) + 2.0d0
   endforall
   Write(*,*) '---Original series---'
   do i=1, N
      Write(*,f) (y(j,i), j=0, M)
      Write(10,f) (y(j,i), j=0, M)
   enddo
   ! Initialise the dft's
   call drs_fftw3_init(N, M+1, N)
   call dft_forward(y,z)

   Write(*,*) '---transformed series---'
   do i=2, N, 2
      Write(*,f) (z(j,i), j=0, M)
   enddo
   do i=3, N, 2
      Write(*,f) (z(j,i), j=0, M)
   enddo
   call dft_backward(z,y1)
   Write(*,*) '---roundtrip series---'
   dx = 0.0d0
   do i=1, N
      Write(*,f) (y1(j,i), j=0, M)
      do j=0, M
         dx = dx + (y(j,i)-y1(j,i))**2
      enddo
   enddo

   Write(*,*) '---r.m.s.---'
   Write(*,*) dsqrt(dx)
   if (dsqrt(dx).lt.1.0d-10) Write(*,*) 'Pass.'
   call drs_fftw3_cleanup()

end program
