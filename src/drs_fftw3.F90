! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> This module abstracts the computation of Fourier and cosinus transforms.
module drs_fftw3
   implicit none
   save
   !> It makes use of fftw3.
   include  'fftw3.f'
   integer*8:: plan_r, plan_pf, plan_pb
   double precision, allocatable:: in_p(:,:), inout_r(:)
   double complex, allocatable::  out_p(:,:)
   integer, private:: drs_fftw3_Nr, drs_fftw3_Np, drs_fftw3_Nt
contains
   !> Initialises all the fftw3 plans for forward and backward Fourier and cosinus transforms.
   subroutine drs_fftw3_init(Nr, Nt, Np)
      implicit none
      integer,intent(in):: Nr !< Number of points in real space for the radial cosinus transforms.
      integer,intent(in):: Np !< Number of points in real space for the azimuthal Fourier transforms.
      integer,intent(in):: Nt !< Perform this many azimuthal Fourier transforms at a time.
      drs_fftw3_Nr   = Nr
      drs_fftw3_Np   = Np
      drs_fftw3_Nt   = Nt
      allocate(inout_r(Nr+1))
      allocate(in_p(Np,0:Nt-1), out_p(Np/2+1,0:Nt-1))
      call dfftw_plan_r2r_1d(plan_r, Nr, inout_r, inout_r, FFTW_REDFT00, FFTW_MEASURE)
      call dfftw_plan_many_dft_r2c(plan_pf, 1, (/Np/), Nt, &
                             in_p,  (/Np/),     1, Np, &
                             out_p, (/Np/2+1/), 1, Np/2+1, &
                             FFTW_MEASURE)
      call dfftw_plan_many_dft_c2r(plan_pb, 1, (/Np/), Nt, &
                             out_p, (/Np/2+1/), 1, Np/2+1, &
                             in_p,  (/Np/),     1, Np, &
                             FFTW_MEASURE)
   end subroutine

   !> Destroies the plans
   subroutine drs_fftw3_cleanup()
      implicit none
      call dfftw_destroy_plan(plan_r)
      call dfftw_destroy_plan(plan_pf)
      call dfftw_destroy_plan(plan_pb)
   end subroutine

   !> The forward real to spectral DFT.
   subroutine dft_forward(input,output)
      implicit none
      double precision, intent(in)::  input (0:drs_fftw3_Nt-1,drs_fftw3_Np)
      double precision, intent(out):: output(0:drs_fftw3_Nt-1,drs_fftw3_Np)
      double complex:: aux(0:drs_fftw3_Nt-1,drs_fftw3_Np/2+1)
      integer:: m, Nt1

      Nt1 = drs_fftw3_Nt-1

      in_p(1:drs_fftw3_Np, 0:Nt1) = transpose(input(0:Nt1, 1:drs_fftw3_Np))
      call dfftw_execute_dft_r2c(plan_pf, in_p, out_p)
      aux = transpose(out_p)
      output(0:Nt1,1) = dble(aux(0:Nt1,1))/drs_fftw3_Np
      forall (m=2:drs_fftw3_Np/2)
         output(0:Nt1,2*m-2) = dble(aux(0:Nt1,m))/drs_fftw3_Np
         output(0:Nt1,2*m-1) = imag(aux(0:Nt1,m))/drs_fftw3_Np
      endforall
      m = drs_fftw3_Np/2 + 1
      output(0:Nt1,2*m-2) = dble(aux(0:Nt1,m))/drs_fftw3_Np
      if (2*m-1.eq.drs_fftw3_Np) then
         output(0:Nt1,2*m-1) = imag(aux(0:Nt1,m))/drs_fftw3_Np
      endif
   end subroutine

   !> The backward, spectral to real DFT.
   subroutine dft_backward(input,output)
      implicit none
      double precision, intent(in)::  input (0:drs_fftw3_Nt-1,drs_fftw3_Np)
      double precision, intent(out):: output(0:drs_fftw3_Nt-1,drs_fftw3_Np)
      integer:: m, Nt1

      Nt1 = drs_fftw3_Nt-1

      out_p(1,0:Nt1) = dcmplx(input(0:Nt1,1),0.0d0)
      forall (m=2:drs_fftw3_Np/2)
         out_p(m,0:Nt1) = dcmplx(input(0:Nt1,2*m-2),input(0:Nt1,2*m-1))
      endforall
      m = drs_fftw3_Np/2 + 1
      if (2*m-1.eq.drs_fftw3_Np) then
         out_p(m,0:Nt1) = dcmplx(input(0:Nt1,2*m-2),input(0:Nt1,2*m-1))
      else
         out_p(m,0:Nt1) = dcmplx(input(0:Nt1,2*m-2),0.0d0)
      endif
      call dfftw_execute_dft_c2r(plan_pb, out_p, in_p)
      output(0:Nt1, 1:drs_fftw3_Np) = transpose(in_p(1:drs_fftw3_Np, 0:Nt1))
   end subroutine

   !> The forward real to spectral cosinus transform
   subroutine cos_r2r_1_r2n(input)
      implicit none
      double precision, intent(inout):: input(drs_fftw3_Nr)
      call dfftw_execute_r2r(plan_r, input, input)
      input(:) = input(:)/(2*(drs_fftw3_Nr-1))
      input(2:drs_fftw3_Nr-1) = input(2:drs_fftw3_Nr-1)*2
   end subroutine
   !> The backward spectral to real cosinus transform
   subroutine cos_r2r_1_n2r(input)
      implicit none
      double precision, intent(inout):: input(drs_fftw3_Nr)
      input(2:drs_fftw3_Nr-1) = input(2:drs_fftw3_Nr-1)/2
      call dfftw_execute_r2r(plan_r, input, input)
   end subroutine

   !> Given a field f1 described at Nr1 points in an interval, 
   !! \a remesh outputs the same field, in the same interval at 
   !! Nr2 points as f2.
   subroutine remesh(Nr1, f1, Nr2, f2)
      implicit none
      integer, intent(in):: Nr1, Nr2 !< The input and output number of points.
      double precision, intent(in)::  f1(Nr1) !< The input field.
      double precision, intent(out):: f2(Nr2) !< The output field.
      double precision:: radarr1(Nr1), radarr2(Nr2)
      integer*8:: plan_Nr1, plan_Nr2
      integer:: i
      if(Nr1 .eq. Nr2) then
         f2 = f1
         return
      endif

      call dfftw_plan_r2r_1d(plan_Nr1, Nr1, radarr1, radarr1, FFTW_REDFT00, FFTW_ESTIMATE)
      call dfftw_plan_r2r_1d(plan_Nr2, Nr2, radarr2, radarr2, FFTW_REDFT00, FFTW_ESTIMATE)
      radarr1 = f1
      radarr2 = 0.0d0

      call dfftw_execute(plan_Nr1)
      radarr1 = radarr1/(Nr1-1)/2

      do i=1, min(Nr1,Nr2)
         radarr2(i) = radarr1(i)
      enddo

      !call dfftw_execute_r2r(plan_Nr2, radarr2, radarr2)
      call dfftw_execute(plan_Nr2)
      f2 = radarr2
      call dfftw_destroy_plan(plan_Nr1)
      call dfftw_destroy_plan(plan_Nr2)
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
