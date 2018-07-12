! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> This module implements the radial domain and operations in it.
module drs_radial
#include "drsDefs.F90"
   use drs_dims
   use drs_mpi
   use drs_Chebyshev
   use drs_fftw3
   implicit none
   save
   private
   double precision, allocatable:: rcoll(:) !< Radial collocation points for Chebychev polynomials
   double precision, allocatable:: rcoll2(:) !< Squares of radial collocation points for Chebychev polynomials
   double precision, allocatable:: drcoll(:) !< Differences for radial collocation points for Chebychev polynomials
   double precision, allocatable, public:: poly(:,:), poly_dr(:,:), poly_ddr(:,:)
   !> Index of the boundaries: b(1)=inner boundary; b(2)=outer boundary
   integer::b(2)

   !> Value of the internal radius.
   double precision:: ri
   !> Value of the external radius.
   double precision:: ro

   public:: rcoll, rcoll2, drcoll, ri, ro, b
   public:: drs_radial_init
   public:: radial_dr_ddr_1D_n2r, radial_dr_ddr_1D_r2r
   public:: radial_derivative_r2r
   public:: radial_dr_ddr_3D_n2r, radial_dr_ddr_3D_r2r, radial_transform_3D_n2r

contains

   !-----------------------------------------------------------------------------
   !> Allocates the required memory for the caches of the radial domain and
   !! Chebyshev polynomials.
   subroutine drs_radial_allocate()
      implicit none
      allocate( rcoll(Nr), rcoll2(Nr), drcoll(Nr) )
      allocate( poly(Nr,Nr), poly_dr(Nr,Nr), poly_ddr(Nr,Nr) )
   end subroutine

   !-----------------------------------------------------------------------------
   !> Initializes the radial domain and the Chebyshev polynomials and their
   !! derivatives.
   subroutine drs_radial_init(riro)
      implicit none
      double precision, intent(in):: riro
      integer:: i
      call drs_radial_allocate()
      call Chebyshev_init(Nr, Nr_s)
      poly     = Chebyshev
      poly_dr  = 2*Chebyshev_dx
      poly_ddr = 4*Chebyshev_ddx
      !-- radial collocation points rcoll(Nr) for Chebychev polynomials:
      do i=1, Nr
         rcoll(i)  = riro/(1.d0-riro) + 0.5d0*(Cheb_x(i)+1)
         rcoll2(i) = rcoll(i)**2
         if (i.gt.1) drcoll(i-1) = rcoll(i-1) - rcoll(i)
      enddo
      drcoll(Nr) = drcoll(Nr-1)
      ! index of the boundaries: b(1)=inner boundary; b(2)=outer boundary
      b = (/Nr, 1/)
      ri = rcoll(b(1))
      ro = rcoll(b(2))
   end subroutine drs_radial_init

   !-----------------------------------------------------------------------------
   !> Evaluates the function in rear space
   subroutine radial_transform_3D_n2r(t)
      implicit none
      integer:: j,l
      double precision, intent(inout):: t(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      jl_do(j,l)
        call Chebyshev_n2x(t(l,j,1:Nr))
      jl_enddo
   end subroutine


   !-----------------------------------------------------------------------------
   !> Computes the coefficients of the radial derivative of \a radarr
   pure subroutine dradial_n2n(radarr, deriv)
      implicit none
      double precision, intent(in):: radarr(Nr)
      double precision, intent(out):: deriv(Nr)
      call Cheb_compute_dx_n2n(radarr, deriv)
   end subroutine

   !-----------------------------------------------------------------------------
   !> Returns the first derivative of \a radarr.
   !! \a radarr is supposed to be given in direct space, derivative is returned
   !! in direct space.
   !! \since 1.6.6
   function radial_derivative_r2r(radarr) result(deriv)
      implicit none
      double precision, intent(in):: radarr(Nr)
      double precision:: deriv(Nr)
      double precision:: aux(Nr)
      aux = radarr
      ! transformation of input: (lmr) -> (lmn)
      call Chebyshev_x2n(aux)
      ! Compute the derivatives in spectral space
      call dradial_n2n(aux, deriv)
      ! transformation of result: (lmn) -> (lmr)
      call Chebyshev_n2x(deriv)
      ! factor of 2 due to the mapping from the radial
      ! coordinate r to the Chebyshev coordinate x, where x runs from -1 to 1.
      ! The interrelation is r=eta/(1-eta)+0.5(x+1)  see the def of rcoll in the
      ! initialization routine. Obviously, d/dr = 2*d/dx, that's this factor of 2.
      deriv = 2.0d0*deriv
   end function radial_derivative_r2r

   !-----------------------------------------------------------------------------
   !> Returns first radial derivative in \a t1, second derivative in \a t2.
   !! Input \a t is supposed to be given in spectral space.
   !! On output both \a t and its derivatives are returned
   !! in real space.
   !! \since 1.6.6
   subroutine radial_dr_ddr_1D_n2r(t, t1, t2)
      implicit none
      double precision, intent(inout):: t(Nr)
      double precision, intent(out):: t1(Nr), t2(Nr)
      call Cheb_compute_dx_ddx_n2x(t, t1, t2)
      call Chebyshev_n2x(t)
      !> A factor of 2 for each derivative is due to the mapping from the radial
      !! coordinate r to the Chebyshev coordinate x, where x runs from -1 to 1.
      !! The interrelation is r=eta/(1-eta)+0.5(x+1)  see the def. of rcoll in the
      !! initialization routine. Obviously, d/dr = 2*d/dx.
      t1(1:Nr) = 2*t1(1:Nr)
      t2(1:Nr) = 4*t2(1:Nr)
   end subroutine radial_dr_ddr_1D_n2r

   !-----------------------------------------------------------------------------
   !> Returns first radial derivative in \a t1, second derivative in \a t2.
   !! Input \a t is supposed to be given in real space.
   !! On output both derivatives are returned in real space.
   !! \since 1.6.6
   subroutine radial_dr_ddr_1D_r2r(t, t1, t2)
      implicit none
      double precision, intent(in):: t(Nr)
      double precision, intent(out):: t1(Nr), t2(Nr)
      double precision:: aux(Nr)
      aux=t
      ! transformation of input: (lmr) -> (lmn)
      call Chebyshev_x2n(aux)
      ! Compute the derivatives with respect to r and output in lmr
      call radial_dr_ddr_1D_n2r(aux, t1, t2)
   end subroutine

   !----------------------------------------------------------------------
   !> Calculates first and second radial derivatives of 3D-array in lmr space.
   !! Includes dealiasing in n.
   !! @param t The original field in lmr space.
   !! @param t1 The first radial drivative in lmr space.
   !! @param t2 The second radial drivative in lmr space.
   !! \since 1.6.6
   subroutine radial_dr_ddr_3D_r2r(t,t1,t2)
      implicit none
      double precision, intent(in)::   t(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision, intent(out):: t1(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision, intent(out):: t2(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      integer:: j,l,m
      t1 = 0.d0
      t2 = 0.d0

      do j=1, blk_ps_size(mpi_rank)
         m = blk_ts_start(j)
         do l=m, Nt_s
            call Cheb_compute_dx_ddx_x2x(t(l,j,1:Nr), t1(l,j,1:Nr), t2(l,j,1:Nr))
            ! This factor of 2 is due to the mapping from the radial
            ! coordinate r to the Chebyshev coordinate x, where x runs from -1 to 1.
            ! The interrelation is r=eta/(1-eta)+0.5(x+1)  see the def of rcoll in the
            ! initialization routine. Obviously, d/dr = 2*d/dx, that's this factor of 2.
            t1(l,j,1:Nr) = 2*t1(l,j,1:Nr)
            t2(l,j,1:Nr) = 4*t2(l,j,1:Nr)
         enddo  ! l=m,Nt_s
      enddo  ! j=1,nadl
   end subroutine radial_dr_ddr_3D_r2r

   !----------------------------------------------------------------------
   !> Calculates first and second radial derivatives of 3D-array in lmn space.
   !! includes dealiasing in n. Transforms the original field to lmr space.
   !! @param t0 The original field. lmn space on entry, lmr space on exit.
   !! @param t1 The first radial drivative in lmr space.
   !! @param t2 The second radial drivative in lmr space.
   !! \since 1.6.6
   subroutine radial_dr_ddr_3D_n2r(t0,t1,t2)
      implicit none
      double precision, intent(inout):: t0(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision, intent(out)::   t1(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision, intent(out)::   t2(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      integer:: j,l,m

      do j=1, blk_ps_size(mpi_rank)
         m = blk_ts_start(j)
         do l=m, Nt_s
            call Cheb_compute_dx_ddx_n2x(t0(l,j,1:Nr), t1(l,j,1:Nr), t2(l,j,1:Nr))
            call Chebyshev_n2x(t0(l,j,1:Nr))
            ! This factor of 2 is due to the mapping from the radial
            ! coordinate r to the Chebyshev coordinate x, where x runs from -1 to 1.
            ! The interrelation is r=eta/(1-eta)+0.5(x+1)  see the def of rcoll in the
            ! initialization routine. Obviously, d/dr = 2*d/dx, that's this factor of 2.
            t1(l,j,1:Nr) = 2*t1(l,j,1:Nr)
            t2(l,j,1:Nr) = 4*t2(l,j,1:Nr)
         enddo  ! l=m,Nt_s
      enddo  ! j=1,nadl
   end subroutine radial_dr_ddr_3D_n2r

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
