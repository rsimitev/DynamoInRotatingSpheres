! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
module drs_legendre
#include "drsDefs.F90"
   use drs_dims
   use drs_mpi
   use drs_fftw3
   implicit none
   save
   double precision, allocatable:: legendre(:,:,:) !< The unnormalised Legendre polynomials
   double precision, allocatable:: leg_neg (:,:,:) !< The unnormalised Legendre polynomials for negative m multiplied by the integration factors.
   double precision, allocatable:: dleg    (:,:,:) !< d Plm(cos(theta))/d theta
   double precision, allocatable:: leg_sin (:,:,:) !< Plm/sin(theta)
   double precision, allocatable:: plmfac  (:,:) !< sqrt( (l+m)!/(l-m)!/(2l+1) ) = sqrt( (l-m+1)*(l-m+2)*...*((l+m)/(2l+1) )
   double precision, allocatable, target:: costheta(:) !< Gauss-Legendre integration points
   double precision, allocatable, target:: sintheta(:) !< Gauss-Legendre integration points
   double precision, allocatable:: w(:) !< Gauss-Legendre integration weights
   integer, allocatable:: llp1(:) !< Table of \f$l(l+1)\f$.

   double precision, private, parameter:: pi  = 3.141592653589793d0
   ! Interface copyed from SHTOOLS module
   interface
      subroutine PlmBar_d1(p, dp, lmax, z, csphase, cnorm)
         integer, intent(in) ::   lmax
         real*8, intent(out) ::   p(:), dp(:)
         real*8, intent(in) ::   z
         integer, intent(in), optional :: csphase, cnorm
      end subroutine PlmBar_d1
      subroutine PLegendreA_d1(p, dp, lmax, z, csphase)
         integer, intent(in) ::   lmax
         real*8, intent(out) ::   p(:), dp(:)
         real*8, intent(in) ::   z
         integer, intent(in), optional :: csphase
      end subroutine PLegendreA_d1
      integer function PlmIndex(l,m)
         integer, intent(in)   :: l, m
      end function PlmIndex
   end interface

contains

   !------------------------------------------------------------------
   !> Allocates module arrays legendre, leg_neg, dleg, plmfac, costheta
   !! sintheta, w and llp1.
   subroutine drs_legendre_allocation()
      implicit none
      integer:: lmax
      lmax = min(Nt_s+2,Nt)
      allocate(legendre( 0:Nt, 0:lmax, 0:blk_ps_size(mpi_rank)))
      allocate(leg_neg ( 0:lmax, 0:Nt, 0:blk_ps_size(mpi_rank)))
      allocate(dleg    ( 0:Nt, 0:lmax, 0:blk_ps_size(mpi_rank)))
      allocate(leg_sin ( 0:Nt, 0:lmax, 0:blk_ps_size(mpi_rank)))
      allocate(plmfac  ( 0:Nt_s+2, 0:Nt_s+2))
      allocate(costheta(0:Nt))
      allocate(sintheta(0:Nt))
      allocate(w(0:Nt))
      allocate(llp1(0:Nt_s))
   end subroutine

   !------------------------------------------------------------------
   !> Initialize the Legendre associated Polynomials and the
   !! Gauss-Legendre co-location points.
   subroutine drs_legendre_init()
      implicit none
      integer:: l
      forall(l=0:Nt_s) llp1(l) = l*(l+1)
      call initNormalization(Normalised, Nt_s, plmfac)
      !-- calculate locations costheta(0:Nt) and weights w(0:Nt)
      !-- for the gauss-legendre integration:
      !-- costheta(0)   = -1.d0 (theta=pi)
      !-- costheta(Nt)  =  1.d0 (theta=0)
      call gauleg(-1.0d0, 1.0d0, costheta(0:Nt), w(0:Nt), Nt+1)
      sintheta = dsqrt(1.0d0-costheta**2)
      !-- logitudial fourier-transformation:
      call legendre_init_new()
   end subroutine drs_legendre_init

   !------------------------------------------------------------------
   !> Computes the normalization factors for the the Legendre associated
   !! polynomials.
   subroutine initNormalization(normType, lmax, norms)
      implicit none
      double precision, intent(out):: norms(0:lmax+2, 0:lmax+2)
      !> normalization types:\n
      !!  - Normalised for normalization to 2.\n
      !!  - UnNormalized for no normalization, that is, \f\[ (P_l^m)^2  = 2*\frac{(l+m)!}{(l-m)!(2l+1)} \f\].
      integer, intent(in):: normType, lmax
      integer:: j, m, l
      double precision:: fac1, fac2
      !-- 14.1.97 M.A.: normalization factors for Plms:
      !-- 24.02.97 M.A.: bug fix for m=0.
      if (normType==Normalised) then
         norms = 0.d0
         do l=0, lmax+2
            do m=0, l
               norms(l,m) = dsqrt(2.0d0)
            enddo
         enddo
      else
         norms = 0.d0
         do l=0, lmax+2
            norms(l,0) = dsqrt(2.0d0)/dsqrt(dble(2*l+1))
            do m=1,l
               fac1 = 1.d0
               fac2 = 1.d0
               !-- sqrt( 2*(l+m)!/(l-m)!/(2l+1) ) = sqrt( 2*(l-m+1)*(l-m+2)*...*(l+m)/(2l+1) )
               !-- splitted in two parts to avoid floating overflow.
               do j=l-m+1,l+m,2
                  fac1 = fac1*dble(j)
                  fac2 = fac2*dble(j+1)
               enddo
               norms(l,m) = dsqrt(2.0d0*fac1)*dsqrt(fac2)/dsqrt(dble(2*l+1))
            enddo
         enddo
      endif
   end subroutine

   !------------------------------------------------------------------
   !> Initializes the tables of Associated Legendre Polynomials.
   subroutine legendre_init_new()
      implicit none
      double precision:: tmp_leg((Nt_s+1)*(Nt_s+2)/2)
      double precision:: tmp_dleg((Nt_s+1)*(Nt_s+2)/2)
      double precision:: fac
      integer:: i,j,k,l,lmax
      integer:: m
      legendre = 0.0d0
      dleg     = 0.0d0
      leg_neg  = 0.0d0
      leg_sin  = 0.0d0
      lmax     = min(Nt_s+2,Nt)
      do i=0, Nt
         call PlmBar_d1(tmp_leg, tmp_dleg, Nt_s, costheta(i), -1, 1)
         do j=1, blk_ps_size(mpi_rank)
            m = blk_ts_start(j)
            do l = m, Nt_s
               k = PlmIndex ( l, m )
               legendre(i, l, j) = tmp_leg(k)
               fac = (1.0d0/plmfac(l,m))**2
               if(m.le.lmax) leg_neg(l, i, j)  = fac*w(i)*tmp_leg(k)
               dleg(i, l, j )    = -tmp_dleg(k)*sintheta(i)
               leg_sin(i, l, j ) =   tmp_leg(k)/sintheta(i)
            enddo
         enddo
      enddo
   end subroutine legendre_init_new

   !------------------------------------------------------------------
   !> Computes the Guass-Legendre quadrature points and weights.
   pure subroutine gauleg(x1,x2,x,w,n)
      implicit none
      integer, intent(in):: n
      double precision, intent(in):: x1, x2
      double precision, intent(out):: x(n), w(n)
      double precision, parameter:: eps=1.0d-15
      double precision:: pp, p1, p2, p3, z, z1
      double precision:: xm, xl
      integer:: i,j, m
      x = 0.0d0
      w = 0.0d0
      m=(n+1)/2         ! The roots are symmetric in the interval, so we
                  ! only have to find half of them.
      xm = (x2 + x1) / 2.0d0    ! midpoint of integration
      xl = (x2 - x1) / 2.0d0    ! Scaling factor between interval of integration, and that of
                                ! the -1 to 1 interval for the Gauss-Legendre interval
      ! Compute roots and weights
      do i=1,m
         z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))    ! Approximation for the ith root
         ! Find the true value using newtons method
         do
             p1=1.0d0
             p2=0.0d0
             do j=1, n     ! determine the legendre polynomial evaluated at z (p1) using
                           ! recurrence relationships
                 p3=p2
                 p2=p1
                 p1=(dble(2*j-1)*z*p2-dble(j-1)*p3)/dble(j)
             enddo

             pp=dble(n)*(z*p1-p2)/(z*z-1.0d0)    ! This is the derivative of the legendre polynomial
                                                 ! using recurrence relationships
             z1=z                    ! This is newtons method here
             z=z1-p1/pp

             if (abs(z-z1) <= eps) exit
         enddo
         x(i)     = xm - xl*z
         x(n+1-i) = xm + xl*z
         w(i)     = 2.0d0*xl/((1.d0-z**2)*pp**2)
         w(n+1-i) = w(i)
      enddo
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
