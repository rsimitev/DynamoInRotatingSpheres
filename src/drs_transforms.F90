! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Spherical transforms.
module drs_transforms
#include "drsDefs.F90"
   use drs_dims
   use drs_mpi
   use drs_fftw3
   use drs_radial
   use drs_legendre
   implicit none
   save

   double precision, private, parameter:: pi  = 3.141592653589793d0

contains

   !--------------------------------------------------------------------------
   !> This routine performs three steps:
   !!~~~~~~
   !!  1. transformation    tt_t(theta,phi) --> tt_t(theta,m)  transposed.
   !!  2. redistribution        dist(theta) --> dist(m)
   !!  3. transformation        tt(theta,m) --> t(l,m)      distrib. in m.
   !!
   !!  input: tt_t(0:(blk_t_max_size-1),Np) (transposed)
   !!         has to be in direct space, in phi-direction t(:,1:Np)
   !!         only one period is stored:
   !!         tt_t(.,j) means position phi=(2pi/m0)*(j-1)/Np.
   !!
   !! output:
   !!         t(0:drs_Nt,drs_Np) dims: (0:Nt_s,Np_s) dealiased in (lm)!
   !!         in m-dir.: real(m=0),real(m=m0),im(m=m0),real(m=2m0),im(m=2*m0),...
   !!         in l-dir.: includes normalization factors!
   !!
   !!~~~~~~
   subroutine ylmt(input, output)
      implicit none
      !-- input  field (transposed):
      double precision, intent(in):: input(0:blk_t_size(mpi_rank)-1,Np)
      double precision, intent(out):: output(0:Nt_s,blk_ps_size(mpi_rank))
      double precision:: aux(0:blk_t_size(mpi_rank)-1,Np)
      double precision:: aux2(0:Nt,blk_ps_size(mpi_rank))  !  2D workarr distrib(m)
      integer:: j

      output = 0.0d0
      ! First a Fourier transform (phi -> m, with normalization 1/Np)
      call dft_forward(input, aux)
      ! Dealiase in m
      if (Np_s+1 .le. Np) aux(:, Np_s+1:Np) = 0.d0
      ! The field is now in input(itheta,j)  (distrib. in theta)
      ! transposition: aux(itheta,j) distrib(theta) --> aux2 distrib(m):
      call transpos_theta2phi(aux, Np, aux2, Nt)
      ! The field is now in aux2(itheta,j=2*(m/m0))  (distrib. in m)
      ! Gauss integration aux2(theta) -> output(l)
      do j=1,blk_ps_size(mpi_rank)
         output(mm(j):Nt_s,j) = matmul( leg_neg(mm(j):Nt_s,0:Nt,j), aux2(0:Nt,j) )
      enddo
      !-- the result is now in output(l,j=2*(m/m0))  (distrib. in m)
   end subroutine ylmt

   !--------------------------------------------------------------------------
   !> This routine performs three steps:
   !!~~~~~~
   !!  1. transformation    tt_t(theta,phi) --> tt_t(theta,m)  transposed.
   !!  2. redistribution        dist(theta) --> dist(m)
   !!  3. transformation        tt(theta,m) --> t(l,m)      distrib. in m.
   !!
   !!  input: tt_t(0:(blk_t_max_size-1),Np) (transposed)
   !!         has to be in direct space, in phi-direction t(:,1:Np)
   !!         only one period is stored:
   !!         tt_t(.,j) means position phi=(2pi/m0)*(j-1)/Np.
   !!
   !! output:
   !!         t(0:drs_Nt,drs_Np) dims: (0:Nt_s,Np_s) dealiased in (lm)!
   !!         in m-dir.: real(m=0),real(m=m0),im(m=m0),real(m=2m0),im(m=2*m0),...
   !!         in l-dir.: includes normalization factors!
   !!
   !!~~~~~~
   subroutine ylmt_3D(input, output)
      implicit none
      !-- input  field (transposed):
      double precision, intent(in):: input(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: output(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision:: aux(0:blk_t_size(mpi_rank)-1,Np)
      double precision:: aux2(0:Nt,blk_ps_size(mpi_rank))  !  2D workarr distrib(m)
      integer:: j, i

      output = 0.0d0
      do i=1, Nr
         ! First a Fourier transform (phi -> m, with normalization 1/Np)
         call dft_forward(input(0:blk_t_size(mpi_rank)-1,1:Np,i), aux)
         ! Dealiase in m
         if (Np_s+1 .le. Np) aux(:, Np_s+1:Np) = 0.d0
         ! The field is now in input(itheta,j)  (distrib. in theta)
         ! transposition: aux(itheta,j) distrib(theta) --> aux2 distrib(m):
         call transpos_theta2phi(aux, Np, aux2, Nt)
         ! The field is now in aux2(itheta,j=2*(m/m0))  (distrib. in m)
         ! Gauss integration aux2(theta) -> output(l)
         do j=1,blk_ps_size(mpi_rank)
            output(mm(j):Nt_s,j,i) = matmul( leg_neg(mm(j):Nt_s,0:Nt,j), aux2(0:Nt,j) )
         enddo
         !-- the result is now in output(l,j=2*(m/m0))  (distrib. in m)
      enddo
   end subroutine ylmt_3D
   ! -----------------------------------------------------------------------------
   !--
   !!-- this routine performs three steps:
   !!--
   !!--  1. transformation        t(l,m) --> tt(theta,m)
   !!--  2. redistribution       dist(m) --> dist(theta)
   !!--  3. transformation   tt(theta,m) --> tt(theta,phi)
   !!--
   !!--
   !!--  input:
   !!--         t(0:Nt_s,Np_s)        ! (dealiased!).
   !!--         t(0:drs_Nt,drs_Np) ! physical dims.
   !!--         in m-dir.: real(m=0),real(m=m0),im(m=m0),real(m=2m0),im(m=2*m0),...
   !!--         in l-dir.: includes normalization factors!
   !!--
   !!-- output:
   !!--         tt_t(0:Ntl,Np)
   !!--         tt_t(0:(blk_t_max_size-1),Np)  !  physical dims.
   !!--
   !!--         in phi-direction tt_t(.,1:Np) only one period is stored.
   !!--         tt(.,j) means position phi=(2pi/m0)*(j-1)/Np.
   !--
   !-- 05/96 M.A.: improved order of indices: t(0:drs_Nt,drs_Np)
   !-- 06/96 M.A.: workarrays are local now
   !-- 09/96 M.A.: opt.: replaced big tranf. (m0*Np) by small (Np).
   !-- 10/96 M.A.: parallelized.
   !-- 24.10.96 M.A.: output field is transposed now.
   !-- 01/97 M.A.: DOACROSS is done in calling routine (rhs).
   subroutine ylmb(input,output)
      implicit none
      !-- input  field, (distrib. in m):
      double precision, intent(in):: input(0:Nt_s,1:blk_ps_size(mpi_rank))
      !-- output field, transposed:
      double precision, intent(out):: output(0:blk_t_size(mpi_rank)-1,Np)
      double precision:: aux(0:Nt,1:blk_ps_size(mpi_rank)) ! (2D, distrib. in m)
      integer:: j

      output = 0.0d0
      aux    = 0.0d0
      !-- first add up all Legendre polynomials  input(l) --> om(theta)  (distrib. in m)
      !-- with dealiasing of input:
      do j=1, blk_ps_size(mpi_rank)
         aux(0:Nt,j) = matmul(legendre(0:Nt,mm(j):Nt_s,j),input(mm(j):Nt_s,j))
      enddo
      !-- result is now in omegaarr(itheta,m)  (distr. in m)
      !-- transposition omegaarr distr(m) --> tt distr(theta):
      call m2phi_2D(aux,output)
      !-- result is now in tt_t(itheta,j)  (distr. in theta).
   end subroutine ylmb

   ! ----------------------------------------------------------------------
   !!--
   !!--    transformation   t(theta,m) --> t_t(theta,phi)
   !!-- and transposition  (distrib. m --> distrib. theta)
   !!--
   !!--  input:
   !!--         t(.,j) means m=m0*(j/2):       distrib. in m)
   !!--         real(m=0),real(m=m0),im(m=m0),real(m=2m0),im(m=2*m0),...
   !!--
   !!-- output:
   !!--         t_t(0:l(maxb-1),Np) transposed
   !!--         t_t(.,j) means position phi=(2pi/m0)*(j-1)/Np.
   !!--
   !-- 05/96 M.A.: improved order of indices: t(0:drs_Nt,drs_Np)
   !-- 06/96 M.A.: workarrays are local now
   !-- 09/96 M.A.: opt.: replaced big tranf. (m0*Np) by small (Np).
   !-- 10/96 M.A.: m2phi_2D is the parallel version of m2phi.
   !-- 20.10.96 M.A.: output array is transposed now.
   subroutine m2phi_2D(input,output)
      implicit none
      !-- input field to transform:
      double precision, intent(in):: input(0:Nt,blk_ps_size(mpi_rank))
      !-- output field (transposed):
      double precision, intent(out):: output(0:blk_t_size(mpi_rank)-1,Np)
      double precision:: aux(0:blk_t_size(mpi_rank)-1,Np)

      !-- transposition input distr(m) --> output distr(theta):
      call transpos_phi2theta(input, Nt, aux, Np)

      !-- result is now in output(theta,m)  (distr. in theta)
      !-- the FFT  output(m) --> output(phi):
      ! Dealiase in m
      if (Np_s+1 .le. Np) aux(:, Np_s+1:Np) = 0.d0

      call dft_backward(aux, output)
      !-- result is now in output(theta,phi)  (distr. in theta).
   end subroutine m2phi_2D

   ! -----------------------------------------------------------------------------
   !> Computes the spectral coefficients of the divergence of the vector field \f$\vec u\f$.\n
   !! On input:
   !! @param vec_r is \f$ u_r*r^2 \f$ in (l,m,r) space.
   !! @param vec_t is \f$ \frac{u_\theta}{r\sin\theta} \f$ in (l,m,r) space.
   !! @param vec_p is \f$ \frac{u_\phi}{r\sin\theta} \f$ in (l,m,r) space.
   !! On output:
   !! @param nonlin is the divergence in (l,m,r) space.
   subroutine my_div(vec_r,vec_t,vec_p,nonlin)
      implicit none

      double precision, intent(in)::  vec_r(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(in)::  vec_t(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(in)::  vec_p(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(out):: nonlin(0:Nt_s,blk_ps_size(mpi_rank),Nr)

      double precision:: deriv(Nr), factor_lp1, factor_lm1
      integer:: i, j, l, m, jg

      nonlin = 0.0d0

      do j=1, blk_ps_size(mpi_rank)
         m  = blk_ts_start(j)
         jg = blk_ps_start(mpi_rank) + j - 1
         do l=m, Nt_s
            deriv = radial_derivative_r2r(vec_r(l,j,1:Nr))
            ! factor_lp1 = l*(l+m+1)/dble(2*l+3)
            factor_lp1 = l*dsqrt( dble((l+m+1)*(l-m+1))/dble((2*l+1)*(2*l+3)) )
            ! factor_lm1 = (l+1)*(l-m)/dble(2*l-1)
            factor_lm1 = (l+1)*dsqrt( dble((l+m)*(l-m))/dble((2*l+1)*(2*l-1)) )
            r:do i=1, Nr
               nonlin(l,j,i) = 1.0d0/rcoll2(i)*deriv(i)
               if (l.eq.Nt_s) then
                  nonlin(l,j,i) = nonlin(l,j,i) + vec_t(l-1,j,i)*factor_lm1
               else
                  if (l.ne.0) then
                     nonlin(l,j,i) = nonlin(l,j,i) + vec_t(l-1,j,i)*factor_lm1 - &
                           vec_t(l+1,j,i)*factor_lp1
                  endif
               endif
               if (jg.gt.1) then
                  if (mod(jg,2).eq.0) then
                     if(j+1.le.blk_ps_size(mpi_rank)) then
                        nonlin(l,j,i) = nonlin(l,j,i) - m*vec_p(l,j+1,i)
                     endif
                  else
                     nonlin(l,j,i) = nonlin(l,j,i) + m*vec_p(l,j-1,i)
                  endif
               endif
            enddo r
         enddo
      enddo
   end subroutine my_div

   !*******************************************************************
   ! imported from my_2rot.F
   ! Revision 1.2  96/05/10  13:56:37  13:56:37  btpa08 (M. Ardes)
   !
   !*******************************************************************
   !--  input: vec_r,vec_t,vec_p  (lmr)
   !-- output: nonlin = rot_r(input)       (lmr)
   !--
   subroutine my_rot(vec_t,vec_p,nonlin)
      implicit none
      double precision, intent(in)::  vec_t(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(in)::  vec_p(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(out):: nonlin(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision:: factor_lp1, factor_lm1

      integer:: i,j,l,m
      integer:: jg

      nonlin = 0.0d0

      do i=1,Nr
         do j=1,blk_ps_size(mpi_rank)
            m  = blk_ts_start(j)
            jg = blk_ps_start(mpi_rank) + j - 1
            do l=m, Nt_s
               ! factor_lp1 = l*(l+m+1)/dble(2*l+3)
               factor_lp1 = l*dsqrt( dble((l+m+1)*(l-m+1))/dble((2*l+1)*(2*l+3)) )
               ! factor_lm1 = (l+1)*(l-m)/dble(2*l-1)
               factor_lm1 = (l+1)*dsqrt( dble((l+m)*(l-m))/dble((2*l+1)*(2*l-1)) )

               if (l.eq.Nt_s) then
                  nonlin(l,j,i) = vec_p(l-1,j,i)*factor_lm1
               elseif (l.eq.0) then
                  nonlin(l,j,i) = 0.d0
               else
                  nonlin(l,j,i) = vec_p(l-1,j,i)*factor_lm1-&
                                  vec_p(l+1,j,i)*factor_lp1
               endif
               if (jg.eq.1) cycle
               if (mod(jg,2).eq.0) then
                  if(j+1.le.blk_ps_size(mpi_rank)) then
                     nonlin(l,j,i) = nonlin(l,j,i) + m*vec_t(l,j+1,i)
                  endif
               else
                  nonlin(l,j,i) = nonlin(l,j,i) - m*vec_t(l,j-1,i)
               endif
            enddo
         enddo
      enddo
   end subroutine  my_rot

   ! ---------------------------------------------------------------
   !--  input:  rot_r,rot_t,rot_p    (lmr)
   !-- output:  nonlin = rotrot_r(input)      (lmr)
   !--
   subroutine my_rotrot(vec_r,vec_t,vec_p,nonlin)
      implicit none

      double precision, intent(in)::  vec_r(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(in)::  vec_t(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(in)::  vec_p(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(out):: nonlin(0:Nt_s,blk_ps_size(mpi_rank),Nr)

      double precision:: radarr(Nr),deriv(Nr)
      integer:: i, j, l, m, jg
      double precision:: factor_lm1, factor_lp1

      nonlin = 0.0d0

      do j=1, blk_ps_size(mpi_rank)
         m  = blk_ts_start(j)
         jg = blk_ps_start(mpi_rank) + j - 1
         do l=m,Nt_s
            ! factor_lp1 = l*(l+m+1)/dble(2*l+3)
            factor_lp1 = l*dsqrt( dble((l+m+1)*(l-m+1))/dble((2*l+1)*(2*l+3)) )
            ! factor_lm1 = (l+1)*(l-m)/dble(2*l-1)
            factor_lm1 = (l+1)*dsqrt( dble((l+m)*(l-m))/dble((2*l+1)*(2*l-1)) )

            if (l.eq.Nt_s) then
               radarr(1:Nr) = vec_t(l-1,j,1:Nr)*factor_lm1
            elseif (l.eq.0) then
               radarr = 0.d0
            else
               radarr(1:Nr) = vec_t(l-1,j,1:Nr)*factor_lm1 - &
                              vec_t(l+1,j,1:Nr)*factor_lp1
            endif
            if (jg.gt.1) then
               if (mod(jg,2).eq.0) then
                  if(j+1.le.blk_ps_size(mpi_rank)) then
                     radarr(1:Nr) = radarr(1:Nr) - m*vec_p(l,j+1,1:Nr)
                  endif
               else
                  radarr(1:Nr) = radarr(1:Nr) + m*vec_p(l,j-1,1:Nr)
               endif
            endif
            radarr(1:Nr) = radarr(1:Nr)*rcoll2(1:Nr)
            deriv = radial_derivative_r2r(radarr)
            do i = 1, Nr
               nonlin(l,j,i) = (llp1(l)/rcoll2(i)*vec_r(l,j,i) + deriv(i))/rcoll2(i)
            enddo
         enddo
      enddo
   end subroutine my_rotrot

   !> Computes the poloidal and toroidal scalar coefficients of a 3D vector field.\n
   !! This is the part of the calculation that is common to both the flow and
   !! magnetic field.
   subroutine vectorField2PolTor_common(vr,vt,vp,pol,tor)
      implicit none
      double precision,intent(inout):: pol(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision,intent(inout):: tor(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision,intent(inout):: vr(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision,intent(inout):: vt(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision,intent(inout):: vp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision:: ur(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      double precision:: ut(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      double precision:: up(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      integer:: t_end
      integer:: i, l, j
      t_end = blk_t_size(mpi_rank)-1

      forall (l=0:t_end, j=1:Np, i=1:Nr)
         vr(l,j,i) = vr(l,j,i)*rcoll2(i)
      endforall
      call ylmt_3D(vr(0:t_end,1:Np,1:Nr), ur(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr))

      forall (l=0:t_end, j=1:Np, i=1:Nr)
         vt(l,j,i) = vt(l,j,i)/(sintheta(blk_t_start(mpi_rank)+l)*rcoll(i))
      endforall
      call ylmt_3D(vt(0:t_end,1:Np,1:Nr), ut(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr))

      forall (l=0:t_end, j=1:Np, i=1:Nr)
         vp(l,j,i) = vp(l,j,i)/(sintheta(blk_t_start(mpi_rank)+l)*rcoll(i))
      endforall
      call ylmt_3D(vp(0:t_end,1:Np,1:Nr), up(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr))

      !-- my_rotrot computes tor(,,) = rotrot_r(u)    (toroidal)
      call my_rotrot(ur,ut,up,tor)
      !-- my_rot computes    pol(,,) = rot_r(u)       (poloidal)
      call my_rot(ut,up,pol)
   end subroutine vectorField2PolTor_common

   !> Multiplies the poloidal and toroidal scalar coefficients by the apropriate factors
   !! of r and l*(l+1) to make them the poloidal and toroidal scalars of the flow
   !! as defined for this program.
   subroutine PolTor_common2PolTor_flow(pol,tor)
      implicit none
      double precision,intent(inout):: pol(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision,intent(inout):: tor(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      integer:: i, l, j
      do i=1, Nr
         jl_do(j,l)
            if (l.ne.0) then
               tor(l,j,i) = -tor(l,j,i)/llp1(l)
            else
               tor(l,j,i) = 0.0d0
            endif
         jl_enddo
      enddo

      do i=1, Nr
         jl_do(j,l)
            if (l.ne.0) then
               pol(l,j,i) = rcoll(i)*pol(l,j,i)/llp1(l)
            else
               pol(l,j,i) = 0.0d0
            endif
         jl_enddo
      enddo
   end subroutine PolTor_common2PolTor_flow

   !> Multiplies the poloidal and toroidal scalar coefficients by the apropriate factors
   !! of r and l*(l+1) to make them the poloidal and toroidal scalars of the
   !! magnetic field as defined for this program.
   subroutine PolTor_common2PolTor_field(pol,tor)
      implicit none
      double precision,intent(inout):: pol(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision,intent(inout):: tor(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      integer:: i, l, j
      do i=1, Nr
         jl_do(j,l)
            if (l.ne.0) then
               tor(l,j,i) = -rcoll2(i)*tor(l,j,i)/llp1(l)
            else
               tor(l,j,i) = 0.0d0
            endif
         jl_enddo
      enddo

      do i=1, Nr
         jl_do(j,l)
            if (l.ne.0) then
               pol(l,j,i) = -rcoll2(i)*pol(l,j,i)/llp1(l)
            else
               pol(l,j,i) = 0.0d0
            endif
         jl_enddo
      enddo
   end subroutine PolTor_common2PolTor_field

   !> Computes the spherical harmonic coefficients of the divergence of a
   !! 3D vector field in real space.
   subroutine vectorField2Divergence(ur, ut, up, div)
      implicit none
      double precision, intent(inout):: ur(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(inout):: ut(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(inout):: up(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: div(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision:: ur_aux(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      double precision:: ut_aux(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      double precision:: up_aux(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      integer:: i,l,j, t_end
      t_end    = blk_t_size(mpi_rank)-1

      ! Transform ur
      forall (i=1:Nr, l=0:t_end, j=1:Np)
         ur(l,j,i) = ur(l,j,i)*rcoll2(i)
      endforall
      call ylmt_3D(ur(0:t_end,1:Np,1:Nr), ur_aux(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr))
      ! Transform ut
      forall (i=1:Nr, l=0:t_end, j=1:Np)
         ut(l,j,i) = ut(l,j,i)/(sintheta(blk_t_start(mpi_rank)+l)*rcoll(i))
      endforall
      call ylmt_3D(ut(0:t_end,1:Np,1:Nr), ut_aux(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr))
      ! Transform up
      forall (i=1:Nr, l=0:t_end, j=1:Np)
         up(l,j,i) = up(l,j,i)/(sintheta(blk_t_start(mpi_rank)+l)*rcoll(i))
      endforall
      call ylmt_3D(up(0:t_end,1:Np,1:Nr), up_aux(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr))
      ! Compute the spectral part of the divergence
      call my_div(ur_aux,ut_aux,up_aux,div)
   end subroutine


end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
