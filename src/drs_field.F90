! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
module drs_field
#include "drsDefs.F90"
   use drs_mpi
   use drs_params
   use drs_dims
   use drs_legendre
   use drs_transforms
   use drs_radial
   use drs_hypDiff
   implicit none
   save

   double precision, allocatable::     field_pol(:,:,:)
   double precision, allocatable::     field_tor(:,:,:)
   double precision, allocatable::  field_pol_dr(:,:,:)
   double precision, allocatable::  field_tor_dr(:,:,:)
   double precision, allocatable:: field_pol_ddr(:,:,:)
   double precision, allocatable:: field_tor_ddr(:,:,:)
   double precision, allocatable:: field_pol_lap(:,:,:)
   double precision, allocatable:: field_tor_lap(:,:,:)
   double precision, allocatable:: field_pol_avg(:,:,:)
   double precision, allocatable:: field_tor_avg(:,:,:)
   double precision, allocatable:: field_r_t(:,:,:)
   double precision, allocatable:: field_t_t(:,:,:)
   double precision, allocatable:: field_p_t(:,:,:)
   double precision, allocatable:: rot_field_r_t(:,:,:)
   double precision, allocatable:: rot_field_t_t(:,:,:)
   double precision, allocatable:: rot_field_p_t(:,:,:)

   logical, private:: field_allocated=.false.
   logical, private:: initialised=.false.
contains

   !------------------------------------------------------------------
   !> Allocates arrays related to the field and its derivatives in
   !! spectral and real space.
   subroutine drs_field_allocation()
      implicit none
      if(field_allocated) return
      allocate(    field_pol(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(    field_tor(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate( field_pol_dr(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate( field_tor_dr(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(field_pol_ddr(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(field_tor_ddr(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(field_pol_lap(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(field_tor_lap(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(field_pol_avg(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(field_tor_avg(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(    field_r_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(    field_t_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(    field_p_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(rot_field_r_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(rot_field_t_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(rot_field_p_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      field_allocated=.true.
   end subroutine

   !---------------------------------------------------------------------------
   !> If the field module was not allocated, allocate and set it to zero.
   !! Compute the derivatives and broadcast the boundary conditions.
   subroutine drs_field_init(error)
      implicit none
      integer, intent(out):: error
      error=0
      ! If we were already initialised, do nothing
      if(initialised) then
         error=-1
         return
      endif
      ! If we have not been previously allocated, allocate and
      ! initialise to 0
      if(.not.field_allocated) then
         call drs_field_allocation()
         field_pol = 0.0d0
         field_tor = 0.0d0
         error=-2 ! The flow was allocated and initialised here.
      endif
      call drs_bcast(magBC, 2)
      call radial_dr_ddr_3D_r2r(field_tor,field_tor_dr,field_tor_ddr)
      call radial_dr_ddr_3D_r2r(field_pol,field_pol_dr,field_pol_ddr)
      call update_field_pol_lap()
      call update_field_tor_lap()
      initialised= .true.
   end subroutine

   !---------------------------------------------------------------------------
   !> Caches the toroidal part of the Laplacian of the field.
   subroutine update_field_tor_lap()
      implicit none
      integer::l,j,i
      do i=1, Nr
         jl_do(j,l)
            field_tor_lap(l, j, i) = field_tor_ddr(l, j, i) - llp1(l)*field_tor(l, j, i)/rcoll2(i)
            call drs_apply_hypDiff(field_tor_lap(l, j, i), l)
         jl_enddo
      enddo
   end subroutine

   !---------------------------------------------------------------------------
   !> Caches the poloidal part of the Laplacian of the field.
   subroutine update_field_pol_lap()
      implicit none
      integer::l,j,i
      do i=1, Nr
         jl_do(j,l)
            field_pol_lap(l, j, i) = field_pol_ddr(l, j, i) - llp1(l)*field_pol(l, j, i)/rcoll2(i)
            call drs_apply_hypDiff(field_pol_lap(l, j, i), l)
         jl_enddo
      enddo
   end subroutine

   !---------------------------------------------------------------------------
   !>
   subroutine drs_field_random_init(noise)
      implicit none
      double precision, intent(in):: noise
      double precision:: ranr
      integer:: ran, j, jg, m, l
      ran=1
      do j=1,min(6,Np_s)
         jg = blk_ps_start(mpi_rank) + j-1
         m  = m0*(jg/2)
         do l=max(1,m), min(4*m0,Nt_s)
            ran  = mod(ran*106+1283,6075)
            ranr = 1.0*ran/6075.0
            field_tor(l,j,1:Nr) = field_tor(l,j,1:Nr) + noise*ranr*sin(pi*(rcoll(1:Nr)-eta/(1-eta)))/(plmfac(l,m)*(l+1)**2*(m+1)**4)
         enddo
      enddo
   end subroutine drs_field_random_init

   !---------------------------------------------------------------------------
   !>  These lines take care of boundary conditions
   !!  If the value at a boundary is different from 0, insert bc/2
   !!  because of the factor 2 between Chebychev and radial differentiation
   !!  If the boundary condition is not on the derivative but the value of the
   !!  function itself (as is the case for the temperature), no factor of 2 is needed
   !!  if the boundary value bc is different from 0.
   subroutine apply_field_pol_BC(pol, l, m)
      implicit none
      double precision, intent(inout):: pol(Nr)
      integer, intent(in):: l, m
      pol(1)  = 0.0d0
      pol(Nr) = 0.0d0
   end subroutine

   !---------------------------------------------------------------------------
   !>  These lines take care of boundary conditions
   !!  If the value at a boundary is different from 0, insert bc/2
   !!  because of the factor 2 between Chebychev and radial differentiation
   !!  If the boundary condition is not on the derivative but the value of the
   !!  function itself (as is the case for the temperature), no factor of 2 is needed
   !!  if the boundary value bc is different from 0.
   subroutine apply_field_tor_BC(tor, l, m)
      implicit none
      double precision, intent(inout):: tor(Nr)
      integer, intent(in):: l, m
      tor(1)  = 0.0d0
      tor(Nr) = 0.0d0

      ! Imposed field at the boundaries for the case of magneto-convection
      ! TODO: This needs to be made configurable!
      if((drs_calc_type.eq.MagnetoConvection).and.(m.eq.0).and.(l.eq.2)) then
         tor(1)  = -38.49*rcoll(1)
         tor(Nr) = -38.49*rcoll(Nr)
      endif
   end subroutine

   !---------------------------------------------------------------------------
   !> This routine computes:
   !!~~~~~
   !!     B = rotrot(rH) + rot(rG)
   !!  rotB = rot(rotrot(rH) + rot(rG))
   !!
   !! the fields are defined as:
   !! field in program   field in equation
   !!-----------------   ---------------
   !!   field_pol     =   rH
   !!   field_tor     =   rG
   !!
   !! difference to calc_u():
   !! a) all equations are multiplied by 1/r
   !! b) all terms which include d/dr are different:
   !!    d/dr(rP) = p + r d/dr(p)  but
   !!    d/dr(rH) = d/dr(field_pol).
   !!
   !!  input: (not modified on output)
   !!     common /Bfields/ field_tor,field_pol    (lmr)
   !!     common /Bderivatives/ field_tor_dr,field_tor_ddr,field_pol_dr,field_pol_ddr     (lmr)
   !!
   !! output: (theta,phi, transposed):
   !!     Br_t,Bt_t,Bp_t,rot_Br_t,rot_Bt_t,rot_Bp_t
   !!
   !!    06.96 M.A.: unique variable names in common blocks
   !! 05.09.96 M.A.: optim.: removed "if(mm(j)>0)".
   !! 04.10.96 M.A. parallelized.
   !!~~~~~~
   subroutine calc_B(Br_t, Bt_t, Bp_t, rot_Br_t, rot_Bt_t, rot_Bp_t)
      implicit none

      double precision,intent(out):: Br_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision,intent(out):: Bt_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision,intent(out):: Bp_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision,intent(out):: rot_Br_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision,intent(out):: rot_Bt_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision,intent(out):: rot_Bp_t(0:blk_t_size(mpi_rank)-1,Np,Nr)

      call calc_field(Br_t, Bt_t, Bp_t)
      call calc_rot_field(rot_Br_t, rot_Bt_t, rot_Bp_t)
   end subroutine calc_B

   !---------------------------------------------------------------------------
   !> Computes the three components (Br, Bt, Bp) of the field in physical space.
   subroutine calc_field(Br, Bt, Bp)
      implicit none
      double precision, intent(out):: Br(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: Bt(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: Bp(0:blk_t_size(mpi_rank)-1,Np,Nr)

      double precision:: workarr(0:Nt_s,1:blk_ps_size(mpi_rank))
      double precision:: Blm(0:Nt,1:blk_ps_size(mpi_rank))

      integer:: i,j,jg,l
      integer:: blk_size, blk_start
      !--
      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)

      do i=1, Nr
         !-- get Br, which is pure poloidal
         jl_do(j,l)
            workarr(l,j) = field_pol(l,j,i)*llp1(l)/rcoll2(i)
         jl_enddo
         do j=1, blk_size
            !-- transf. Br(l,mm(j)) --> Br(theta,mm(j))
            Blm(0:Nt,j) = matmul(legendre(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo
         !-- transposition dist(mm(j))--> dist(theta) and FFT mm(j)-->phi:
         call m2phi_2D(Blm, Br(:,:,i))
      enddo

      do i=1, Nr
         !-- get Btheta
         !-- first poloidal part:
         jl_do(j,l)
            workarr(l,j) = field_pol_dr(l,j,i)/rcoll(i)
         jl_enddo
         do j=1, blk_size   !  all mm(j)'s.
            Blm(0:Nt,j) = matmul(dleg(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo
         !-- then toroidal part:
         !CC    do mm(j)=m0,mmax,m0  ! only mm(j)>0, but distinction between real and imag.
         do j=1, blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j-1
               if(mod(jg,2).eq.0) then
                  !-- get real parts
                  if(j+1.le.blk_size) then
                     workarr(mm(j):Nt_s,j) = -mm(j)*field_tor(mm(j):Nt_s,j+1,i)/rcoll(i)
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else
                  !-- get the imaginary parts
                  workarr(mm(j):Nt_s, j) = mm(j)*field_tor(mm(j):Nt_s,j-1,i)/rcoll(i)
               endif
            endif  !  mm(j)>0
         enddo
         do j=1,blk_size     !  only mm(j)>0.
            if(mm(j).gt.0) then
               Blm(0:Nt,j) = Blm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo
         call m2phi_2D(Blm, Bt(:,:,i))
      enddo

      do i=1, Nr
         !-- get Bphi
         !-- first toroidal part:
         do j=1,blk_size  !  all mm(j)'s.
            Blm(0:Nt,j) =  -matmul(dleg(0:Nt,mm(j):Nt_s,j),field_tor(mm(j):Nt_s,j,i))/rcoll(i)
         enddo

         !-- then poloidal part:
         !CC    do mm(j)=m0,mmax,m0  !  only mm(j)>0, but dist. between real and imag
         do j=1, blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j -1
               if(mod(jg,2).eq.0) then
               !-- get real parts
                  if(j+1.le.blk_size) then
                     workarr(mm(j):Nt_s,j)= -mm(j)*field_pol_dr(mm(j):Nt_s,j+1,i)/rcoll(i)
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else
                  !-- get the imaginary parts
                  workarr(mm(j):Nt_s,j) = mm(j)*field_pol_dr(mm(j):Nt_s,j-1,i)/rcoll(i)
               endif
            endif  !  mm(j)>0
         enddo

         do j=1,blk_size
            if(mm(j).gt.0) then  !  only mm(j)>0
               Blm(0:Nt,j) = Blm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo
         call m2phi_2D(Blm, Bp(:,:,i))
      enddo  !  i=1,Nr
   end  subroutine  calc_field

   !---------------------------------------------------------------------------
   !> Computes the three components rot(Br, Bt, Bp) of the curl of the field in
   !! physical space.
   subroutine calc_rot_field(rotB_r, rotB_t, rotB_p)
      implicit none
      double precision, intent(out):: rotB_r(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: rotB_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: rotB_p(0:blk_t_size(mpi_rank)-1,Np,Nr)

      double precision:: workarr(0:Nt_s,1:blk_ps_size(mpi_rank))
      double precision:: Blm(0:Nt,1:blk_ps_size(mpi_rank))

      integer:: i,j,jg,l
      integer:: blk_size, blk_start

      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)

      do i=1, Nr
         !-- get rotB_r, which is pure toroidal
         Blm = 0.0d0
         jl_do(j,l)
            workarr(l,j) = field_tor(l,j,i)*llp1(l)/rcoll2(i)
         jl_enddo
         do j=1,blk_size
            Blm(0:Nt,j) = matmul(legendre(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo
         call m2phi_2D(Blm,rotB_r(:,:,i))
      enddo

      do i=1, Nr
      !-- get rotB_theta
         Blm = 0.0d0
      !-- first toroidal part:
         jl_do(j,l)
            workarr(l,j) = field_tor_dr(l,j,i)/rcoll(i)
         jl_enddo
         do j=1,blk_size   !  all mm(j)'s.
            Blm(0:Nt,j) = matmul(dleg(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo

         !-- then poloidal part:
         !CC   do mm(j)=m0,mmax,m0  ! only mm(j)>0, but dist. between real and imag
         do j=1,blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j -1
               if(mod(jg,2).eq.0) then
                  !-- get real parts
                  if(j+1.le.blk_size) then
                     forall(l=mm(j):Nt_s) workarr(l,j) = field_pol_lap(l,j+1,i)*mm(j)/rcoll(i)
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else
                  !-- get the imaginary parts
                  forall(l=mm(j):Nt_s) workarr(l,j) = -field_pol_lap(l,j-1,i)*mm(j)/rcoll(i)
               endif
            endif  !  mm(j)>0
         enddo  !  j=1,blk_size
         do j=1, blk_size
            if(mm(j).gt.0) then  ! only mm(j)>0.
               Blm(0:Nt,j) = Blm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo
         call m2phi_2D(Blm, rotB_t(:,:,i))
      enddo

      do i=1, Nr
         !-- get rotB_phi
         Blm = 0.0d0
         !-- first poloidal part:
         jl_do(j,l)
            workarr(l,j) = field_pol_lap(l,j,i)/rcoll(i)
         jl_enddo
         do j=1, blk_size  !  all mm(j)'s.
            Blm(0:Nt,j) = matmul(dleg(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo

         !-- then toroidal part:
         !CC    do mm(j)=m0,mmax,m0  !  only mm(j)>0, but distinction between real and imag.
         do j=1,blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j -1
               if(mod(jg,2).eq.0) then
                  !-- get real parts
                  if(j+1.le.blk_size) then
                     workarr(mm(j):Nt_s,j) = -field_tor_dr(mm(j):Nt_s,j+1,i)*mm(j)/rcoll(i)
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else
                  !-- get the imaginary parts
                  workarr(mm(j):Nt_s,j) = field_tor_dr(mm(j):Nt_s,j-1,i)*mm(j)/rcoll(i)
               endif
            endif  !  mm(j)>0
         enddo  !  j=1,blk_size

         do j=1,blk_size   ! only mm(j)>0.
            if(mm(j).gt.0) then
               Blm(0:Nt,j) = Blm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo
         call m2phi_2D(Blm,rotB_p(:,:,i))
      enddo  !  i=1,Nr
   end subroutine  calc_rot_field

   !---------------------------------------------------------------------------
   !> Computes the l-spectrum of the magnetic field.\n
   !! @param Bspec is an array of Nt_s+1 elements starting from 0 corresponding
   !! to the power of the l.th SH degree of the field.
   subroutine calc_field_lspec(Bspec)
      implicit none
      double precision, intent(out):: bspec(0:Nt_s)
      integer:: i,j,jg,l,m,mmax
      double precision:: plm
      integer:: blk_size,blk_start

      blk_size = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)
      mmax = m0*(Np_s/2)

      bspec = 0.0
      do l=0,Nt_s
         do m=0,min(l,mmax),m0
            plm = plmfac(l,m)
            do i=1,Nr-1
               jg=2*(m/m0)+1
               j=jg-blk_start+1
               if(j.ge.1 .and. j.le.blk_size) then
                  bspec(l) = bspec(l) + ( &
                        ( llp1(l)/rcoll2(i)*  plm*field_pol(l,j,i)   )**2 + &
                        ( llp1(l)/rcoll2(i+1)*plm*field_pol(l,j,i+1) )**2   &
                        )*0.5d0*drcoll(i)
               endif
               if (m.gt.0) then
                  j = j - 1
                  if(j.ge.1 .and. j.le.blk_size) then
                     bspec(l) = bspec(l) + ( &
                        ( llp1(l)/rcoll2(i)*  plm*field_pol(l,j,i)   )**2 + &
                        ( llp1(l)/rcoll2(i+1)*plm*field_pol(l,j,i+1) )**2   &
                        )*0.5d0*drcoll(i)
                  endif
               endif
            enddo
         enddo
      enddo
      call sum_over_all_cpus(bspec)
   end subroutine

   !---------------------------------------------------------------------------
   !> Computes the m-spectrum of the magnetic field.\n
   !! @param Bspec is an array of m0*Np_s+1 elements starting from 1 corresponding
   !! to the power of the m.th SH order of the field.
   subroutine calc_field_mspec(Bspec)
      implicit none
      double precision, intent(out):: Bspec(m0*Np_s+1)
      integer:: i,j,jg,l,m,mmax
      double precision:: plm
      integer:: blk_size,blk_start

      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)

      mmax = m0*(Np_s/2)
      Bspec = 0.0
      do m=0, mmax, m0
         jg = 2*(m/m0)+1
         j  = jg - blk_start + 1
         if(j.gt.0 .and. j.le.blk_size) then
            do l=m,Nt_s
               plm = plmfac(l,m)
               do i=1,Nr-1
                  !-- imag. part
                  if(j.le.blk_size) then
                     Bspec(jg) = Bspec(jg) + ( &
                        ( llp1(l)/rcoll2(i)*  plm*field_pol(l,j,i)   )**2 +  &
                        ( llp1(l)/rcoll2(i+1)*plm*field_pol(l,j,i+1) )**2    &
                        )*0.5*drcoll(i)
                  endif      !  (j.le.blk_size)
                  if (m.gt.0) then
                     !-- add the real part
                     Bspec(jg) = Bspec(jg) + ( &
                        ( llp1(l)/rcoll2(i)*  plm*field_pol(l,j,i)   )**2 + &
                        ( llp1(l)/rcoll2(i+1)*plm*field_pol(l,j,i+1) )**2   &
                        )*0.5*drcoll(i)
                  endif
               enddo         !  i=1,Nr-1
            enddo            !  l=m,Nt_s
         endif               !  jj.gt.0
      enddo                  !  m=0,mmax
      call sum_over_all_cpus(bspec)
   end subroutine

   !---------------------------------------------------------------------------
   !> Computes the n-spectrum of the magnetic field.\n
   !! @param Bspec is an array of Nr_s elements starting from 1 corresponding
   !! to the power of the n.th Chebychev degree of the field.
   subroutine calc_field_nspec(Bspec)
      implicit none
      double precision, intent(out):: Bspec(Nr_s)
      integer:: i,j,jg,jj,l,m,mmax
      double precision:: plm
      double precision:: radarr(Nr)
      double precision:: aux(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      integer:: blk_size,blk_start

      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)
      mmax = m0*(Np_s/2)
      bspec(1:Nr) = 0.0
      do j=1,blk_size
         m = m0*((j+blk_start)/2)
         do l=m, Nt_s
            radarr(1:Nr) = llp1(l)*field_pol(l,j,1:Nr)/rcoll2(1:Nr)
            call cos_r2r_1_r2n(radarr)
            aux(l,j,1:Nr) = radarr(1:Nr)
         enddo
      enddo
      do i=1,Nr_s
         do m=0, mmax, m0
            jg = 2*(m/m0) + 1
            jj = jg - blk_start + 1
            if(jj.gt.0 .and. jj.le.blk_size) then
               do l=m, Nt_s
                  plm = plmfac(l,m)
                  !-- imag. part
                  j = jg - blk_start+1
                  if(j.le.blk_size) then
                     bspec(i) = bspec(i) + (plm*aux(l,j,i))**2
                  endif      !  (j.le.blk_size)
                  if (m.gt.0) then
                     !-- add the real part
                     j = j-1
                     bspec(i) = bspec(i) + (plm*aux(l,j,i))**2
                  endif
               enddo         !  l=m,Nt_s
            endif            !  jj.gt.0
         enddo               !  m=0,mmax
      enddo                  !  i=1,Nr
      call sum_over_all_cpus(bspec)
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
