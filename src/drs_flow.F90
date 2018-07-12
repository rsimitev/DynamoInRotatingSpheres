! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
module drs_flow
#include "drsDefs.F90"
   use drs_mpi
   use drs_params
   use drs_dims
   use drs_legendre
   use drs_transforms
   use drs_radial
   implicit none
   save

   double precision, allocatable::     flow_pol(:,:,:)
   double precision, allocatable::     flow_tor(:,:,:)
   double precision, allocatable::  flow_pol_dr(:,:,:)
   double precision, allocatable::  flow_tor_dr(:,:,:)
   double precision, allocatable:: flow_pol_ddr(:,:,:)
   double precision, allocatable:: flow_tor_ddr(:,:,:)
   double precision, allocatable:: flow_pol_lap(:,:,:)
   double precision, allocatable:: flow_tor_lap(:,:,:)
   double precision, allocatable:: flow_pol_avg(:,:,:)
   double precision, allocatable:: flow_tor_avg(:,:,:)
   double precision, allocatable:: flow_r_t(:,:,:)
   double precision, allocatable:: flow_t_t(:,:,:)
   double precision, allocatable:: flow_p_t(:,:,:)
   double precision, allocatable:: rot_flow_r_t(:,:,:)
   double precision, allocatable:: rot_flow_t_t(:,:,:)
   double precision, allocatable:: rot_flow_p_t(:,:,:)
   logical, private:: initialised = .false.
   logical, private:: flow_allocated = .false.

contains
   !---------------------------------------------------------------------------
   !> Allocates all the variables needed to deal with a varying flow.
   subroutine drs_flow_allocation()
      implicit none
      if(flow_allocated) return
      allocate(    flow_pol(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(    flow_tor(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate( flow_pol_dr(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate( flow_tor_dr(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(flow_pol_ddr(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(flow_tor_ddr(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(flow_pol_lap(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(flow_tor_lap(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(flow_pol_avg(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(flow_tor_avg(0:Nt_s,blk_ps_size(mpi_rank),Nr))
      allocate(    flow_r_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(    flow_t_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(    flow_p_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(rot_flow_r_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(rot_flow_t_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(rot_flow_p_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      flow_allocated=.true.
   end subroutine

   !---------------------------------------------------------------------------
   !> Initialises the flow derivatives.
   !! Makes sure everyone knows the boundary conditions.
   subroutine drs_flow_init(error)
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
      if(.not.flow_allocated) then
         call drs_flow_allocation()
         flow_pol = 0.0d0
         flow_tor = 0.0d0
         error=-2 ! The flow was allocated and initialised here.
      endif

      call drs_bcast(flowBC, 2)
      call radial_dr_ddr_3D_r2r(flow_pol,flow_pol_dr,flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor,flow_tor_dr,flow_tor_ddr)
      call update_flow_pol_lap()
      call update_flow_tor_lap()
      initialised = .true.
   end subroutine

   !---------------------------------------------------------------------------
   !> Caches the toroidal part of the laplacian of the flow in lmr space
   subroutine update_flow_tor_lap()
      implicit none
      integer::l,j,i
      flow_tor_lap = 0.0d0
      do i=1, Nr
         jl_do(j,l)
            flow_tor_lap(l,j,i) = flow_tor_ddr(l, j, i) + &
                                  dble(2-llp1(l))/rcoll2(i)*flow_tor(l, j, i) + &
                                  4*flow_tor_dr(l, j, i)/rcoll(i)
         jl_enddo
      enddo
   end subroutine

   !---------------------------------------------------------------------------
   !> Caches the poloidal part of the laplacian of the flow in lmr space
   subroutine update_flow_pol_lap()
      implicit none
      integer::l,j,i
      flow_pol_lap = 0.0d0
      do i=1, Nr
         jl_do(j,l)
            flow_pol_lap(l,j,i) = flow_pol_ddr(l, j, i) + &
                                  2*flow_pol_dr(l, j, i)/rcoll(i)  - &
                                  dble(llp1(l))/rcoll2(i)*flow_pol(l, j, i)
         jl_enddo
      enddo
   end subroutine

   !---------------------------------------------------------------------------
   !>  Implement the boundary conditions on the RHS of the poloidal equation
   subroutine apply_flow_pol_BC(pol)
      implicit none
      double precision, intent(inout):: pol(Nr)
      ! Inner boundary
      pol(Nr) = 0.0d0
      ! Outer boundary
      pol(1)  = 0.0d0
   end subroutine

   !---------------------------------------------------------------------------
   !>  Implement the boundary conditions on the RHS of the toroidal equation
   subroutine apply_flow_tor_BC(tor)
      implicit none
      double precision, intent(inout):: tor(Nr)
      ! Inner boundary
      tor(Nr) = 0.0d0
      ! Outer boundary
      tor(1)  = 0.0d0
   end subroutine

   !---------------------------------------------------------------------------
   !> Abstracts computing the flow and its curl in real space.
   subroutine calc_u(ur, ut, up, rotu_r, rotu_t, rotu_p)
      implicit none

      double precision, intent(out):: ur(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: ut(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: up(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: rotu_r(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: rotu_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: rotu_p(0:blk_t_size(mpi_rank)-1,Np,Nr)

      call calc_flow(ur, ut, up)
      call calc_rot_flow(rotu_r, rotu_t, rotu_p)
   end subroutine calc_u

   !---------------------------------------------------------------------------
   !> This routine computes the spherical components of the flow (ur, ut, up)
   !! in real space:
   !!
   !! \f[
   !!  \vec{u} = \vec{\nabla}\times(\vec{\nabla}\times(\vec{r} flow\_pol))
   !!  + \vec{\nabla}\times(r \vec{r} flow\_tor)
   !! \f]
   subroutine calc_flow(ur_t, ut_t, up_t)
      implicit none
      double precision, intent(out):: ur_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: ut_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: up_t(0:blk_t_size(mpi_rank)-1,Np,Nr)

      double precision:: workarr(0:Nt_s,1:blk_ps_size(mpi_rank))
      double precision:: ulm(0:Nt,1:blk_ps_size(mpi_rank))

      integer:: i,j,l,jg
      integer:: blk_size, blk_start

      blk_start = blk_ps_start(mpi_rank)
      blk_size  = blk_ps_size(mpi_rank)

      ! call calc_component_r(flow_pol, ur_t)
      do i=1, Nr
         !-- get ur (which is pure poloidal)
         workarr = 0.0d0
         jl_do(j,l)
            workarr(l,j) = flow_pol(l,j,i)*llp1(l)/rcoll(i)
         jl_enddo
         !-- workarr(l,j) = ur in (lm)-space
         do j=1, blk_size
            !-- transf. ur(l,mm(j)) --> ur(theta,mm(j))
            ulm(0:Nt,j) = matmul(legendre(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo
         !-- transpos.  distrib(mm(j))--> distrib(theta)  and transform. mm(j)-->phi:
         call m2phi_2D(ulm, ur_t(:,:,i))
      enddo

      do i=1, Nr
         !-- get u_theta
         !-- first poloidal part:
         workarr = 0.0d0
         jl_do(j,l)
            workarr(l,j) = flow_pol_dr(l,j,i) + flow_pol(l,j,i)/rcoll(i)
         jl_enddo

         do j=1, blk_size
            !-- utheta(theta,mm(j),r) (pol.)
            ulm(0:Nt,j) = matmul(dleg(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo
            !-- then toroidal part, here we must distinguish between real and imag. part:

         workarr = 0.0d0
         do j=1, blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j - 1
               if(mod(jg,2).eq.0) then
                  !-- get real parts
                  if(j+1.le.blk_size) then
                     do l=mm(j), Nt_s
                        workarr(l,j) = -mm(j)*rcoll(i)*flow_tor(l,j+1,i)
                     enddo
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else
                  !-- get the imaginary parts
                  do l=mm(j), Nt_s
                     workarr(l,j) = mm(j)*rcoll(i)*flow_tor(l,j-1,i)
                  enddo
              endif
            endif  !  mm(j)>0
         enddo  !  j=1,blk_size

         do j=1, blk_size
            if(mm(j).gt.0) then         ! only mm(j)>0, see above
               !--  add utheta(theta,mm(j),r) (tor., mm(j)>0)
               ulm(0:Nt,j) = ulm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo
         call m2phi_2D(ulm, ut_t(:,:,i))
      enddo


      do i=1, Nr
         !-- get u_phi
         !-- first toroidal part:
         workarr = 0.0d0
         jl_do(j,l)
            workarr(l,j) = -rcoll(i)*flow_tor(l,j,i)
         jl_enddo

         do j=1,blk_size
            ulm(0:Nt,j) = matmul(dleg(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo

         !-- then poloidal part, here we must distinguish between real and imag. part:
         !CC        do mm(j)=m0,mmax,m0  ! only need mm(j)>0 because there's no pol. part mm(j)=0!
         do j=1, blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j - 1
               if(mod(jg,2).eq.0) then !-- get real parts
                  if(j+1 .le. blk_size) then
                     workarr(mm(j):Nt_s,j) = -mm(j)*( flow_pol_dr(mm(j):Nt_s,j+1,i) + flow_pol(mm(j):Nt_s,j+1,i)/rcoll(i) )
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else !-- get the imaginary parts
                  workarr(mm(j):Nt_s,j) = mm(j)*( flow_pol_dr(mm(j):Nt_s,j-1,i) + flow_pol(mm(j):Nt_s,j-1,i)/rcoll(i) )
               endif
            endif  !  mm(j)>0
         enddo  !  j=1,blk_size

         do j=1, blk_size
            if(mm(j).gt.0) then    ! only mm(j)>0, see above
               ulm(0:Nt,j) = ulm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo
         call m2phi_2D(ulm, up_t(:,:,i))
      enddo
   end subroutine calc_flow

   !---------------------------------------------------------------------------
   !> This routine computes the spherical components of the curl of the flow
   !! in spherical coordinates:
   !!
   !! \f[
   !!  \nabla\times\vec{u} = rot(rotrot(rP) + rot(rT))
   !! \f]
   subroutine calc_rot_flow(rotu_r, rotu_t, rotu_p)
      implicit none
      double precision, intent(out):: rotu_r(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: rotu_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(out):: rotu_p(0:blk_t_size(mpi_rank)-1,Np,Nr)

      double precision:: workarr(0:Nt_s,1:blk_ps_size(mpi_rank))
      double precision:: ulm(0:Nt,1:blk_ps_size(mpi_rank))

      integer:: i,j,l,jg
      integer:: blk_size, blk_start

      blk_start = blk_ps_start(mpi_rank)
      blk_size  = blk_ps_size(mpi_rank)

      rotu_r = 0.0d0
      rotu_t = 0.0d0
      rotu_p = 0.0d0
      do i=1, Nr  !  only one big r-loop for calc_u now.
         !-- get rotu_r, which is pure toroidal:
         jl_do(j,l)
            workarr(l,j) = llp1(l)*flow_tor(l,j,i)
         jl_enddo
         do j=1, blk_size
            ulm(0:Nt,j) = matmul(legendre(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo
         call m2phi_2D(ulm, rotu_r(:,:,i))
      enddo

      do i=1, Nr
         !-- get rotu_theta:
         !-- first the toroidal part:
         jl_do(j,l)
            workarr(l,j) = rcoll(i)*flow_tor_dr(l,j,i) + 2*flow_tor(l,j,i)
         jl_enddo
         do j=1, blk_size     ! all mm(j)'s!
            ulm(0:Nt,j) = matmul(dleg(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo

         !-- then the poloidal part:
         do j=1, blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j - 1
               if(mod(jg,2).eq.0) then
                  !-- get real parts
                  if(j+1.le.blk_size) then
                     forall(l=mm(j):Nt_s) workarr(l,j) = mm(j)*flow_pol_lap(l,j+1,i)
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else
                  !-- get the imaginary parts
                  forall(l=mm(j):Nt_s) workarr(l,j) = -mm(j)*flow_pol_lap(l,j-1,i)
               endif
            endif  !  mm(j)>0
         enddo  !  j=1,blk_size

         do j=1, blk_size
            if(mm(j).gt.0) then    !   only mm(j)>0, see above.
               ulm(0:Nt,j) = ulm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo  !  j=1,blk_size

         call m2phi_2D(ulm, rotu_t(:,:,i))
      enddo

      do i=1, Nr
         !-- get rotu_phi
         ulm = 0.0d0
         !-- first poloidal part:
         jl_do(j,l)
            workarr(l,j) = flow_pol_lap(l,j,i)
         jl_enddo

         do j=1, blk_size       !  all mm(j)'s.
            ulm(0:Nt,j) = matmul(dleg(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
         enddo

         !-- then toroidal part:
         do j=1, blk_size
            if(mm(j).gt.0) then
               jg = blk_ps_start(mpi_rank) + j - 1
               if(mod(jg,2).eq.0) then !-- get real parts
                  if(j+1.le.blk_size) then
                     workarr(mm(j):Nt_s,j) = -mm(j)*( rcoll(i)*flow_tor_dr(mm(j):Nt_s,j+1,i) + 2*flow_tor(mm(j):Nt_s,j+1,i) )
                  else
                     workarr(mm(j):Nt_s,j) = 0.0d0
                  endif
               else !-- get the imaginary parts
                  workarr(mm(j):Nt_s,j)= mm(j)*(rcoll(i)*flow_tor_dr(mm(j):Nt_s,j-1,i) + 2*flow_tor(mm(j):Nt_s,j-1,i))
               endif
            endif  !  mm(j)>0
         enddo  !  j=1,blk_size

         do j=1, blk_size
            if(mm(j).gt.0) then  !  only mm(j)>0.
               ulm(0:Nt,j) = ulm(0:Nt,j) + matmul(leg_sin(0:Nt,mm(j):Nt_s,j),workarr(mm(j):Nt_s,j))
            endif
         enddo
         call m2phi_2D(ulm, rotu_p(:,:,i))
      enddo  !  i=1,Nr
   end subroutine calc_rot_flow

   !---------------------------------------------------------------------------
   !> Computes the l-spectrum of the radial flow.
   subroutine calc_flow_lspec(uspec)
      implicit none
      double precision, intent(out):: uspec(0:Nt_s)
      integer:: i,j,jg,l,m,mmax
      double precision:: plm
      integer:: blk_size,blk_start

      blk_size = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)
      mmax = m0*(Np_s/2)
      uspec = 0.0
      do l=0,Nt_s
         do m=0,min(l,mmax),m0
            plm = plmfac(l,m)
            do i=1,Nr-1
               jg=2*(m/m0)+1
               j=jg-blk_start+1
               if(j.ge.1 .and. j.le.blk_size) then
                  uspec(l)=uspec(l)+((llp1(l)/rcoll(i)*plm*flow_pol(l,j,i))**2+&
                     (llp1(l)/rcoll(i+1)*plm*flow_pol(l,j,i+1))**2)*0.5*&
                     drcoll(i)
               endif
               if (m.gt.0) then
                  j=j-1
                  if(j.ge.1 .and. j.le.blk_size) then
                     uspec(l)=uspec(l)+((llp1(l)/rcoll(i)*plm*flow_pol(l,j,i))**2+&
                        (llp1(l)/rcoll(i+1)*plm*flow_pol(l,j,i+1))**2)*0.5*&
                        drcoll(i)
                  endif
               endif
            enddo
         enddo
      enddo
      call sum_over_all_cpus(uspec)
   end subroutine

   !---------------------------------------------------------------------------
   !> Computes the m-spectrum of the radial flow.
   subroutine calc_flow_mspec(uspec)
      implicit none
      double precision, intent(out):: uspec(m0*Np_s+1)
      integer:: i,j,jg,l,m,mmax
      double precision:: plm
      integer:: blk_size,blk_start

      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)

      mmax  = m0*(Np_s/2)
      uspec = 0.0
      do m=0, mmax, m0
         jg = 2*(m/m0)+1
         j  = jg - blk_start + 1
         if(j.gt.0 .and. j.le.blk_size) then
            do l=m,Nt_s
               plm = plmfac(l,m)
               do i=1,Nr-1
                  !-- imag. part
                  if(j.le.blk_size) then
                     uspec(jg) = uspec(jg) + ( &
                        ( llp1(l)/rcoll(i)*  plm*flow_pol(l,j,i)  )**2 + &
                        ( llp1(l)/rcoll(i+1)*plm*flow_pol(l,j,i+1))**2   &
                        )*0.5* drcoll(i)
                  endif      !  (j.le.blk_size)
                  if (m.gt.0) then
                     !-- add the real part
                     uspec(jg) = uspec(jg) + ( &
                        ( llp1(l)/rcoll(i)*  plm*flow_pol(l,j,i)   )**2 + &
                        ( llp1(l)/rcoll(i+1)*plm*flow_pol(l,j,i+1) )**2   &
                        )*0.5* drcoll(i)
                  endif
               enddo         !  i=1,Nr-1
            enddo            !  l=m,Nt_s
         endif               !  jj.gt.0
      enddo                  !  m=0,mmax
      call sum_over_all_cpus(uspec)
   end subroutine

   !---------------------------------------------------------------------------
   !> Computes the n-spectrum of the radial flow.
   subroutine calc_flow_nspec(uspec)
      implicit none
      double precision,intent(out):: uspec(Nr_s)
      integer:: i,j,jg,jj,l,m,mmax
      double precision:: plm
      double precision:: radarr(Nr)
      double precision:: aux(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      integer:: blk_size,blk_start

      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)
      mmax = m0*(Np_s/2)

      !-- power spectrum of Chebysheff coefficients
      uspec = 0.0

      do j=1,blk_size
         m = m0*((j+blk_start)/2)
         do l=m, Nt_s
            !c----- radial velocity
            radarr(1:Nr) = llp1(l)*flow_pol(l,j,1:Nr)/rcoll(1:Nr)
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
                     uspec(i) = uspec(i) + (plm*aux(l,j,i))**2
                  endif      !  (j.le.blk_size)
                  if (m.gt.0) then
                     !-- add the real part
                     j = j-1
                     uspec(i) = uspec(i) + (plm*aux(l,j,i))**2
                  endif
               enddo         !  l=m,Nt_s
            endif            !  jj.gt.0
         enddo               !  m=0,mmax
      enddo                  !  i=1,Nr
      call sum_over_all_cpus(uspec)
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
