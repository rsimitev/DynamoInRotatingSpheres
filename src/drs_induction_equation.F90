! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2015
!> Takes care of evolving the induction equation.
#include "drsDefs.F90"
module drs_induction_equation
   use drs_mpi
   use drs_flow
   use drs_field
   use drs_time
   use drs_dims
   use drs_params
   implicit none
   !> The Crank-Nicholson inverse operator for the toroidal flow.
   double precision, allocatable:: inv_lhs_IE_tor(:,:,:)
   !> The Crank-Nicholson inverse operator for the poloidal flow.
   double precision, allocatable:: inv_lhs_IE_pol(:,:,:)
   double precision, allocatable:: rhs_IE_tor_old(:,:,:),rhs_IE_tor(:,:,:)
   double precision, allocatable:: rhs_IE_pol_old(:,:,:),rhs_IE_pol(:,:,:)
   logical, private:: initialised=.false.

contains

   !---------------------------------------------------------------
   !> Allocates memory for the Crank-Nicholson inverse operators.
   !! Ensures all dependences are initialised.
   !! @param error Error code to be forwarded.
   subroutine drs_induction_equation_init(error)
      implicit none
      integer, intent(out):: error
      error=0
      if(initialised) then
         error=-1
         return
      endif
      ! Need to make sure the dimensions are already set.
      call mpi_dims_init(Nt, Np_s, m0, error)
      if(error.gt.0) return
      ! Need to make sure the field is initialised
      call drs_field_init(error)
      if(error.gt.0) return

      if(flow_present) then
         call drs_flow_init(error)
         if(error.gt.0) return
      endif

      allocate(rhs_IE_tor(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(rhs_IE_pol(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(rhs_IE_tor_old(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(rhs_IE_pol_old(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(inv_lhs_IE_tor(Nr,Nr,0:Nt_s))
      allocate(inv_lhs_IE_pol(Nr,Nr,0:Nt_s))
      rhs_IE_tor_old = 0.0d0
      rhs_IE_pol_old = 0.0d0
      initialised = .true.
      ! Generate the CN matrices for the field with factor Pm.
      call cnBp_init(inv_lhs_IE_pol, h, Pm, error)
      if(error.gt.0) return
      call cnBt_init(inv_lhs_IE_tor, h, Pm, error)
      if(error.gt.0) return
   end subroutine

   !---------------------------------------------------------------
   !> Computes matrices for the Crank-Nicholson scheme on the toroidal field.
   !! @par dinv is the resulting matrix.
   !! @par dt is the timestep currently being used
   !! @par fact is a factor for the time step, tipically, the magnetic Prantl number.
   !! @param error Error code to be forwarded.
   subroutine cnBt_init(dinv, dt, fac,error)
      implicit none
      double precision, target, intent(out):: dinv(Nr,Nr,0:Nt_s)
      integer, intent(out):: error
      double precision, intent(in):: dt, fac
      double precision, pointer:: aux(:,:)
      double precision:: dto2fac
      integer:: indx(Nr)
      integer:: i,j,l
      error=0

      if(.not.initialised) then
         error=1
         dinv=0.0d0
         return
      endif

      dto2fac = dt/(2*fac)
      do l=0, Nt_s
         aux => dinv(1:Nr,1:Nr,l)
         do j=1,Nr ! mode
            do i=1,Nr ! radius
               aux(i,j) = - poly_ddr(i,j) + llp1(l)*poly(i,j)/rcoll2(i)
               call drs_apply_hypDiff(aux(i,j), l)
               aux(i,j) = aux(i,j)*dto2fac + poly(i,j)
            enddo
         enddo

         ! Implement boundary conditions: fixed values at boundaries
         if ((magBC(1).eq.Vacuum).or.(magBC(1).eq.PseudoVacuum)) then
            aux(Nr,1:Nr) = poly(Nr,1:Nr)
         endif
         if ((magBC(2).eq.Vacuum).or.(magBC(2).eq.PseudoVacuum)) then
            aux(1, 1:Nr) = poly(1, 1:Nr)
         endif

         ! Dealiase
         if (Nr_s+1 .le. Nr) aux(:, Nr_s+1:Nr) = 0.d0

         call ludcmp(aux,Nr,Nr,indx)

         !-- invert all cn-matrixes for tor. field:
         call matinv(aux, indx, Nr, Nr)
      enddo
   end subroutine cnBt_init

   !---------------------------------------------------------------
   !> Computes matrices for the Crank-Nicholson scheme on the Poloidal field.
   !! @par dinv is the resulting matrix.
   !! @par dt is the timestep currently being used
   !! @par fact is a factor for the time step, tipically, the magnetic Prantl number.
   !! @param error Error code to be forwarded.
   subroutine cnBp_init(dinv, dt, fac, error)
      implicit none

      double precision, target, intent(out):: dinv(Nr,Nr,0:Nt_s)
      integer, intent(out):: error
      double precision, intent(in):: dt, fac
      double precision, pointer:: aux(:,:)
      double precision:: dto2fac
      integer:: indx(Nr)
      integer:: i,j,l
      error=0

      if(.not.initialised) then
         error=1
         dinv = 0.0d0
         return
      endif

      dto2fac = dt/(2.0d0*fac)
      do l=0, Nt_s
         aux => dinv(1:Nr,1:Nr,l)
         do j=1,Nr ! Mode
            do i=1,Nr ! Radius
               aux(i,j) = - poly_ddr(i,j) + llp1(l)*poly(i,j)/rcoll2(i)
               call drs_apply_hypDiff(aux(i,j), l)
               aux(i,j) = aux(i,j)*dto2fac + poly(i,j)
            enddo
         enddo

         if(magBC(1).eq.Vacuum) then
            ! implement boundary cond.: vacuum at inner and outer bound.
            ! inner bound.:  (d/dr - l/r) field_pol = 0
            ! outer bound.:  (d/dr + (l+1)/r) field_pol = 0
            ! note: i=1 corresponds to r_o, i=Nr corresponds to r_i
            aux(Nr,1:Nr_s) = -dble(l+1)/rcoll(Nr)
            do j=1, Nr_s
               aux(Nr,j) = (-1.0d0)**(j+1)*aux(Nr,j) + poly_dr(Nr,j)
            enddo
         elseif(magBC(1).eq.PseudoVacuum) then
            ! implement boundary cond.: pseudo-vacuum at inner and outer bound.
            ! inner bound.:  d/dr field_pol = 0
            ! outer bound.:  d/dr field_pol = 0
            ! note: i=1 corresponds to r_o, i=Nr corresponds to r_i
            do j=1, Nr_s
               aux(Nr,j) = poly_dr(Nr,j)
            enddo
         endif
         if(magBC(2).eq.Vacuum) then
            ! implement boundary cond.: vacuum at inner and outer bound.
            ! inner bound.:  (d/dr - l/r) field_pol = 0
            ! outer bound.:  (d/dr + (l+1)/r) field_pol = 0
            ! note: i=1 corresponds to r_o, i=Nr corresponds to r_i
            aux(1, 1:Nr_s) =  dble(l)/rcoll(1)
            do j=1, Nr_s
               aux(1, j) = aux(1, j) + poly_dr(1, j)
            enddo
         elseif(magBC(2).eq.PseudoVacuum) then
            ! implement boundary cond.: pseudo-vacuum at inner and outer bound.
            ! inner bound.:  d/dr field_pol = 0
            ! outer bound.:  d/dr field_pol = 0
            ! note: i=1 corresponds to r_o, i=Nr corresponds to r_i
            do j=1, Nr_s
               aux(1, j) = poly_dr(1, j)
            enddo
         endif

         ! Dealiase
         if (Nr_s+1 .le. Nr) aux(:, Nr_s+1:Nr) = 0.d0

         call ludcmp(aux,Nr,Nr,indx)

         !-- invert all cn-matrixes for pol. field:
         call matinv(aux, indx, Nr, Nr)
      enddo
   end subroutine

   !---------------------------------------------------------------
   !> Computes the non-linear terms of the induction equation.
   subroutine InductionEquation(rhs_IE_tor,rhs_IE_pol)
      implicit none
      double precision, intent(out):: rhs_IE_tor(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision, intent(out):: rhs_IE_pol(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision:: u_field_r(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision:: u_field_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision:: u_field_p(0:blk_t_size(mpi_rank)-1,Np,Nr)

      rhs_IE_tor = 0.0d0
      rhs_IE_pol = 0.0d0

      !-- induction term B x u:
      if(flow_present) then
         u_field_r = field_t_t*flow_p_t - field_p_t*flow_t_t
         u_field_t = field_p_t*flow_r_t - field_r_t*flow_p_t
         u_field_p = field_r_t*flow_t_t - field_t_t*flow_r_t
         !-- B now contains (B x u)    (theta,phi,r)

         ! Prepare to compute poloidal and toroidal terms
         call vectorField2PolTor_common(u_field_r,u_field_t, u_field_p, rhs_IE_pol,rhs_IE_tor)
         call PolTor_common2PolTor_field(rhs_IE_pol,rhs_IE_tor)
      endif
   end subroutine InductionEquation

   !---------------------------------------------------------------
   !> Applies boundary conditions and updates radial derivatives.
   subroutine just_apply_mag_BC()
      implicit none
      double precision:: field_tor_new(Nr)
      double precision:: field_pol_new(Nr)
      integer:: j, l
      jl_do(j,l)
         field_pol_new(1:Nr) = field_pol(l, j, 1:Nr)
         call apply_field_pol_BC(field_pol_new, l, mm(j))
         field_pol(l,j,1:Nr) = field_pol_new(1:Nr)
      jl_enddo
      call radial_dr_ddr_3D_r2r(field_pol, field_pol_dr, field_pol_ddr)

      ! Toroidal field
      jl_do(j,l)
         field_tor_new(1:Nr) = field_tor(l, j, 1:Nr)
         call apply_field_tor_BC(field_tor_new, l, mm(j))
         field_tor(l,j,1:Nr) = field_tor_new(1:Nr)
      jl_enddo
      call radial_dr_ddr_3D_r2r(field_tor, field_tor_dr, field_tor_ddr)
   end subroutine

   !---------------------------------------------------------------
   !> Updates the magnetic field using a hybrid Crank-Nicholson/Addams-Bashford implicit
   !! integration scheme.
   !! @param h is the time step;
   !! @param h1 and @param h2 are the Addams-Bashford weights;
   !! @param Pm is the magnetic Prandtl number;
   !! @param error Error code to be forwarded.
   subroutine update_field(h, h1, h2, Pm, error)
      implicit none
      double precision, intent(in):: h, h1, h2, Pm
      integer, intent(out)::error
      double precision:: ho2Pm
      double precision:: field_tor_new(Nr)
      double precision:: field_pol_new(Nr)
      integer:: l, j, i

      error=0
      if(.not.initialised) then
         error=1
         return
      endif

      if(.not.field_evolves) return

      ho2Pm = h/(2*Pm)
      ! The poloidal field
      jl_do(j,l)
         do i=1, Nr
            field_pol_new(i) = field_pol(l, j, i) + field_pol_lap(l, j, i)*ho2Pm &
                             + rhs_IE_pol(l, j, i)*h1 - rhs_IE_pol_old(l, j, i)*h2
         enddo
         call apply_field_pol_BC(field_pol_new, l, mm(j))
         ! multiply with inverted CN-matrices:
         ! This leaves all the fields in lmn space
         field_pol(l,j,1:Nr) = matmul(inv_lhs_IE_pol(1:Nr, 1:Nr, l),field_pol_new(1:Nr))
      jl_enddo
      rhs_IE_pol_old = rhs_IE_pol
      !-- all fields are given in (lmn) space
      !-- new derivatives are in (lmr).
      ! All fields are dealised in radial_dr_ddr_3D_n2r.
      call radial_dr_ddr_3D_n2r(field_pol, field_pol_dr, field_pol_ddr)

      ! The toroidal field
      jl_do(j,l)
         do i=1, Nr
            field_tor_new(i) = field_tor(l, j, i) + field_tor_lap(l, j, i)*ho2Pm &
                             + rhs_IE_tor(l, j, i)*h1 - rhs_IE_tor_old(l, j, i)*h2
         enddo
         call apply_field_tor_BC(field_tor_new, l, mm(j))
         ! multiply with inverted CN-matrices:
         ! This leaves all the fields in lmn space
         field_tor(l,j,1:Nr) = matmul(inv_lhs_IE_tor(1:Nr, 1:Nr, l),field_tor_new(1:Nr))
      jl_enddo
      rhs_IE_tor_old = rhs_IE_tor
      !-- all fields are given in (lmn) space
      !-- new derivatives are in (lmr).
      ! All fields are dealised in radial_dr_ddr_3D_n2r.
      call radial_dr_ddr_3D_n2r(field_tor, field_tor_dr, field_tor_ddr)
      call update_field_pol_lap()
      call update_field_tor_lap()
   end subroutine update_field

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
