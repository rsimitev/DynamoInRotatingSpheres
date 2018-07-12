! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2015
!> Takes care of evolving the composition equation.
#include "drsDefs.F90"
module drs_composition_equation
   use drs_params
   use drs_dims
   use drs_flow
   use drs_comp
   use drs_time
   use drs_mpi
   use drs_hypDiff
   use drs_radial
   implicit none
   !> The Crank-Nicholson inverse operator for the composition.
   double precision, allocatable:: inv_lhs_CE(:,:,:)
   double precision, allocatable:: rhs_CE_old(:,:,:)
   !> CE for Composition Equation
   double precision, allocatable:: rhs_CE(:,:,:)
   logical, private:: initialised=.false.

contains

   !---------------------------------------------------------------
   !> Ensures all requirements to solve the composition equation are
   !! initialised.
   subroutine drs_composition_equation_init(error)
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

      call drs_comp_init(error)
      if(error.gt.0) return

      ! Need to make sure the flow is initialised
      if(flow_present) then
         call drs_flow_init(error)
         if(error.gt.0) return
      endif

      allocate(rhs_CE(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(rhs_CE_old(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(inv_lhs_CE(Nr,Nr,0:Nt_s))
      rhs_CE_old = 0.0d0
      initialised=.true.

      ! It generates matrices for the composition with factor Pc.
      if(comp_evolves) call cncomp_init(inv_lhs_CE, h, Pc, error)
      if(error.gt.0) return
   end subroutine

   !---------------------------------------------------------------
   !> Computes matrices for the Crank-Nicholson scheme on the Composition.
   !! @par dinv is the resulting matrix.
   !! @par dt is the timestep currently being used
   !! @par fact is a factor for the time step. In the current implementation
   !! this is the compositional Prandtl number.
   subroutine CNComp_init(dinv, dt, fac, error)
      implicit none
      double precision, target, intent(out)::dinv(Nr,Nr,0:Nt_s)
      double precision, intent(in):: dt,fac
      integer, intent(out):: error
      double precision, pointer:: aux(:,:)
      double precision:: dto2fac
      integer:: indx(Nr)
      integer:: i,j,l
      error=0

      if(.not.initialised) then
         error=1
         return
      endif

      dto2fac = dt/(2.0d0*fac)
      do l=0,Nt_s
         aux => dinv(1:Nr,1:Nr,l)
         do j=1,Nr ! Mode
            do i=1,Nr ! Radius
               ! The poly/SH representation of the laplacian
               aux(i,j) = - poly_ddr(i,j) + llp1(l)*poly(i,j)/rcoll2(i) - 2.0d0/rcoll(i)*poly_dr(i,j)
               call drs_apply_hypDiff(aux(i,j), l)
               ! The poly representation of the time derivative
               aux(i,j) = aux(i,j)*dto2fac + poly(i,j)
            enddo
         enddo

         ! implement boundary conditions: fixed values at boundaries
         if (compBC(1) == FixedComposition) aux(1, 1:Nr) = poly(1, 1:Nr)
         if (compBC(2) == FixedComposition) aux(Nr,1:Nr) = poly(Nr,1:Nr)
         ! Fixed flux at boundaries
         if (tempBC(1) == FixedChemFlux) aux(1, 1:Nr) = poly_dr(1, 1:Nr)
         if (tempBC(2) == FixedChemFlux) aux(Nr,1:Nr) = poly_dr(Nr,1:Nr)

         ! Dealiase
         if (Nr_s+1 .le. Nr) aux(:,Nr_s+1:Nr) = 0.d0

         call ludcmp(aux,Nr,Nr,indx)

         !-- invert all cn-matrixes for pol. flow:
         call matinv(aux, indx, Nr, Nr)
      enddo
   end subroutine CNComp_init

   !---------------------------------------------------------------
   !> Computes the composition change due to advection
   subroutine CompositionEquation(rhs_CE, error)
      implicit none
      double precision, intent(out):: rhs_CE(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      integer, intent(out):: error
      double precision:: ur_comp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision:: ut_comp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision:: up_comp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      integer:: t_end
      integer:: i,l,j
      error=0
      if(.not.initialised) then
         error=1
         return
      endif

      rhs_CE = 0.0d0
      t_end  = blk_t_size(mpi_rank)-1
      if (flow_present) then
         if(drs_calc_type.ne.LinearCompositionalOnset) then
            ur_comp = flow_r_t*comp_t
            ut_comp = flow_t_t*comp_t
            up_comp = flow_p_t*comp_t
            call vectorField2Divergence(ur_comp, ut_comp, up_comp, rhs_CE)
         endif

         forall (i=1:Nr,l=0:Nt_s,j=1:blk_ps_size(mpi_rank), l.ge.blk_ts_start(j) )
            rhs_CE(l,j,i) = -rhs_CE(l,j,i) - flow_pol(l,j,i)*llp1(l)/rcoll(i)*comp_profile_dr(i)
         endforall
      endif
   end subroutine CompositionEquation

   !---------------------------------------------------------------
   !> Applies boundary conditions and recomputes radial derivatives.
   subroutine just_apply_comp_BC()
      implicit none
      double precision:: comp_new(Nr)
      integer:: j, l
      jl_do(j,l)
         comp_new(1:Nr) = comp(l, j, 1:Nr)
         call apply_comp_BC(comp_new)
         comp(l,j,1:Nr) = comp_new(1:Nr)
      jl_enddo
      call radial_dr_ddr_3D_r2r(comp, comp_dr, comp_ddr)
   end subroutine

   !---------------------------------------------------------------
   !> Updates the composition using a hybrid Crank-Nicholson/Addams-Bashford implicit
   !! integration scheme.
   !! @param h is the time step;
   !! @param h1 and @param h2 are the Addams-Bashford weights;
   !! @param Pc is the compositional Prandtl number;
   subroutine update_comp(h, h1, h2, Pc, error)
      implicit none
      double precision, intent(in):: h, h1, h2, Pc
      integer, intent(out):: error
      double precision:: comp_new(Nr)
      double precision:: ho2Pc
      integer:: l, j, i
      error=0
      if(.not.initialised) then
         error=1
         return
      endif

      if(.not.comp_evolves) return

      ho2Pc = h/(2*Pc)
      jl_do(j,l)
         do i=1, Nr
            comp_new(i) = comp(l, j, i) + comp_lap(l, j, i)*ho2Pc + rhs_CE(l, j, i)*h1 - rhs_CE_old(l, j, i)*h2
         enddo
         call apply_comp_BC(comp_new)
         ! multiply with inverted CN-matrices:
         ! This leaves all the fields in lmn space
         comp(l,j,1:Nr) = matmul(inv_lhs_CE(1:Nr, 1:Nr, l), comp_new(1:Nr))
      jl_enddo
      rhs_CE_old = rhs_CE
      !-- all fields are given in (lmn) space
      !-- new derivatives are in (lmr).
      ! All fields are dealised in radial_dr_ddr_3D_n2r.
      if(drs_calc_type.eq.LinearCompositionalOnset) then
         comp(:,1,:) = 0.0d0
      endif
      call radial_dr_ddr_3D_n2r(comp, comp_dr, comp_ddr)
   end subroutine update_comp

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
