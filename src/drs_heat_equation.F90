! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2015
!> Takes care of evolving the heat equation.
#include "drsDefs.F90"
module drs_heat_equation
   use drs_params
   use drs_dims
   use drs_flow
   use drs_temp
   use drs_time
   use drs_mpi
   use drs_hypDiff
   use drs_radial
   implicit none
   !> The Crank-Nicholson inverse operator for the temperature.
   double precision, allocatable:: inv_lhs_TE(:,:,:)
   !> TE for Temperature Equation
   double precision, allocatable:: rhs_TE(:,:,:)
   double precision, allocatable:: rhs_TE_old(:,:,:)
   logical, private:: initialised=.false.

contains

   !---------------------------------------------------------------
   !> Ensures all requirements to solve the heat equation are
   !! initialised.
   subroutine drs_heat_equation_init(error)
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

      call drs_temp_init(error)
      if(error.gt.0) return

      ! Need to make sure the flow is initialised
      if(flow_present) then
         call drs_flow_init(error)
         if(error.gt.0) return
      endif

      allocate(rhs_TE(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(rhs_TE_old(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(inv_lhs_TE(Nr,Nr,0:Nt_s))
      rhs_TE_old = 0.0d0
      initialised=.true.

      ! It generates matrices for the temperature with factor Pt.
      if(temp_evolves) call cntemp_init(inv_lhs_TE, h, Pt, error)
      if(error.gt.0) return
   end subroutine

   !---------------------------------------------------------------
   !> Computes matrices for the Crank-Nicholson scheme on the Temperature.
   !! @par dinv is the resulting matrix.
   !! @par dt is the timestep currently being used
   !! @par fact is a factor for the time step. In the current implementation
   !! this is the thermal Prandtl number.
   subroutine CNTemp_init(dinv, dt, fac, error)
      implicit none
      double precision, target, intent(out):: dinv(Nr,Nr,0:Nt_s)
      double precision, intent(in):: dt, fac
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
         do j=1,Nr ! mode
            do i=1,Nr ! radius
               ! The poly/SH representation of the laplacian
               aux(i,j) = - poly_ddr(i,j) + llp1(l)*poly(i,j)/rcoll2(i) - 2.0d0/rcoll(i)*poly_dr(i,j)
               call drs_apply_hypDiff(aux(i,j), l)
               ! The poly representation of the time derivative
               aux(i,j) = aux(i,j)*dto2fac + poly(i,j)
            enddo
         enddo

         ! implement boundary conditions: fixed values at boundaries
         if (tempBC(1) == FixedTemperature) aux(1, 1:Nr) = poly(1, 1:Nr)
         if (tempBC(2) == FixedTemperature) aux(Nr,1:Nr) = poly(Nr,1:Nr)
         ! Fixed flux at boundaries
         if (tempBC(1) == FixedHeatFlux) aux(1, 1:Nr) = poly_dr(1, 1:Nr)
         if (tempBC(2) == FixedHeatFlux) aux(Nr,1:Nr) = poly_dr(Nr,1:Nr)

         ! Dealiase
         if (Nr_s+1 .le. Nr) aux(:,Nr_s+1:Nr) = 0.d0

         call ludcmp(aux,Nr,Nr,indx)

         !-- invert all cn-matrixes for tor. flow:
         call matinv(aux, indx, Nr, Nr)
      enddo
   end subroutine CNTemp_init

   !---------------------------------------------------------------
   !> Computes the temperature change due to advection.
   subroutine TemperatureEquation(rhs_TE, error)
      implicit none
      double precision, intent(out):: rhs_TE(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      integer, intent(out):: error
      double precision:: ur_temp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision:: ut_temp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision:: up_temp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      integer:: t_end
      integer:: i,l,j
      error=0
      if(.not.initialised) then
         error=1
         return
      endif

      rhs_TE = 0.0d0
      t_end  = blk_t_size(mpi_rank)-1
      if(flow_present) then
         if(drs_calc_type.ne.LinearThermalOnset) then
            ur_temp = flow_r_t*temp_t
            ut_temp = flow_t_t*temp_t
            up_temp = flow_p_t*temp_t
            call vectorField2Divergence(ur_temp, ut_temp, up_temp, rhs_TE)
         endif

         forall (i=1:Nr,l=0:Nt_s,j=1:blk_ps_size(mpi_rank), l.ge.blk_ts_start(j) )
            rhs_TE(l,j,i) = -rhs_TE(l,j,i) - flow_pol(l,j,i)*llp1(l)/rcoll(i)*temp_profile_dr(i)
         endforall
      endif
   end subroutine TemperatureEquation

   !---------------------------------------------------------------
   !> Applies boundary conditions and updates the derivatives.
   subroutine just_apply_temp_BC()
      implicit none
      double precision:: temp_new(Nr)
      integer:: j, l
      jl_do(j,l)
         temp_new(1:Nr) = temp(l, j, 1:Nr)
         call apply_temp_BC_RHS(temp_new)
         temp(l,j,1:Nr) = temp_new(1:Nr)
      jl_enddo
      call radial_dr_ddr_3D_r2r(temp, temp_dr, temp_ddr)
   end subroutine

   !---------------------------------------------------------------
   !> Updates the temperature using a hybrid Crank-Nicholson/Addams-Bashford implicit
   !! integration scheme.
   !! @param h is the time step;
   !! @param h1 and @param h2 are the Addams-Bashford weights;
   !! @param Pt is the Thermal Prandtl number;
   subroutine update_temp(h, h1, h2, Pt, error)
      implicit none
      double precision, intent(in):: h, h1, h2, Pt
      integer, intent(out):: error
      double precision:: temp_new(Nr)
      double precision:: ho2Pt
      integer:: l, j, i
      error=0
      if(.not.initialised) then
         error=1
         return
      endif

      if(.not.temp_evolves) return

      ho2Pt = h/(2*Pt)
      jl_do(j,l)
         do i=1, Nr
            temp_new(i) = temp(l, j, i) + temp_lap(l, j, i)*ho2Pt + rhs_TE(l, j, i)*h1 - rhs_TE_old(l, j, i)*h2
         enddo
         call apply_temp_BC_RHS(temp_new)
         ! multiply with inverted CN-matrices:
         ! This leaves all the fields in lmn space
         temp(l,j,1:Nr) = matmul(inv_lhs_TE(1:Nr, 1:Nr, l), temp_new(1:Nr))
      jl_enddo
      rhs_TE_old = rhs_TE
      !-- all fields are given in (lmn) space
      !-- new derivatives are in (lmr).
      ! All fields are dealised in radial_dr_ddr_3D_n2r.
      if(drs_calc_type.eq.LinearThermalOnset) then
         temp(:,1,:) = 0.0d0
      endif
      call radial_dr_ddr_3D_n2r(temp, temp_dr, temp_ddr)
      call update_temp_lap()
   end subroutine update_temp
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
