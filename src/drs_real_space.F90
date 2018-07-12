! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Takes care of contructing the nonlinear terms of all equations and other
!! quantities in real space.
module drs_real_space
#include "drsDefs.F90"
   use drs_time
   use drs_dims
   use drs_mpi
   use drs_params
   use drs_flow
   use drs_temp
   use drs_radial
   use drs_probes
   use drs_debug
   implicit none
   save


contains

   !-----------------------------------------------------------------------
   !> Computes in real space the quantities that need to be computed in real space.
   subroutine evaluate_real_space()
      implicit none
      if(flow_present.and.flow_evolves) then
         !--     for the kinematic case the velocity field must be calculated
         !--     only at the first call of rhs()
         !--      in: /fields/t,p /derivatives/ flow_tor_dr,..,flow_pol_dr,flow_pol_ddr,..         (l,m,r)
         !--     out:         flow_r_t,flow_t_t,flow_p_t,..,rot_flow_p_t (theta,phi,r)
         call calc_u(flow_r_t, flow_t_t, flow_p_t, rot_flow_r_t, rot_flow_t_t, rot_flow_p_t)
         
         if(drs_calc_type.eq.KinematicDynamo) flow_evolves = .FALSE.
      endif

      if(field_present.and.field_evolves) then
         call calc_B(field_r_t, field_t_t, field_p_t, rot_field_r_t, rot_field_t_t, rot_field_p_t)
      endif

      if(temp_present.and.temp_evolves) call calc_temp(temp_t)

#ifdef COMP
      if(comp_present.and.comp_evolves) call calc_comp(comp, comp_t)
#endif
   end subroutine evaluate_real_space
   
   !--------------------------------------------------------------------------
   !> Computes the CFL numbers that we are going to use to set the time-step.\n
   !! The first three elements are just the Glatzmaier's CFL numbers (JCP 55 (1984), eq. 8)\n
   !! A fourth element was added to account for the Coriolis force.
   !! All Prandtl numbers are taken into account 
   subroutine drs_real_space_compute_cfl()
      implicit none
      integer, parameter::ncfl=4
      double precision:: aux
      integer:: i, j, l

      cfl = 1.D-20

      ! CFL for diffusion of the flow
      ! Radial 
      cfl(1) = 0.25d0*(1.0d0/Nr)**2

      ! CFL for the radial diffusion of the remaining quantities
      aux = 1.0d0
      if(temp_present  .and. Pt.lt.aux) aux = Pt
#ifdef COMP
      if(comp_present  .and. Pc.lt.aux) aux = Pc
#endif
      if(field_present .and. Pm.lt.aux) aux = Pm
      ! Chose the smallest Prandtl (that is the quantity with the highest
      ! diffusivity) to determine the diffusion step limit.
      cfl(1) = cfl(1)*aux

      ! CFL for the flow
      do i=1, Nr
         do j=1, Np
            do l=0, blk_t_size(mpi_rank)-1
               aux = abs( flow_r_t(l,j,i)/drcoll(i) )
               if (aux .gt. cfl(2)) cfl(2) = aux
               aux = sqrt( flow_t_t(l,j,i)**2 + flow_p_t(l,j,i)**2 )/rcoll(i)
               if (aux .gt. cfl(3)) cfl(3) = aux
            enddo
         enddo
      enddo
      cfl(2) = 1./cfl(2)
      cfl(3) = 1./( cfl(3)*sqrt( dble(Np*(Nt+1)) ) )

      ! Minimum step size for the Coriolis force
      cfl(4) = 0.5D0*sqrt(1.D0/Ta)

      cfl = cfl/2.0d0
      ! CPU0 collects the minimum...
      call drs_minimize(cfl)
      ! ...and sends it to everyone else.
      call drs_bcast(cfl, ncfl)
   end subroutine drs_real_space_compute_cfl
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
