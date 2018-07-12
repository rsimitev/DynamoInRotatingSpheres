module drs_renderers
#include "drsDefs.F90"
#define ERR_UNSUPPORTED_OPTION 33
   use drs_dims
   use drs_mpi
   use drs_params
   use drs_flow
   use drs_field
   use drs_temp
#ifdef COMP
   use drs_comp
#endif
   use drs_legendre
   use drs_transforms
   implicit none
   double precision, allocatable:: render_out(:,:,:)
   double precision, allocatable:: XX(:,:,:)
   double precision, allocatable:: YY(:,:,:)
   double precision, allocatable:: ZZ(:,:,:)

contains
   subroutine drs_renderers_allocation(what)
      implicit none
      integer, intent(in):: what
      if((what-10*(what/10)).eq.4) then
         if(.not.allocated(XX)) allocate(XX(0:Nt, Np, Nr))
         if(.not.allocated(YY)) allocate(YY(0:Nt, Np, Nr))
         if(.not.allocated(ZZ)) allocate(ZZ(0:Nt, Np, Nr))
      else
         if(.not.allocated(render_out)) allocate(render_out(0:Nt, Np, Nr))
      endif
   end subroutine drs_renderers_allocation

   !> Makes a decision about what to render.\n
   !! Numbers are coded as:\n
   !!~~~~~~
   !! a b c d e
   !! | | | | |>e - component 1, 2 or 3 for vectors, irrelevant for scalars
   !! | | | |> d - coordinate system or stream lines
   !! | | |> c - quantity to be ploted
   !! | |> b - curl, gradient or divergence or 0
   !! |> a - scalar product with selection or 0
   !!
   !!
   !!  e = 1, 2 or 3 for first second or third coordinate or meridional, azimuthal and poloidal streamlines
   !!      1 or 2 for total or anomaly scalar fiels
   !!      4 for all three coordinates
   !!  d = 1, 2 or 3 for cartesian (x,y,x), spherical (r,t,p) or cyllindrical (s, p, z) components respectively, 4 for streamlines, 0 for none
   !!  c = 1 for the flow
   !!      2 for the magetic field
   !!      3 for the temperature field
   !!      4 for the composition field
   !!      5 for the magetic field outside the core (up to ro+1)
   !!  b = 1 for the curl
   !!      2 for the gradient
   !!      3 for the divergence
   !!      0 for nothing
   !!  a = 1 for scalar product with flow
   !!      2 for scalar product with field
   !!      0 for nothing
   !!~~~~~~~
   !!
   !! For example, if I want the meridional (spherical coordinates) component of the curl of the flow,
   !! a=0, b=1, c=1, d=2, e=2 so @param what = 01122
   subroutine render(what)
      implicit none
      integer, intent(in)::what
      select case (what)
         case(113) ! z-velocity
            call render_uz()
         case(121) ! radial velocity:
            call render_ur()
         case(122) ! meridional velocity
            call render_ut()
         case(123) ! azimuthal velocity
            call render_up()
         case(124)
            call render_u()
         case(141) !-- flow meridional streamlines r*sinth* d/dth v (veloc.)
            call render_streamlines_t()
         case(143) !-- flow poloidal streamlines
            call render_poloidal_streamlines()
         case(213) ! z-field
            call render_Bz()
         case(221) ! radial field:
            call render_Br()
         case(222) ! meridional field
            call render_Bt()
         case(223) ! azimuthal field
            call render_Bp()
         case(224)
            call render_B()
         case(301) ! temperature
            call render_temperature()
         case(302) !-- temperature perturbation:
            call render_temperature_perturbation()
#ifdef COMP
         case(401) ! composition
            call render_composition()
         case(402) !-- composition perturbation:
            call render_composition_perturbation()
#endif
         case(524)
            call render_B_outside()
         case(1121) ! radial vorticity:
            call render_rotu_r()
         case(1122)
            call render_rotu_t()
         case(1123) !-- omega phi:
            call render_rotu_p()
         case(1124) !-- omega phi:
            call render_rotu()
         case(1113) ! z vorticity
            call render_rotu_z()
         case(11100) !-- helicity: u.curl(u)
            call render_helicity()
         case default
            spew "We do not support this option yet: ", what
            stop ERR_UNSUPPORTED_OPTION
      end select
   end subroutine

   subroutine render_ur()
      implicit none
      double precision:: ut(0:Nt,Np,Nr)
      double precision:: up(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_flow(render_out, ut, up )
   end subroutine
   subroutine render_u()
      implicit none
      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_flow(XX,YY,ZZ)
   end subroutine

   subroutine render_Br()
      implicit none
      double precision:: Bt(0:Nt,Np,Nr)
      double precision:: Bp(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(field_pol, field_pol_dr, field_pol_ddr)
      call radial_dr_ddr_3D_r2r(field_tor, field_tor_dr, field_tor_ddr)
      call calc_field(render_out, Bt, Bp )
   end subroutine
   subroutine render_Bt()
      implicit none
      double precision:: Br(0:Nt,Np,Nr)
      double precision:: Bp(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(field_pol, field_pol_dr, field_pol_ddr)
      call radial_dr_ddr_3D_r2r(field_tor, field_tor_dr, field_tor_ddr)
      call calc_field(Br, render_out, Bp )
   end subroutine

   !> Renders the azimuthal component of the magnetic field
   subroutine render_Bp()
      implicit none
      double precision:: Bt(0:Nt,Np,Nr)
      double precision:: Br(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(field_pol, field_pol_dr, field_pol_ddr)
      call radial_dr_ddr_3D_r2r(field_tor, field_tor_dr, field_tor_ddr)
      call calc_field(Br, Bt,render_out )
   end subroutine

   !> Renders the z component of the magnetic field.
   subroutine render_Bz()
      implicit none
      double precision:: Br(0:Nt,Np,Nr)
      double precision:: Bt(0:Nt,Np,Nr)
      double precision:: Bp(0:Nt,Np,Nr)
      double precision:: c, c2
      integer:: l

      call radial_dr_ddr_3D_r2r(field_pol, field_pol_dr, field_pol_ddr)
      call radial_dr_ddr_3D_r2r(field_tor, field_tor_dr, field_tor_ddr)
      call calc_field(Br, Bt, Bp)

      do l=0, Nt
         c  = costheta(l)
         c2 = cos(acos(costheta(l)) + pi/2)
         render_out(l,:,:) = Br(l,:,:)*c + Bt(l,:,:)*c2
      enddo
   end subroutine

   !> Render all three spherical components of the magnetic field.
   subroutine render_B()
      implicit none
      call radial_dr_ddr_3D_r2r(field_pol, field_pol_dr, field_pol_ddr)
      call radial_dr_ddr_3D_r2r(field_tor, field_tor_dr, field_tor_ddr)
      call calc_field(XX,YY,ZZ)
   end subroutine

   !> Render all three spherical components of the magnetic field outside the
   !! outer core.
   subroutine render_B_outside()
      implicit none
      integer:: l,j,i
      double precision:: rout, rin
      double precision:: rcoll_old(Nr), rcoll2_old(Nr), drcoll_old(Nr)
      ! Save the old radioal coordinates
      rcoll_old  = rcoll
      rcoll2_old = rcoll2
      drcoll_old = drcoll
      ! Redefine the radial coordinates
      rin  = rcoll(Nr)
      rout = rcoll(1)
      rcoll = rout + 2.0d0*(rcoll - rin)
      do i=1, Nr
         rcoll2(i) = rcoll(i)**2
         if (i.gt.1) drcoll(i-1) = rcoll(i-1) - rcoll(i)
      enddo
      drcoll(Nr) = drcoll(Nr-1)
      ! reset the toroidal field
      field_tor     = 0.0d0
      field_tor_dr  = 0.0d0 
      field_tor_ddr = 0.0d0
      ! Upward continue the poloidal field
      forall(l=0:Nt_s, j=1:blk_ps_size(mpi_rank), i=1:Nr)
          field_pol(l,j,i) = field_pol(l,j,1)*(rout/rcoll(i))**l
      endforall

      call radial_dr_ddr_3D_r2r(field_pol, field_pol_dr, field_pol_ddr)
      call calc_field(XX,YY,ZZ)
      ! Restore the radial coordinates
      rcoll  = rcoll_old
      rcoll2 = rcoll2_old
      drcoll = drcoll_old
   end subroutine

   !> Renders the radial component of the curl of the flow.
   subroutine render_rotu_r()
      implicit none
      double precision:: rotu_t(0:Nt,Np,Nr)
      double precision:: rotu_p(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_rot_flow(render_out, rotu_t, rotu_p )
   end subroutine

   !> Renders all three spherical components of the curl of the flow (vorticity).
   subroutine render_rotu()
      implicit none
      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_rot_flow(XX, YY, ZZ)
   end subroutine

   !> u_phi:
   subroutine render_up()
      implicit none
      double precision:: ur(0:Nt,Np,Nr)
      double precision:: ut(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_flow(ur, ut, render_out)
   end subroutine

   !> rot(u)_phi:
   subroutine render_rotu_p()
      implicit none
      double precision:: rotu_r(0:Nt,Np,Nr)
      double precision:: rotu_t(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_rot_flow(rotu_r, rotu_t, render_out)
   end subroutine

   !> u_theta:
   subroutine render_ut()
      implicit none
      double precision:: ur(0:Nt,Np,Nr)
      double precision:: up(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_flow(ur, render_out, up)
   end subroutine

   !> rot(u)_theta:
   subroutine render_rotu_t()
      implicit none
      double precision:: rotu_r(0:Nt,Np,Nr)
      double precision:: rotu_p(0:Nt,Np,Nr)

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_rot_flow(rotu_r, render_out, rotu_p)
   end subroutine

   !> u_z
   subroutine render_uz()
      implicit none
      double precision:: ur(0:Nt,Np,Nr)
      double precision:: ut(0:Nt,Np,Nr)
      double precision:: up(0:Nt,Np,Nr)
      double precision:: c, c2
      integer:: l

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_flow(ur, ut, up)

      do l=0, Nt
         c  = costheta(l)
         c2 = cos(acos(costheta(l)) + pi/2)
         render_out(l,:,:) = ur(l,:,:)*c + ut(l,:,:)*c2
      enddo
   end subroutine

   !> rot(u)_z:
   subroutine render_rotu_z()
      implicit none
      double precision:: rotu_r_t(0:Nt,Np,Nr)
      double precision:: rotu_t_t(0:Nt,Np,Nr)
      double precision:: rotu_p_t(0:Nt,Np,Nr)
      double precision:: c, c2
      integer:: l

      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_rot_flow(rotu_r_t, rotu_t_t, rotu_p_t)
      !-- u,rotu are in (theta,phi, r) now
      do l=0, Nt
         c  = costheta(l)
         c2 = cos(acos(costheta(l))+pi/2)
         render_out(l,:,:) = rotu_r_t(l,:,:)*c + rotu_t_t(l,:,:)*c2
      enddo
   end subroutine

   !> Renders the temperature perturbation
   subroutine render_temperature_perturbation()
      implicit none
      !-- temp(l,m,r) -> render_out(theta,phi,r)
      call calc_temp(render_out)
   end subroutine

   !> Renders the total temperature.
   subroutine render_temperature()
      implicit none
      integer:: i, j, l
      !-- temp(l,m,r)
      call render_temperature_perturbation()
      !-- temp(theta,phi,r)
      forall(i=1:Nr,l=0:Nt,j=1:Np)
          render_out(l,j,i) = render_out(l,j,i) + temp_profile(i)
      endforall
   end subroutine

#ifdef COMP
   !> Renders the concentration perturbation
   subroutine render_composition_perturbation()
      implicit none
      !-- comp(l,m,r) -> render_out(theta,phi,r)
      call calc_comp(comp, render_out)
   end subroutine

   !> Renders the total concentration.
   subroutine render_composition()
      implicit none
      integer:: i, j, l
      !-- comp(l,m,r)
      call render_composition_perturbation()
      !-- render_out(theta,phi,r)
      forall(i=1:Nr,l=0:Nt,j=1:Np)
          render_out(l,j,i) = render_out(l,j,i) + comp_profile(i)
      endforall
   end subroutine
#endif

   !> Renders helicity
   subroutine render_helicity()
      implicit none
      double precision:: ur(0:Nt,Np,Nr)
      double precision:: ut(0:Nt,Np,Nr)
      double precision:: up(0:Nt,Np,Nr)
      double precision:: rotu_r(0:Nt,Np,Nr)
      double precision:: rotu_t(0:Nt,Np,Nr)
      double precision:: rotu_p(0:Nt,Np,Nr)
      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call calc_flow(ur, ut, up)
      call calc_flow(rotu_r, rotu_t, rotu_p)

      !-- u,rotu are in (r,theta,phi) now
      render_out(:,:,:) = ur(:,:,:)*rotu_r(:,:,:)+ &
                          ut(:,:,:)*rotu_t(:,:,:)+ &
                          up(:,:,:)*rotu_p(:,:,:)
   end subroutine

   subroutine render_temprature_grad_r()
      implicit none
      integer:: l,j
      call ylmb(temp, render_out)
      !-- temp(theta,phi,r)
      do l=0, Nt
         do j=1, Np
            render_out(l,j,1:Nr) = radial_derivative_r2r(render_out(l,j,1:Nr))
        enddo
      enddo
   end subroutine

   subroutine render_streamlines_t()
      implicit none
      double precision:: s
      integer:: l,j,i
      double precision:: aux(0:Nt,Np_s)
      do i=1, Nr
         do l=0, Nt
            s = sqrt(1.0d0-costheta(l)**2)
            do j=1, Np_s
               aux(l,j) = sum(rcoll(i)*s*flow_pol(m0*(j/2):Nt_s,j,i)*dleg(l,m0*(j/2):Nt_s,j))
            enddo
         enddo
         ! only transform m->phi because l->lat is already done in do loop:
         call m2phi_2d(aux, render_out(:,:,i))
      enddo
   end subroutine

   !> Renders the poloidal flow streamlines
   subroutine render_poloidal_streamlines()
      implicit none
      double precision:: aux(0:Nt,Np,Nr)
      integer:: m, l, j, i
      do i=1, Nr
         do l=0, Nt_s
            do m = m0, Np_s
               j = 2*(m/m0)
               aux(l,j,i)   =  m*flow_pol(l,j+1,i)*rcoll(i)
               aux(l,j+1,i) = -m*flow_pol(l,j,i)*rcoll(i)
            enddo
         enddo
      enddo
      !-- temp(r,l,m)
      call ylmb(aux, render_out)
      !-- temp(r,theta,phi)
   end subroutine

   !> radial stream function for the flow
   subroutine render_radial_streamfunction()
      implicit none
      call ylmb(flow_tor, render_out)
   end subroutine render_radial_streamfunction

end module drs_renderers

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab

