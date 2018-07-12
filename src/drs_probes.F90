! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> This module implements some probing facilities for the running models.
module drs_probes
#include "drsDefs.F90"
   use drs_mpi
   use drs_dims
   use drs_time
   use drs_params
   use drs_flow
   use drs_field
   use drs_temp
#ifdef COMP
   use drs_comp
#endif
   use drs_hypDiff
   use drs_io_units
   use drs_error_codes
   implicit none
   save
   double precision, allocatable:: ur_avg(:), ut_avg(:), up_avg(:)
   double precision, allocatable:: up2(:), ut2(:)
   double precision, allocatable:: adv_avg(:), t2_avg(:)
   double precision, allocatable:: tspec_avg(:), urspec_avg(:)
   double precision, allocatable:: Brspec_avg(:)
   ! Growth rate
   double precision:: groth,Ekin, EB
   double precision, allocatable:: dOmega(:,:) !< Weights for volume integration.
   double precision:: Rm = 1.0d0 !< The magnetic Reynolds number.

contains

   !---------------------------------------------------------------------------
   !> Allocate quantities to be probed
   subroutine drs_probes_allocation()
      implicit none
      allocate( ur_avg(Nr), ut_avg(Nr), up_avg(Nr) )
      allocate( adv_avg(Nr), t2_avg(Nr))
      allocate( tspec_avg(0:Nt_s), urspec_avg(0:Nt_s) )
      allocate( Brspec_avg(0:Nt_s) )
      allocate( up2(Nr), ut2(Nr) )
      allocate( dOmega(0:blk_t_size(mpi_rank)-1, Nr) )
   end subroutine

   !---------------------------------------------------------------------------
   !> Initialise arrays and integration factors.
   subroutine drs_probes_init(time)
      implicit none
      double precision:: time
      integer:: i,l, lg

      ur_avg = 0.0d0
      ut_avg = 0.0d0
      up_avg = 0.0d0
      adv_avg  = 0.0d0
      t2_avg = 0.0d0
      ut2 = 0.0d0
      up2 = 0.0d0

      tspec_avg  = 0.0d0
      urspec_avg = 0.0d0

      if(field_evolves) Brspec_avg = 0.0d0

      dOmega = 0.0d0
      do l=0, blk_t_size(mpi_rank)-1
         lg = blk_t_start(mpi_rank) + l
         dOmega(l,1) = 2*pi/Np*rcoll2(1)*w(lg)*0.5*(rcoll(1)-rcoll(2))
         do i=2,Nr-1
            dOmega(l,i) = 2*pi/Np*rcoll2(i)*w(lg)*0.5*(rcoll(i-1)-rcoll(i+1))
         enddo
         dOmega(l,Nr) = 2*pi/Np*rcoll2(Nr)*w(lg)*0.5*(rcoll(Nr-1)-rcoll(Nr))
      enddo

      call update_time_last_sample(time)
   end subroutine

   !-----------------------------------------------------------------------------
   !> Performs the integration in theta and phy of a function f raised to the
   !! power n.
   !! as \f[ F(r) = \int\int f(r,\theta,\phi)^n \sin\theta d\theta d\phi \f]
   subroutine integrate_power_surf(f, n, f_int)
      implicit none
      !> The function to be integrated.
      double precision, intent(in):: f(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      !> The power it should be raised to.
      integer, intent(in):: n
      !> Integral as a function of r.
      double precision, intent(out):: f_int(Nr)
      integer:: i, j
      integer:: Ntl, Nts

      Ntl = blk_t_size(mpi_rank)-1
      Nts = blk_t_start(mpi_rank)

      if(n==1) then
         do i=1, Nr
            f_int(i) = 0.d0
            do j=1, Np
               f_int(i) = f_int(i) + sum(w(Nts:Nts+Ntl)*f(0:Ntl,j,i))
            enddo
         enddo
      else
         do i=1, Nr
            f_int(i) = 0.d0
            do j=1, Np
               f_int(i) = f_int(i) + sum(w(Nts:Nts+Ntl)*f(0:Ntl,j,i)**n)
            enddo
         enddo
      endif
      call sum_over_all_cpus(f_int)
      f_int = f_int/(2*Np)
   end subroutine

   !-----------------------------------------------------------------------------
   !> Computes the total integrated energy of a vector field based on its components.\n
   !! Only root contains the solution.
   double precision function energy(vr,vt,vp)
      implicit none
      double precision, intent(in):: vr(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(in):: vt(0:blk_t_size(mpi_rank)-1,Np,Nr)
      double precision, intent(in):: vp(0:blk_t_size(mpi_rank)-1,Np,Nr)
      integer::j, Ntl

      Ntl = blk_t_size(mpi_rank)-1
      energy  = 0.d0
      do j=1, Np
         energy = energy + sum( (vr(0:Ntl,j,1:Nr)**2 + &
                                 vt(0:Ntl,j,1:Nr)**2 + &
                                 vp(0:Ntl,j,1:Nr)**2)*dOmega(0:Ntl,1:Nr) )
      enddo
      energy = 0.5d0*energy
      call sum_over_all_cpus(energy)
   end function energy

   !--------------------------------------------------------------------------
   !> Encapsulate saving quantities in real and spectral space.
   subroutine save_stuff(nsample)
      implicit none
      integer, intent(inout):: nsample

      time_since_last_sample = time - time_last_sample
      call measure_lm()
      !-- (p,t,temp,field_tor,field_pol are not mod., they are still in (l,m,r) space)
      !-- measure: compute some integrals in direct space:
      !-- attention: measure() needs flow_r_t,flow_t_t,..,temp_t
      !-- (so the non-blocking comm. has to be finished here!)
      call measure(temp_t,&
                  rot_flow_r_t, rot_flow_t_t, rot_flow_p_t)
      call update_time_last_sample(time)
      nsample = nsample + 1
   end subroutine

   !-----------------------------------------------------------------------------
   !> Computes the components of the kinetic energy
   !! energies: equations 13,14,15,16 of Sun, Schubert, Glatzmaier 1992, p. 99
   !! compute the non-axisymmetric energies:
   !! here we have to take into account the normalization factors,
   !! because we are interested in the real energies!\n
   !! pol   =   P\n
   !! tor   =   T/r
   subroutine compute_KE_components(nkes, nkea, etors, etora, drkes, drkea, mckes, mckea)
      implicit none
      double precision, intent(out):: nkes
      double precision, intent(out):: nkea
      double precision, intent(out):: etors
      double precision, intent(out):: etora
      double precision, intent(out):: drkes
      double precision, intent(out):: drkea
      double precision, intent(out):: mckes
      double precision, intent(out):: mckea
      double precision:: energs(Nr)
      double precision:: energa(Nr)
      double precision:: ztors(Nr)
      double precision:: ztora(Nr)
      integer:: i, j, l, m, jg, jj
      double precision:: w1,wr,z,r,r2,p11,p22


      do i=1,Nr
         energs(i)=0.d0
         energa(i)=0.d0
         ztors(i)=0.d0
         ztora(i)=0.d0
         r  = rcoll(i)
         r2 = rcoll2(i)
         do l=1,Nt_s
            do m=1,l
               jg = 2*(m/m0)  !  global index of real part
               jj = jg-blk_ps_start(mpi_rank)+1   !   local index of real part
               if(jj.ge.1 .and. jj+1.le.blk_ps_size(mpi_rank) .and. mod(m,m0).eq.0) then
                  !-- W**2 = (rP)**2 = (rp)**2:
                  p11 = plmfac(l,m)*flow_pol(l,jj,i)
                  p22 = plmfac(l,m)*flow_pol(l,jj+1,i)
                  w1  = (p11**2 + p22**2)*2.d0*pi*r2
                  !-- (d/dr W)**2 = (d/dr(rP))**2 = (r*d/dr(flow_pol) + flow_pol)**2:
                  p11 = plmfac(l,m)*(flow_pol_dr(l,jj,i)*r   + flow_pol(l,jj,i))
                  p22 = plmfac(l,m)*(flow_pol_dr(l,jj+1,i)*r + flow_pol(l,jj+1,i))
                  wr  = (p11**2 + p22**2)*2.d0*pi
                  !-- Z**2 = (rT)**2 = (r**2*flow_tor)**2:
                  p11 = plmfac(l,m)*flow_tor(l,jj,i)
                  p22 = plmfac(l,m)*flow_tor(l,jj+1,i)
                  z   = (p11**2 + p22**2)*2.d0*pi*r2**2
                  !-- eq. (13):
                  if(mod(l+m,2) .eq. 0) then
                     energs(i) = energs(i) + (llp1(l)/r2*w1 + wr + z)*llp1(l)
                     ztors(i)  = ztors(i)  + z*llp1(l)
                  else
                     energa(i) = energa(i) + (llp1(l)/r2*w1 + wr + z)*llp1(l)
                     ztora(i)  = ztora(i)  + z*llp1(l)
                  endif
               endif
            enddo  !  m=1,l
         enddo  !  l=1,Nt_s
      enddo   !  i=1,Nr
      !-- nke : nonaxisymmetric kinetic energy
      !-- etor: toroidal nonaxisymmetric energy
      nkes  = integrate_r(energs)
      nkea  = integrate_r(energa)
      etors = integrate_r(ztors)
      etora = integrate_r(ztora)
      call sum_over_all_cpus(nkes)
      call sum_over_all_cpus(nkea)
      call sum_over_all_cpus(etors)
      call sum_over_all_cpus(etora)

      !-- drke: axisymetric (m=0) differential rotation kinetic energy
      jg = 1
      m  = m0*(jg/2)
      j  = jg - blk_ps_start(mpi_rank) + 1
      if(j.eq.1) then  !  true only for mpi_rank=0.
         energs = 0.d0
         energa = 0.d0
         do l=1,Nt_s
            !-- Z0**2 = (rT)**2 = (r**2*flow_tor)**2:
            if(mod(l+m,2) .eq. 0) then
               do i=1,Nr
                  energs(i) = energs(i) + flow_tor(l,1,i)**2*2.d0*pi/(2*l+1)*llp1(l)*rcoll(i)**4
               enddo
            else
               do i=1,Nr
                  energa(i) = energa(i) + flow_tor(l,1,i)**2*2.d0*pi/(2*l+1)*llp1(l)*rcoll(i)**4
               enddo
            endif
         enddo
         drkes = integrate_r(energs)
         drkea = integrate_r(energa)
      endif

      !-- mcke: mean (m=0) meridional circulation kinetic energy
      jg = 1
      m  = m0*(jg/2)
      j  = jg - blk_ps_start(mpi_rank) + 1
      if(j.eq.1) then  !  true only for mpi_rank=0.
         energs=0.d0
         energa=0.d0
         do l=1,Nt_s
            !-- l(l+1)*[ l(l+1)/r**2*W0**2 + (d/drW0)**2] =
            !-- l(l+1)*[ l(l+1)*flow_pol**2       + (r*d/drp + flow_pol)**2]
            if(mod(l+m,2) .eq. 0) then
            do i=1,Nr
               energs(i) = energs(i) + &
                           2.d0*pi/(2*l+1)*llp1(l)*( llp1(l)*flow_pol(l,1,i)**2 + &
                           (flow_pol_dr(l,1,i)*rcoll(i) + flow_pol(l,1,i))**2 )
            enddo
            else
            do i=1,Nr
               energa(i) = energa(i) + &
                           2.d0*pi/(2*l+1)*llp1(l)*( llp1(l)*flow_pol(l,1,i)**2 + &
                           (flow_pol_dr(l,1,i)*rcoll(i)+flow_pol(l,1,i))**2 )
            enddo
            endif
         enddo
            mckes = integrate_r(energs)
            mckea = integrate_r(energa)
      endif
   end subroutine

   !-----------------------------------------------------------------------------
   !>computes the components of the magnetic energy.
   subroutine compute_BE_components(Bnkes, Bnkea, Betors, Betora, Bdrkes, Bdrkea, Bmckes, Bmckea)
      implicit none
      double precision, intent(out):: Bnkes
      double precision, intent(out):: Bnkea
      double precision, intent(out):: Betors
      double precision, intent(out):: Betora
      double precision, intent(out):: Bdrkes
      double precision, intent(out):: Bdrkea
      double precision, intent(out):: Bmckes
      double precision, intent(out):: Bmckea
      double precision:: energs(Nr)
      double precision:: energa(Nr)
      double precision:: ztors(Nr)
      double precision:: ztora(Nr)
      integer:: i, j, l, m, jg, jj
      double precision:: w1,wr,z,r2,p11,p22

      !-- Bnke : all nonaxisymmetric energy
      !-- Betor: toroidal nonaxisymmetric energy
      do i=1,Nr
         energs(i) = 0.d0
         energa(i) = 0.d0
         ztors(i)  = 0.d0
         ztora(i)  = 0.d0
         r2 = rcoll2(i)
         do l=1,Nt_s
            do m=1,l
               jg = 2*(m/m0)  !  global index of real part
               jj = jg - blk_ps_start(mpi_rank)+1   !   local index of real part
               if(jj.ge.1 .and. jj+1.le.blk_ps_size(mpi_rank) .and. mod(m,m0).eq.0) then
                  !-- W**2 = (rH)**2 = field_pol**2:
                  p11 = plmfac(l,m)*field_pol(l,jj,i)
                  p22 = plmfac(l,m)*field_pol(l,jj+1,i)
                  w1  = (p11*p11+p22*p22)*2.d0*pi
                  !-- (d/drW)**2 = (d/dr(rH))**2 = (d/dr(field_pol))**2:
                  p11 = plmfac(l,m)*field_pol_dr(l,jj,i)
                  p22 = plmfac(l,m)*field_pol_dr(l,jj+1,i)
                  wr  = (p11*p11 + p22*p22)*2.d0*pi
                  !-- Z**2 = (rG)**2 = field_tor**2:
                  p11 = plmfac(l,m)*field_tor(l,jj,i)
                  p22 = plmfac(l,m)*field_tor(l,jj+1,i)
                  z   = (p11*p11+p22*p22)*2.d0*pi
                  !-- eq. (13):
                  if(mod(l+m,2) .eq. 0) then
                     energs(i) = energs(i) + (llp1(l)/r2*w1+wr+z)*llp1(l)
                     ztors(i)  = ztors(i)  + z*llp1(l)
                  else
                     energa(i) = energa(i) + (llp1(l)/r2*w1+wr+z)*llp1(l)
                     ztora(i)  = ztora(i)  + z*llp1(l)
                  endif
               endif
            enddo  !  m=1,l
         enddo  !  l=1,Nt_s
      enddo
      Bnkes  = integrate_r(energs)
      Bnkea  = integrate_r(energa)
      Betors = integrate_r(ztors)
      Betora = integrate_r(ztora)
      call sum_over_all_cpus(Bnkes)
      call sum_over_all_cpus(Bnkea)
      call sum_over_all_cpus(Betors)
      call sum_over_all_cpus(Betora)

      !-- Bdrke: axisymetric (m=0) differential rotation kinetic energy
      jg=1
      m = m0*(jg/2)
      j = jg - blk_ps_start(mpi_rank) + 1
      if(j.eq.1) then  !  true only for mpi_rank=0.
         do i=1, Nr
            energs(i) = 0.d0
            energa(i) = 0.d0
            do l=1,Nt_s
               !-- Z0**2 = (rG)**2 = field_tor**2:
               if(mod(l+m,2) .eq. 0) then
                  energs(i) = energs(i) + field_tor(l,1,i)**2*4.d0*pi/dble(2*l+1)*llp1(l)*0.5d0
               else
                  energa(i) = energa(i) + field_tor(l,1,i)**2*4.d0*pi/dble(2*l+1)*llp1(l)*0.5d0
               endif
            enddo
         enddo
         Bdrkes = integrate_r(energs)
         Bdrkea = integrate_r(energa)
      endif

      !-- Bmcke: mean (m=0) meridional circulation kinetic energy
      jg=1
      m = m0*(jg/2)
      j = jg - blk_ps_start(mpi_rank)+1
      if(j.eq.1) then
         do i=1, Nr
            energs(i) = 0.d0
            energa(i) = 0.d0
            do l=1,Nt_s
               !-- l(l+1)*[ l(l+1)/r**2*W0**2 + (d/drW0)**2] =
               !-- l(l+1)[ l(l+1)H**2      + (d/dr(rH))**2 ] =
               !-- l(l+1)[ l(l+1)(field_pol/r)**2 + (d/dr(field_pol))**2 ]
               if(mod(l+m,2) .eq. 0) then
                  energs(i) = energs(i) + 2.d0*pi/(2*l+1)*llp1(l)*&
                              (llp1(l)*(field_pol(l,1,i)/rcoll(i))**2+field_pol_dr(l,1,i)**2)
               else
                  energa(i) = energa(i) + 2.d0*pi/(2*l+1)*llp1(l)*&
                              (llp1(l)*(field_pol(l,1,i)/rcoll(i))**2+field_pol_dr(l,1,i)**2)
               endif
            enddo
         enddo
         Bmckes = integrate_r(energs)
         Bmckea = integrate_r(energa)
      endif
   end subroutine

   !-----------------------------------------------------------------------------
   !> Measures and saves quantities of interest (in (l,m,r)).
   !! All inputs have to be in (l,m,r) space
   subroutine measure_lm()
      implicit none
      integer:: i,j,l,mmax

      mmax  = m0*(Np_s/2)

      if(temp_evolves) then
         !-- get average radial temperature gradient
         temp_dr_avg(1:Nr) = temp_dr_avg(1:Nr) + temp_dr(0,1,1:Nr)*time_since_last_sample

         !-- get time average of temperature field
         do i=1,Nr
            jl_do(j,l)
               temp_avg(l,j,i) = temp_avg(l,j,i) + temp(l,j,i)*time_since_last_sample
            jl_enddo
         enddo

         !-- get the temperature l-power spectrum in adimensional units and its time average
         call average_unnormalised_scalar_l_spectrum(temp, tspec_avg)
      endif

      if(flow_evolves) then
         !-- get time average of toroidal, poloidal flows
         do i=1,Nr
            jl_do(j,l)
               flow_tor_avg(l,j,i) = flow_tor_avg(l,j,i) + flow_tor(l,j,i)*time_since_last_sample
               flow_pol_avg(l,j,i) = flow_pol_avg(l,j,i) + flow_pol(l,j,i)*time_since_last_sample
            jl_enddo
         enddo
         !-- get the l power spectrum of radial velocity
         call average_unnormalised_flow_l_spectrum(urspec_avg)
         call save_flow_dissipation(mmax)
         call save_flow_coeffs()
      endif

      if(field_evolves) then
         !-- get the l power spectrum of radial magnetic field
         call average_unnormalised_field_l_spectrum(Brspec_avg)
         call save_magnetic_dissipation(mmax)
         call save_field_coeffs()
      endif
   end subroutine measure_lm

   !---------------------------------------------------------------------------
   !> Computes the unnormalised flow l-spectrum and recomputes the average.
   subroutine average_unnormalised_flow_l_spectrum(urspec_avg)
      implicit none
      double precision, intent(inout):: urspec_avg(0:Nt_s)
      double precision:: spec
      integer:: l, m, i, j, jg, mmax
      mmax  = m0*(Np_s/2)
      do l=0, Nt_s
         spec=0.d0
         do m=0,min(l,mmax),m0
            do i=1,Nr-1
               jg=2*(m/m0)+1
               j = jg - blk_ps_start(mpi_rank) + 1
               if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
                  spec = spec + ( (llp1(l)/rcoll2(i)*  flow_pol(l,j,i))**2 + &
                                  (llp1(l)/rcoll2(i+1)*flow_pol(l,j,i+1))**2 )*0.5*drcoll(i)
               endif
               if(m.gt.0) then
                  j = j - 1
                  if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
                     spec = spec + ( (llp1(l)/rcoll2(i)*  flow_pol(l,j,i))**2 + &
                                     (llp1(l)/rcoll2(i+1)*flow_pol(l,j,i+1))**2 )*0.5*drcoll(i)
                  endif
               endif
            enddo
         enddo
         urspec_avg(l) = urspec_avg(l) + spec*time_since_last_sample
      enddo
   end subroutine average_unnormalised_flow_l_spectrum

   !---------------------------------------------------------------------------
   !> Computes the unnormalised field l-spectrum and recomputes the average.
   subroutine average_unnormalised_field_l_spectrum(Brspec_avg)
      implicit none
      double precision, intent(inout):: Brspec_avg(0:Nt_s)
      double precision:: spec
      integer:: l, m, i, j, jg, mmax
      mmax  = m0*(Np_s/2)
      do l=0,Nt_s
         spec=0.d0
         do m=0,min(l,mmax),m0
            do i=1,Nr-1
               jg = 2*(m/m0)+1
               j  = jg - blk_ps_start(mpi_rank) + 1
               if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
                  spec = spec + ( (llp1(l)/rcoll(i)**4*  field_pol(l,j,i))**2 + &
                                  (llp1(l)/rcoll(i+1)**4*field_pol(l,j,i+1))**2 )*0.5*drcoll(i)
               endif
               if (m.gt.0) then
                  j=j-1
                  if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
                     spec = spec + ( (llp1(l)/rcoll(i)**4*  field_pol(l,j,i) )**2 + &
                                     (llp1(l)/rcoll(i+1)**4*field_pol(l,j,i+1))**2 )*0.5*drcoll(i)
                  endif
               endif
            enddo
         enddo
         Brspec_avg(l) = Brspec_avg(l) + spec*time_since_last_sample
      enddo
   end subroutine average_unnormalised_field_l_spectrum

   !---------------------------------------------------------------------------
   !> Computes the unnormalised l-spectrum of a scalar field and recomputes the average.
   subroutine average_unnormalised_scalar_l_spectrum(scalar, scalar_spec_avg)
      implicit none
      double precision, intent(in):: scalar(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      double precision, intent(inout):: scalar_spec_avg(0:Nt_s)
      double precision:: spec
      integer:: l, m, i, j, jg, mmax
      mmax  = m0*(Np_s/2)
      do l=0,Nt_s
         spec=0.d0
         do m=0,min(l,mmax),m0
            do i=1,Nr-1
               jg = 2*(m/m0) + 1
               j  = jg - blk_ps_start(mpi_rank) + 1
               if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
                  spec = spec + ( scalar(l,j,i)**2 + scalar(l,j,i+1)**2 )*0.5*drcoll(i)
               endif
               if (m.gt.0) then
                  j=j-1
                  if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
                     spec = spec + ( scalar(l,j,i)**2 + scalar(l,j,i+1)**2 )*0.5*drcoll(i)
                  endif
               endif
            enddo
         enddo
         scalar_spec_avg(l) = scalar_spec_avg(l) + spec*time_since_last_sample
      enddo
   end subroutine average_unnormalised_scalar_l_spectrum

   !-----------------------------------------------------------------------------
   !> Computes the normalized power spectrum with respect to l of a scalar field.
   !! @param field is in lmr space.
   !! @param spec the l spectrum of the scalar field \a field.
   subroutine l_spec_of_scalar_field(field,spec)
      implicit none
      double precision, intent(out):: spec(0:Nt_s)
      double precision, intent(in):: field(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      integer:: i,j,jg,l,m,mmax
      double precision:: plm
      integer:: blk_size,blk_start

      blk_size = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)
      mmax = m0*(Np_s/2)

      spec = 0.0
      do l=0,Nt_s
         do m=0,min(l,mmax),m0
            plm = plmfac(l,m)
            do i=1,Nr-1
               jg=2*(m/m0)+1
               j=jg-blk_start+1
               if(j.ge.1 .and. j.le.blk_size) then
                  spec(l) = spec(l) + ( &
                     ( plm*field(l,j,i)   )**2 + &
                     ( plm*field(l,j,i+1) )**2   &
                     )*0.5* drcoll(i)
               endif
               if (m.gt.0) then
                  j=j-1
                  if(j.ge.1 .and. j.le.blk_size) then
                     spec(l) = spec(l) + ( &
                        ( plm*field(l,j,i)   )**2 + &
                        ( plm*field(l,j,i+1) )**2   &
                        )*0.5* drcoll(i)
                  endif
               endif
            enddo
         enddo
      enddo
      call sum_over_all_cpus(spec)
   end subroutine l_spec_of_scalar_field

   !-----------------------------------------------------------------------------
   !>  Calculates the normalized power spectrum of a scalar field with respect to m.
   !! @param field is in lmr space.
   !! @param spec the m spectrum of the scalar field \a field.
   subroutine m_spec_of_scalar_field(field, spec)
      implicit none
      double precision, intent(in):: field(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(out):: spec(m0*Np_s+1)
      integer:: i,j,jg,l,m,mmax
      double precision:: plm
      integer:: blk_size,blk_start

      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)

      mmax = m0*(Np_s/2)
      spec = 0.0

      do i=1,Nr-1
         do m=0, mmax, m0
            jg = 2*(m/m0)+1
            j  = jg - blk_start + 1
            if(j.gt.0 .and. j.le.blk_size) then
               do l=m,Nt_s
                  plm = plmfac(l,m)
                  !-- imag. part
                  if(j.le.blk_size) then
                     spec(jg) = spec(jg) + ( &
                        (plm*field(l,j,i))**2 + (plm*field(l,j,i+1))**2 &
                        )*0.5*drcoll(i)
                  endif      !  (j.le.blk_size)
                  if (m.gt.0) then
                     !-- add the real part
                     spec(jg) = spec(jg) + ( &
                        (plm*field(l,j,i))**2 + (plm*field(l,j,i+1))**2 &
                        )*0.5*drcoll(i)
                  endif
               enddo
            endif
         enddo
      enddo
      call sum_over_all_cpus(spec)
   end subroutine

   !-----------------------------------------------------------------------------
   !> Calculates the normalized power spectrum of a scalar quantity \f$f\f$ with respect
   !! to the Chebyshev polynomials.
   !! \f[ R_n = \sum_{l,m} N_l^m ({f_{nl}^m})^2 \f]
   !! @param field is in lmr space.
   !! @param spec the n spectrum of the scalar field \a field.
   subroutine n_spec_of_scalar_field(field, spec)
      implicit none
      double precision, intent(in):: field(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, intent(out):: spec(Nr_s)
      integer:: i,j,jg,l,m,mmax
      double precision:: radarr(Nr)
      double precision:: aux(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      integer:: blk_size,blk_start

      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)
      mmax = m0*(Np_s/2)

      !-- power spectrum of Chebyshev coefficients
      spec = 0.0
      do j=1,blk_size
         m = m0*((j+blk_start)/2)
         do l=m, Nt_s
            !c----- temperature field
            radarr(1:Nr) = field(l,j,1:Nr)
            call cos_r2r_1_r2n(radarr)
            aux(l,j,1:Nr) = radarr(1:Nr)
         enddo
      enddo
      do i=1,Nr_s
         do m=0, mmax, m0
            jg = 2*(m/m0) + 1
            j = jg - blk_start + 1
            if(j.gt.0 .and. j.le.blk_size) then
               do l=m, Nt_s
                  !-- imag. part
                  if(j.le.blk_size) then
                     spec(i) = spec(i) + (plmfac(l,m)*aux(l,j,i))**2
                  endif      !  (j.le.blk_size)
                  if (m.gt.0) then
                     !-- add the real part
                     if (j.gt.1) &
                       spec(i) = spec(i) + (plmfac(l,m)*aux(l,j-1,i))**2
                  endif
               enddo         !  l=m,Nt_s
            endif            !  jj.gt.0
         enddo               !  m=0,mmax
      enddo                  !  i=1,Nr
      call sum_over_all_cpus(spec)
   end subroutine n_spec_of_scalar_field

   !-----------------------------------------------------------------------------
   !> Performs the integration of the 1d real array @par input in the radial direction.
   pure function integrate_r(input) result(aux)
      implicit none
      double precision, intent(in):: input(Nr)
      double precision:: aux
      integer:: i

      aux = 0.0d0
      do i=1, Nr-1
         aux = aux + ( input(i) + input(i+1) )*0.5d0*drcoll(i)
      enddo
   end function

   !-----------------------------------------------------------------------------
   !> Performs the integration of the 1d complex array @par input in the radial direction.
   pure function c_integrate_r(input) result(aux)
      implicit none
      double complex, intent(in):: input(Nr)
      double complex:: aux
      integer:: i

      aux = 0.0d0
      do i=1, Nr-1
         aux = aux + ( input(i) + input(i+1) )*0.5d0*drcoll(i)
      enddo
   end function

   !-----------------------------------------------------------------------------
   !> Computes the magnetic dissipation truncated up to degree @par mmax.
   subroutine save_magnetic_dissipation(mmax)
      implicit none
      integer, intent(in):: mmax
      double complex:: psi, psi1, psi2, phi, phi1, phi2, phi3
      double complex:: psicc, phicc, phi1cc
      double complex:: diss_tmp, dissip(Nr)
      complex(kind=8), parameter:: ci=(0.0d0,1.0d0)
      double precision:: Btt,Btt0,Bpp,Bpp0
      double precision:: p3(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      integer:: i, j, l, m
      double precision:: r


      ! B-dissipation, tor-tor, m<>0
      do i=1, Nr
         dissip(i) = 0.0d0
         r = rcoll(i)
         do l=1,Nt_s
            do m=m0,min(l,mmax),m0
               j = 2*(m/m0) + 1
               j = j-blk_ps_start(mpi_rank) - 1
               if((j.ge.1) .and. (j.lt.blk_ps_size(mpi_rank))) then
                  psi   = field_tor(l,j,i)
                  psi1  = field_tor_dr(l,j,i)
                  psi2  = field_tor_ddr(l,j,i)
                  psicc = field_tor(l,j,i)
                  if(m.lt.mmax) then
                     psi   = psi   + ci*field_tor( l,j+1,i)
                     psi1  = psi1  + ci*field_tor_dr(l,j+1,i)
                     psi2  = psi2  + ci*field_tor_ddr(l,j+1,i)
                     psicc = psicc - ci*field_tor( l,j+1,i)
                  endif
                  diss_tmp = 2*pi*llp1(l)*plmfac(l,m)**2*psicc*(psi2-llp1(l)*psi/r**2)
                  call drs_apply_hypDiff(diss_tmp, l)
                  dissip(i)=dissip(i)+diss_tmp
               endif
            enddo
         enddo
      enddo
      Btt = dble(c_integrate_r(dissip)*2)

      ! B-dissipation, tor-tor, m=0
      if(blk_ps_start(mpi_rank).eq.0) then
         do i=1,Nr
            dissip(i)=0
            do l=1,Nt_s
               diss_tmp = 4*pi*llp1(l)/(2*l+1)*field_tor(l,1,i)*field_tor_lap(l,1,i)
               dissip(i) = dissip(i) + diss_tmp
            enddo
         enddo
         Btt0 = dble(c_integrate_r(dissip)*2)
      endif

      do l=0,Nt_s
         do j=1,blk_ps_size(mpi_rank)
            p3(l,j,1:Nr) = radial_derivative_r2r(field_pol_ddr(l,j,1:Nr))
         enddo
      enddo

      ! B-dissipation, pol-pol, m<>0
      !! \todo Is giving me NaN's. Check why!
      do i=1,Nr
         dissip(i)=0
         r=rcoll(i)
         do l=1,Nt_s
            do m=m0,min(l,mmax),m0
               j=2*(m/m0)+1-1
               j=j-blk_ps_start(mpi_rank)
               if((j.ge.1) .and. (j.lt.blk_ps_size(mpi_rank))) then
                  phi=field_pol(l,j,i)
                  phi1=field_pol_dr(l,j,i)
                  phi2=field_pol_ddr(l,j,i)
                  phi3=p3(l,j,i)
                  phicc=field_pol(l,j,i)
                  phi1cc=field_pol_dr(l,j,i)
                  if(m.lt.mmax) then
                     phi=field_pol(l,j,i)+ci*field_pol(l,j+1,i)
                     phi1=field_pol_dr(l,j,i)+ci*field_pol_dr(l,j+1,i)
                     phi2=field_pol_ddr(l,j,i)+ci*field_pol_ddr(l,j+1,i)
                     phi3=p3(l,j,i)+ci*p3(l,j+1,i)
                     phicc=field_pol(l,j,i)-ci*field_pol(l,j+1,i)
                     phi1cc=field_pol_dr(l,j,i)-ci*field_pol_dr(l,j+1,i)
                  endif
                  diss_tmp=2*pi*llp1(l)*plmfac(l,m)**2*(llp1(l)/r**2*phicc*(phi2-llp1(l)/r**2*phi)&
                           +phi1cc*(phi3-llp1(l)/r**2*phi1+2*llp1(l)/r**3*phi))
                  call drs_apply_hypDiff(diss_tmp, l)
                  dissip(i)=dissip(i)+diss_tmp
               endif
            enddo
         enddo
      enddo
      Bpp = dble( c_integrate_r(dissip)*2 )

      ! B-dissipation, pol-pol, m=0
      if(blk_ps_start(mpi_rank).eq.0) then
         do i=1,Nr
            dissip(i)=0
            r=rcoll(i)
            do l=1,Nt_s
               phi=field_pol(l,1,i)
               phi1=field_pol_dr(l,1,i)
               phi2=field_pol_ddr(l,1,i)
               phi3=p3(l,1,i)
               phicc=field_pol(l,1,i)
               phi1cc=field_pol_dr(l,1,i)
               diss_tmp=4*pi*llp1(l)/(2*l+1)*&
      (llp1(l)/r**2*phicc*(phi2-llp1(l)/r**2*phi)+phi1cc*(phi3-&
      llp1(l)/r**2*phi1+2*llp1(l)/r**3*phi))
                call drs_apply_hypDiff(diss_tmp, l)
               dissip(i)=dissip(i)+diss_tmp
            enddo
         enddo
         Bpp0 = dble(c_integrate_r(dissip))
      endif


      call sum_over_all_cpus(Btt)
      call sum_over_all_cpus(Btt0)
      call sum_over_all_cpus(Bpp)
      call sum_over_all_cpus(Bpp0)

      if(mpi_rank.eq.0) then
         write(unit_dissB,'(6(D17.9,X))') time, Btt, Btt0, Bpp, Bpp0, Btt+Btt0+Bpp+Bpp0
      endif

   end subroutine save_magnetic_dissipation

   !-----------------------------------------------------------------------------
   ! A.T. & E.G. (Sep 1997)
   ! E.G. (Aug 1999): Bug Fix for hyper-diffusivity
   !> Computes the viscous dissipation.
   !! \todo Separate computing from writing.
   subroutine save_flow_dissipation(mmax)
      implicit none
      integer, intent(in):: mmax
      complex(kind=8), parameter:: ci=(0.0d0,1.0d0)
      double complex:: psi,psi1,psi2,phi,phi1,phi2,phi3
      double complex:: diss_tmp,dissip(Nr)
      double complex:: psicc,phicc,phi1cc
      double precision:: vtt,vtt0,vpp,vpp0, vttl,vtt0l,vppl,vpp0l
      double precision:: p3(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      integer:: i, j, l, m
      double precision:: r, r2

      vttl  = 0.0
      vtt0l = 0.0
      vppl  = 0.0
      vpp0l = 0.0

      ! v-dissipation, tor-tor, m<>0
      do i=1,Nr
         dissip(i) = 0
         r  = rcoll(i)
         r2 = rcoll2(i)
         do l=1,Nt_s
            do m=m0,min(l,mmax),m0
               j=2*(m/m0)+1-1
               j=j-blk_ps_start(mpi_rank)
               if((j.ge.1) .and. (j.lt.blk_ps_size(mpi_rank))) then
                   if(m.lt.mmax) then
                      psi   = r2*( flow_tor(l,j,i)     + ci*flow_tor(l,j+1,i) )
                      psicc = r2*( flow_tor(l,j,i)     - ci*flow_tor(l,j+1,i) )
                      psi1  = r2*( flow_tor_dr(l,j,i)  + ci*flow_tor_dr(l,j+1,i) ) &
                            + 2*r*( flow_tor(l,j,i)    + ci*flow_tor(l,j+1,i) )
                      psi2  = r2*( flow_tor_ddr(l,j,i) + ci*flow_tor_ddr(l,j+1,i) ) &
                            + 4*r*( flow_tor_dr(l,j,i) + ci*flow_tor_dr(l,j+1,i) ) &
                            + 2*( flow_tor(l,j,i) + ci*flow_tor(l,j+1,i) )
                   else
                      psi   = r2*flow_tor(l,j,i)
                      psicc = psi
                      psi1  = r2*flow_tor_dr(l,j,i)  + 2*r*flow_tor(l,j,i)
                      psi2  = r2*flow_tor_ddr(l,j,i) + 4*r*flow_tor_dr(l,j,i) + 2*flow_tor(l,j,i)
                   endif
                   diss_tmp = 2*pi*llp1(l)*plmfac(l,m)**2*psicc*( psi2 - llp1(l)*psi/r2 )
                   call drs_apply_hypDiff(diss_tmp, l)
                   dissip(i) = dissip(i) + diss_tmp
               endif
            enddo
         enddo
      enddo
      vttl=0.0
      do i=1,Nr-1
         vttl = vttl + dble(dissip(i) + dissip(i+1))*0.5d0*drcoll(i)
      enddo
      vttl = vttl*2

      ! v-dissipation, tor-tor, m=0
      if(blk_ps_start(mpi_rank).eq.0) then
         do i=1,Nr
            dissip(i)=0
            r=rcoll(i)
            r2=rcoll2(i)
            do l=1,Nt_s
               psi   = r2*flow_tor(l,1,i)
               psi1  = r2*flow_tor_dr(l,1,i)  + 2*r*flow_tor(l,1,i)
               psi2  = r2*flow_tor_ddr(l,1,i) + 4*r*flow_tor_dr(l,1,i) + 2*flow_tor(l,1,i)
               psicc = psi
               diss_tmp = 4*pi*llp1(l)/(2*l+1)*psicc*( psi2 - llp1(l)/r2*psi )
               call drs_apply_hypDiff(diss_tmp, l)
               dissip(i)=dissip(i)+diss_tmp
            enddo
         enddo
         vtt0l=0.0
         do i=1,Nr-1
            vtt0l = vtt0l + 0.5d0*dble(dissip(i) + dissip(i+1))*drcoll(i)
         enddo
      endif


      do l=0,Nt_s
         do j=1,blk_ps_size(mpi_rank)
            p3(l,j,1:Nr) = radial_derivative_r2r(flow_pol_ddr(l,j,1:Nr))
         enddo
      enddo

      ! v-dissipation, pol-pol, m<>0
      do i=1,Nr
         dissip(i) = 0
         r  = rcoll(i)
         r2 = rcoll2(i)
         do l=1,Nt_s
            do m=m0,min(l,mmax),m0
               j = 2*(m/m0) + 1 - 1
               j = j - blk_ps_start(mpi_rank)
               if((j.ge.1) .and. (j.lt.blk_ps_size(mpi_rank))) then
                  if(m.lt.mmax) then
                     phi    = r*( flow_pol(l,j,i)     + ci*flow_pol(l,j+1,i) )
                     phi1   = r*( flow_pol_dr(l,j,i)  + ci*flow_pol_dr(l,j+1,i) ) &
                            +   (flow_pol(l,j,i)      + ci*flow_pol(l,j+1,i))
                     phi2   = r*( flow_pol_ddr(l,j,i) + ci*flow_pol_ddr(l,j+1,i)) &
                            + 2*(flow_pol_dr(l,j,i)   + ci*flow_pol_dr(l,j+1,i))
                     phi3   = r*( p3(l,j,i)           + ci*p3(l,j+1,i))           &
                            + 3*(flow_pol_ddr(l,j,i)  + ci*flow_pol_ddr(l,j+1,i))
                     phicc  = r*( flow_pol(l,j,i)     - ci*flow_pol(l,j+1,i))
                     phi1cc = r*( flow_pol_dr(l,j,i)  - ci*flow_pol_dr(l,j+1,i))  &
                            + (flow_pol(l,j,i)        - ci*flow_pol(l,j+1,i))
                  else
                     phi    = r*flow_pol(l,j,i)
                     phi1   = r*flow_pol_dr(l,j,i) + flow_pol(l,j,i)
                     phi2   = r*flow_pol_ddr(l,j,i) + 2*flow_pol_dr(l,j,i)
                     phi3   = r*p3(l,j,i) + 3*flow_pol_ddr(l,j,i)
                     phicc  = r*flow_pol(l,j,i)
                     phi1cc = r*flow_pol_dr(l,j,i) + flow_pol(l,j,i)
                  endif
                  diss_tmp = 2*pi*llp1(l)*plmfac(l,m)**2* &
                            ( llp1(l)/r2*phicc*(phi2-llp1(l)/r2*phi) + &
                                          phi1cc*(phi3-llp1(l)/r2*phi1 + &
                                          2*llp1(l)/r**3*phi) )
                  call drs_apply_hypDiff(diss_tmp, l)
                  dissip(i) = dissip(i) + diss_tmp
               endif
            enddo
         enddo
      enddo
      vppl = 0.0
      do i=1,Nr-1
         vppl=vppl+dble(dissip(i)+dissip(i+1))*0.5*drcoll(i)
      enddo
      vppl = vppl*2

      ! v-dissipation, pol-pol, m=0
      if(blk_ps_start(mpi_rank).eq.0) then
         do i=1,Nr
            dissip(i) = 0
            r  = rcoll(i)
            r2 = rcoll2(i)
            do l=1,Nt_s
               phi    = r*flow_pol(l,1,i)
               phi1   = flow_pol(l,1,i) + r*flow_pol_dr(l,1,i)
               phi2   = 2*flow_pol_dr(l,1,i)  + r*flow_pol_ddr(l,1,i)
               phi3   = 3*flow_pol_ddr(l,1,i) + r*p3(l,1,i)
               phicc  = r*flow_pol(l,1,i)
               phi1cc = flow_pol(l,1,i) + r*flow_pol_dr(l,1,i)
               diss_tmp = 4*pi*llp1(l)/(2*l+1)*(llp1(l)/r2*phicc*(phi2-llp1(l)/r2*phi) + &
                                                          phi1cc*(phi3-llp1(l)/r2*phi1 + &
                                                          2*llp1(l)/r**3*phi))
               call drs_apply_hypDiff(diss_tmp, l)
               dissip(i) = dissip(i) + diss_tmp
            enddo
         enddo
         vpp0l=0
         do i=1,Nr-1
            vpp0l=vpp0l+dble(dissip(i)+dissip(i+1))*0.5*drcoll(i)
         enddo
      endif

      Vtt  = vttl
      vtt0 = vtt0l
      vpp  = vppl
      vpp0 = vpp0l
      call sum_over_all_cpus(vtt)
      call sum_over_all_cpus(vtt0)
      call sum_over_all_cpus(vpp)
      call sum_over_all_cpus(vpp0)
      if(mpi_rank.eq.0) then
         write(unit_dissu,'(6D20.11)') time,vtt,vtt0,vpp,vpp0,vtt+vtt0+vpp+vpp0
      endif
   end subroutine save_flow_dissipation

   !-----------------------------------------------------------------------------
   !> Saves some flow coefficients at the present instant.
   !! \todo Should take a list of l's and m's and reply with a list of values
   !! \todo Writing should be moved to io.
   subroutine save_flow_coeffs()
      implicit none
      integer:: i, m, j, jg
      double precision::wert(5)
      integer:: rank(5)
      i = Nr/2
      wert(1) = flow_tor(1,1,i)                      !  wert1 = tor(l=1,m=0)
      rank(1) = 0

      m  = m0
      jg = 2*(m/m0)
      j  = jg - blk_ps_start(mpi_rank) + 1
      if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
         !  wert(2) = re_pol(l=m0,m=m0)
         !  wert(3) = im_pol(l=m0,m=m0)
         wert(2) = flow_pol(m0,j,i)
         wert(3) = flow_pol(m0,j+1,i)
         rank(2:3) = mpi_rank
      endif

      m  = 2*m0
      jg = 2*(m/m0)
      j  = jg - blk_ps_start(mpi_rank) + 1
      if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
         !  wert4 = re_pol(l=2*m0,m=2*m0)
         !  wert5 = re_tor(l=2*m0+1,m=2*m0)
         wert(4) = flow_pol(2*m0,j,i)
         wert(5) = flow_tor(2*m0+1,j,i)
         rank(4:5) = mpi_rank
      endif

      call drs_gather_vars(rank,wert)
      if(mpi_rank.eq.0) then
         !-- file.koeu:
         write(unit_koeu,'(D17.9,5D13.5)') time, (wert(i), i=1,5)
      endif  !  root
   end subroutine save_flow_coeffs

   !-----------------------------------------------------------------------------
   !> Saves some field coefficients at the present instant.
   !! \todo Should take a list of l's and m's and reply with a list of values
   !! \todo Writing should be moved to io.
   subroutine save_field_coeffs()
      implicit none
      integer:: i, m, j, jg
      double precision:: wert(10)
      integer:: rank(10)
      i = Nr/2
      wert(1)  = field_tor(1,1,i)                 !  wert1 = tor(l=1,m=0)
      rank(1) = 0

      m  = m0
      jg = 2*(m/m0)
      j  = jg - blk_ps_start(mpi_rank)+1
      if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
         wert(2) = field_pol(m0,   j,  i)          !  wert2 = re_pol(l=m0,m=m0)
         wert(3) = field_pol(m0,   j+1,i)          !  wert3 = im_pol(l=m0,m=m0)
         wert(4) = field_pol(m0+1, j,  i)          !  wert4 = re_pol(l=m0+1,m=m0)
         rank(2:4) = mpi_rank

         wert(7) = field_tor(m0,   j,  i)          !  wert7 = re_tor(l=m0,m=m0)
         wert(8) = field_tor(m0+1, j,  i)          !  wert8 = re_tor(l=m0+1,m=m0)
         rank(7:8) = mpi_rank
      endif

      m  = 2*m0
      jg = 2*(m/m0)
      j  = jg - blk_ps_start(mpi_rank)+1
      if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
         wert(5) = field_pol(2*m0,j,i)            !  wert5 = re_pol(l=2*m0,m=2*m0)
         wert(6) = field_pol(2*m0+1,j,i)          !  wert6 = re_pol(l=2*m0+1,m=2*m0)
         rank(5:6) = mpi_rank
      endif

      jg = 1
      j  = jg - blk_ps_start(mpi_rank)+1
      if(j.ge.1 .and. j.le.blk_ps_size(mpi_rank)) then
         wert(9)  = field_pol(1,j,i)          !  wert9 = re_pol(l=1,m=0)  (dipole)
         wert(10) = field_pol(2,j,i)          !  wert10 = re_pol(l=2,m=0) (quadr.)
         rank(9:10) = mpi_rank
      endif

      call drs_gather_vars(rank,wert)
      if(mpi_rank.eq.0) then
         write(unit_koeb,'(D17.9,X,10(D17.9,X))') time, (wert(i), i=1, 10)
      endif  !  root
   end subroutine save_field_coeffs

   !-----------------------------------------------------------------------------
   !> Checks resolution of Hartmann layers.
   subroutine check_resolution_Hartman(Rm,error)
      implicit none
      double precision, intent(in):: Rm
      integer, intent(out):: error
      double precision:: delta
      double precision:: dr, dt, dp
      integer, save:: counter = 0

      error = 0
      delta = 1./sqrt(Rm)
      dr = 0.5/(Nr-1)
      dt = 0.5*pi/Nt
      dp = pi/(m0*(Np-1))

      if( delta.lt. dr ) then
         spew 'Radial resolution must be at least 2 times better'
         spew 'to resolve Hartmann layers!'
         error = WARN_RES_TOO_LOW_R
      endif
      if( delta.lt.dt ) then
         spew 'Meridional resolution must be at least 2 times better'
         spew 'to resolve Hartmann layers!'
         error = WARN_RES_TOO_LOW_TH
      endif
      if( delta.lt.dp ) then
         spew 'Azimuthal resolution must be at least 2 times better'
         spew 'to resolve Hartmann layers!'
         error = WARN_RES_TOO_LOW_PH
      endif
      if (error.ne.0) then
         counter = counter + 1
      else
         counter = 0
      endif
      if (counter.gt.5) then
         spew "I hit this problem too many times. Aborting!"
         error = ERR_RES_TOO_LOW_FOR_HARTMAN
      endif
   end subroutine check_resolution_Hartman

   !-----------------------------------------------------------------------------
   !> Kinetic energy of the current flow in real space.
   double precision function KineticEnergy()
      implicit none
      KineticEnergy = energy(flow_r_t, flow_t_t, flow_p_t)
   end function

   !-----------------------------------------------------------------------------
   !> Magnetic energy of the current field in real space.
   double precision function MagneticEnergy()
      implicit none
      MagneticEnergy = energy(field_r_t, field_t_t, field_p_t)
   end function

   !----------------------------------------------------------------
   !> Saves the nusselt number to file with unit \a unit_nu
   !!
   !! This corresponds to the files with extension nu
   !! The file structure is:\n
   !! time, Nu_ICB, Nu_CMB
   subroutine save_nusselt_number()
      implicit none
      double precision:: Nu_CMB, Nu_ICB
      Nu_ICB = nusselt(Nr)
      Nu_CMB = nusselt(1)
      write (unit_nu,'(3(D17.9,X))') time, Nu_ICB, Nu_CMB       !  file.nu
   end subroutine

   !-----------------------------------------------------------------------------
   !> measures and saves quantities of interest (in (theta,phi,r)).
   !! All inputs have to be in (theta,phi,r) space
   !! common /energies/ are calculates in measure_lm.
   !!
   subroutine measure(temp2_t, rotu_r_t, rotu_theta_t, rotu_phi_t)
      implicit none
      !-------------------------------------------------------------
      !-- definition of energies:
      !-- a) program variables
      !--    ending with s: (l+m) even
      !--    ending with a: (l+m) odd.
      !-- b) convection:
      !--    e-symmetric: pol-s,tor-a: mckes,drkea,nkes-etors,etora.
      !--    e-asymm.   : pol-a,tor-s: mckea,drkes,nkea-etora,etors.
      !-- c) magnetic field:
      !--    dipol      : pol-a,tor-s: Bmckea,Bdrkes,Bnkea-Betora,Betors.
      !--    quadrupol  : pol-s,tor-a: Bmckes,Bdrkea,Bnkes-Betors,Betora.
      !-------------------------------------------------------------

      double precision, intent(in):: temp2_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: rotu_r_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: rotu_theta_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: rotu_phi_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
!-- local workarray:
      double precision:: advect(Nr)

      double precision:: helicity_north, helicity_south
      double precision:: enstr, am_x, am_y, am_z
      double precision:: fstr
      double precision:: te2(Nr), ur2(Nr)
      !> energies from measure_lm:
      double precision:: nkes,nkea,etors,etora,drkes,drkea,mckes,mckea
      double precision:: Bnkes,Bnkea,Betors,Betora,Bdrkes,Bdrkea,Bmckes,Bmckea

      double precision, save:: Ekin1 = 0.0d0, Ekin2 = 0.0d0, EB1 = 0.0d0, EB2 = 0.0d0
      integer:: Ntl

      !-- compute the total kinetic energy
      Ekin = KineticEnergy()
      if (Ekin.ne.Ekin) call drs_abort(ERR_KENERGY_IS_NAN)

      ! Growth rate
      if(mpi_rank.eq.0 .and. (drs_calc_type.eq.LinearThermalOnset.or. &
                              drs_calc_type.eq.LinearCompositionalOnset)) then
         !-- grothrate:
         !-- f(flow_tor) = exp(lt)(+const)
         !-- f'(flow_tor)= l*exp(lt)
         if(Ekin1.le.0.D0) Ekin1=0.2*Ekin
         if(Ekin2.le.0.D0) Ekin2=0.1*Ekin
         fstr  = (Ekin-Ekin1)/time_since_last_sample
         groth = 2*fstr/(Ekin+Ekin1)
         Ekin2 = Ekin1
         Ekin1 = Ekin
      endif

      !-- compute the total magnetic energy
      if(field_present.and.field_evolves) then
         ! magnetic Reynolds number Rm = Pm * sqrt(2Ekin/V)
         Rm = Pm*sqrt(2.*Ekin/(4./3.*pi*(1.-eta**3)/(1.-eta)**3))

         EB  = MagneticEnergy()

         if(mpi_rank.eq.0 .and. drs_calc_type.eq.KinematicDynamo) then
            !-- grothrate:
            !-- f(flow_tor) = exp(lt)(+const)
            !-- f'(flow_tor)= l*exp(lt)
            if(EB1.le.0D0) EB1 = 2*EB
            if(EB2.le.0D0) EB2 = 3*EB
            fstr  = (EB1-EB )/time_since_last_sample
            groth = fstr/(EB+EB1)*2.d0
            EB2 = EB1
            EB1 = EB
         endif
         call compute_BE_components(Bnkes, Bnkea, Betors, Betora, Bdrkes, Bdrkea, Bmckes, Bmckea)
      endif  !  if(field_evolves)

      if(flow_evolves) then
         call compute_angular_momentum(flow_t_t, flow_p_t,am_x, am_y, am_z)
         call save_angular_momentum(am_x, am_y, am_z)

         ! compute the enstrophy Int( rot_u**2 )
         enstr = 2.0d0*energy(rotu_r_t, rotu_theta_t, rotu_phi_t)


         call integrate_power_surf(flow_r_t, 2, ur2)
         call integrate_power_surf(flow_t_t, 2, ut2)
         call integrate_power_surf(flow_p_t, 2, up2)
         ur_avg = ur_avg + ur2*time_since_last_sample
         ut_avg = ut_avg + ut2*time_since_last_sample
         up_avg = up_avg + up2*time_since_last_sample

         call compute_KE_components(nkes, nkea, etors, etora, drkes, drkea, mckes, mckea)
         call compute_helicities( flow_r_t, flow_t_t, flow_p_t, &
                                rotu_r_t, rotu_theta_t, rotu_phi_t, &
                                helicity_north, helicity_south)
         if(mpi_rank==0) then
	    !-- save kinetic energies:
            !-- time,Ekin,meanp_sy,meant_sy,flp_sy,flt_sy,
            !--           meanp_as,meant_as,flp_as,flt_as,hel,hel2
            write (unit_ek,'(12(D17.9,X))')&
                   time,Ekin,mckes,drkea,nkes-etors,etora,&     ! e-symm.
                             mckea,drkes,nkea-etora,etors,&      ! e-asymm.
                             helicity_north, helicity_south
	 endif

      endif

      if(temp_evolves) then
         call integrate_power_surf(temp2_t, 2, te2)
         t2_avg = t2_avg + te2*time_since_last_sample

         call compute_advection(flow_r_t, temp2_t, advect)
         adv_avg = adv_avg + advect*time_since_last_sample

         !-----
         !-- compute the Nusselt numbers at both boundaries and write to file
         !-----
         if(mpi_rank.eq.0) then        !  only l=m=0.
            call save_nusselt_number()
         endif                   !  mpi_rank.eq.0
      endif

      !-- some output on stdout:
      if(mpi_rank.eq.0) then
         select case (drs_calc_type)
            case (LinearThermalOnset, LinearCompositionalOnset)
               write(* ,'(A26,X,D17.9,X,I7,X,2(D17.9,X))') 'time,step,Ekin,groth',time,steps,Ekin,groth
            case (KinematicDynamo,KinematicDynamo+Compositional)
               write(* ,'(A26,X,D17.9,X,I7,X,4(D17.9,X))') 'time,step,Ekin,EB,Rm,groth',time,steps,Ekin,EB,Rm,groth
            case (NonlinearConvection,NonlinearConvection+Compositional)
               write(* ,'(A26,X,D17.9,X,I7,X,2(D17.9,X))') 'time,step,Ekin',time,steps,Ekin
            case (NonlinearDynamo,MagnetoConvection,NonlinearDynamo+Compositional,MagnetoConvection+Compositional)
#ifdef STAR
               write(* ,'(A26,X,D17.9,X,I7,X,2(D17.9,X))') 'time,step,Rm,Ra_t',time,steps,Rm,Ra_t
#else
               write(* ,'(A26,X,D17.9,X,I7,X,3(D17.9,X))') 'time,step,Ekin,EB,Rm',time,steps,Ekin,EB,Rm
#endif
         endselect


         Ntl = blk_t_size(mpi_rank)-1
         if(flow_evolves) then
            !-- save a few ur's to determine drift etc.
            write (unit_u_mid,'(5(D17.9,X))') time,flow_r_t(Ntl/2,1,Nr/2), flow_r_t(Ntl/2,1,Nr/2),&
                                                   flow_r_t(Ntl/2,1,Nr/2), flow_r_t(Ntl*5/8,1,Nr/2)

         endif  ! if(flow_evolves)

         if(field_evolves) then
            write (unit_eb,'(10(D17.9,X))')&
               time,EB,Bmckea,Bdrkes,Bnkea-Betora,Betors,&  ! dipol
                       Bmckes,Bdrkea,Bnkes-Betors,Betora   ! quadr.
         endif  ! if(field_evolves)
      endif  ! if(mpi_rank.eq.0)
   end subroutine measure

   !-----------------------------------------------------------------------------
   !> Saves the three components of the angular momentum.
   subroutine save_angular_momentum(am_x, am_y, am_z)
      implicit none
      double precision, intent(in):: am_x, am_y, am_z
      !-- save angular momentum
      if(mpi_rank.eq.0) Write (unit_am,'(4(D17.9,X))') time, am_x, am_y, am_z
   end subroutine

   !-----------------------------------------------------------------------------
   !> Computes the north and south hemisphere helicities
   subroutine compute_helicities(ur, ut, up, rotu_r, rotu_t, rotu_p, helicity_south, helicity_north)
      implicit none
      double precision, intent(in):: ur(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: ut(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: up(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: rotu_r(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: rotu_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: rotu_p(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(out):: helicity_south, helicity_north
      integer:: j, l, lg, Ntl, Nts
      Ntl = blk_t_size(mpi_rank)-1
      Nts = blk_t_start(mpi_rank)

      ! compute the helicity_north Int( u*rot_u ) in the northern hemisphere (Ntl assumed even)
      helicity_north = 0.d0
      do j=1,Np
         do l=0,Ntl
            lg=l+Nts
            if(lg.le.Nt/2-1) &
               helicity_north = helicity_north + &
                          sum( (ur(l,j,1:Nr)*rotu_r(l,j,1:Nr) + &
                                ut(l,j,1:Nr)*rotu_t(l,j,1:Nr) + &
                                up(l,j,1:Nr)*rotu_p(l,j,1:Nr) )*dOmega(l,1:Nr) )
         enddo
      enddo
      call sum_over_all_cpus(helicity_north)

      ! compute the helicity_north Int( u*rot_u ) in the southern hemisphere (Ntl assumed even)
      helicity_south = 0.d0
      do j=1,Np
         do l=0,Ntl
            lg=l+Nts
            if(lg.ge.Nt/2+1) &
               helicity_south = helicity_south + &
                          sum( (ur(l,j,1:Nr)*rotu_r(l,j,1:Nr) + &
                                ut(l,j,1:Nr)*rotu_t(l,j,1:Nr) + &
                                up(l,j,1:Nr)*rotu_p(l,j,1:Nr) )*dOmega(l,1:Nr) )
         enddo
      enddo
      call sum_over_all_cpus(helicity_south)
   end subroutine

   !-----------------------------------------------------------------------------
   !> Computes the heat transported by advection.
   !! as \f[ Q(r) = \int\int
   !! u_r(r,\theta,\phi)*(\Theta(r,\theta,\phi)+T_S(r)) \sin\theta d\theta d\phi
   !! \f]
   subroutine compute_advection(ur, temp, advect)
      implicit none
      double precision, intent(in):: ur(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: temp(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(out):: advect(Nr)
      double precision:: aux(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      integer:: i, j, l, Ntl
      Ntl = blk_t_size(mpi_rank)-1

      do i=1, Nr
         do j=1, Np
            do l=0, Ntl
               aux(l,j,i) = ur(l,j,i)*( temp(l,j,i) + temp_profile(i) )
            enddo
         enddo
      enddo
      call integrate_power_surf(aux, 1, advect)
   end subroutine

   !-----------------------------------------------------------------------------
   !> Computes the Nusselt number, that is, the the ratio between the convective
   !! and the diffusive heat fluxes.
   !!
   !! \f[ Nu = \frac{\partial (T + \Theta)/ \partial r}{\partial T / \partial r} \f]
   !! \todo It only works for serial runs. Needs to be parallelized.
   pure function nusselt(r)
      implicit none
      double precision:: nusselt
      integer, intent(in):: r !< radius

      nusselt = temp_dr(0,1,r)/temp_profile_dr(r) + 1
   end function nusselt

   !-----------------------------------------------------------------------------
   !> Computes and saves the three cartesian components of the total angular momentum
   subroutine compute_angular_momentum(u_t, u_p, ang_mom_x, ang_mom_y, ang_mom_z)
      implicit none
      !> Theta component of the flow in real (tpr) space.
      double precision, intent(in):: u_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      !> Phi component of the flow in real (tpr) space.
      double precision, intent(in):: u_p(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      !> Cartesian components of the total angular momentum.
      double precision, intent(out):: ang_mom_x, ang_mom_y, ang_mom_z
      double precision:: phi, cosp(Np), sinp(Np)
      integer:: j, l, lg
      integer:: Ntl, Nts

      Ntl = blk_t_size(mpi_rank)-1
      Nts = blk_t_start(mpi_rank)

      do j=1, Np
         phi = 2*pi*(j-1)/dble(m0*Np)
         cosp(j) = dcos(phi)
         sinp(j) = dsin(phi)
      enddo
      ! Lz = Int( r*Uphi*sin(theta) )
      ang_mom_z=0.d0
      do j=1,Np
         do l=0,Ntl
            lg  = l + Nts
            ang_mom_z = ang_mom_z + sum( rcoll(1:Nr)*sintheta(lg)*u_p(l,j,1:Nr)*dOmega(l,1:Nr) )
         enddo
      enddo
      call sum_over_all_cpus(ang_mom_z)

      ! Lx = Int( -r*(Uphi*cos(phi)*cos(theta)+Utheta*sin(phi)) )
      ang_mom_x=0.d0
      do j=1, Np
         do l=0, Ntl
            lg = l + Nts
            ang_mom_x = ang_mom_x - sum( (costheta(lg)*cosp(j)*u_p(l,j,1:Nr) + &
                                                       sinp(j)*u_t(l,j,1:Nr))*dOmega(l,1:Nr) )
         enddo
      enddo
      call sum_over_all_cpus(ang_mom_x)

      ! Ly = Int( r*(Utheta*cos(phi)-Uphi*sin(phi)*cos(theta)) )
      ang_mom_y=0
      do j=1,Np
         do l=0,Ntl
            lg = l + Nts
            ang_mom_y = ang_mom_y + sum(            (cosp(j)*u_t(l,j,1:Nr) - &
                                        costheta(lg)*sinp(j)*u_p(l,j,1:Nr) )*dOmega(l,1:Nr) )
         enddo
      enddo
      call sum_over_all_cpus(ang_mom_y)
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
