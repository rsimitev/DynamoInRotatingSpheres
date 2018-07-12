! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2015
!> Takes care of evolving the momentum equation.
#include "drsDefs.F90"
module drs_momentum_equation
   use drs_mpi
   use drs_params
   use drs_dims
   use drs_time
   use drs_legendre, only:llp1, gauleg
   use drs_hypDiff
   use drs_radial
   use drs_flow
   use drs_field
   use drs_temp
#ifdef COMP
   use drs_comp
#endif
   implicit none
   !> The Crank-Nicholson inverse operator for the toroidal field.
   double precision, allocatable:: inv_lhs_NS_tor(:,:,:)
   !> The Crank-Nicholson inverse operator for the poloidal field.
   double precision, allocatable:: inv_lhs_NS_pol(:,:,:)
   !> NS for Navier-Stokes
   double precision, allocatable:: rhs_NS_tor(:,:,:)
   double precision, allocatable:: rhs_NS_pol(:,:,:)
   double precision, allocatable:: rhs_NS_tor_old(:,:,:)
   double precision, allocatable:: rhs_NS_pol_old(:,:,:)
   double precision, allocatable:: green1(:,:), green2(:,:)
   double precision, allocatable:: greenD1(:,:), greenD2(:,:)
   double precision, allocatable:: greenS1(:,:), greenS2(:,:)
   double precision, allocatable:: pinv(:,:,:)
   logical, private:: initialised=.false.

contains

   !---------------------------------------------------------------
   !> Allocates memory for the Crank-Nicholson inverse operators
   !! @param error Error code to be forwarded.
   subroutine drs_momentum_equation_init(error)
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
      ! Need to make sure the flow is initialised
      call drs_flow_init(error)
      if(error.gt.0) return

      if(field_present) then
         call drs_field_init(error)
         if(error.gt.0) return
      endif

      if(temp_present) then
         call drs_temp_init(error)
         if(error.gt.0) return
      endif

#ifdef COMP
      if(comp_present) then
         call drs_comp_init(error)
         if(error.gt.0) return
      endif
#endif

      ! Allocate variables needed for the time-stepping
      allocate(rhs_NS_tor_old(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(rhs_NS_pol_old(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(  green1(Nr, 0:Nt_s),  green2(Nr, 0:Nt_s) )
      allocate( greenD1(Nr, 0:Nt_s), greenD2(Nr, 0:Nt_s) )
      allocate( greenS1(Nr, 0:Nt_s), greenS2(Nr, 0:Nt_s) )
      allocate(  inv_lhs_NS_tor(Nr, Nr, 0:Nt_s) )
      allocate(  inv_lhs_NS_pol(Nr, Nr, 0:Nt_s) )
      allocate(            pinv(Nr, Nr, 0:Nt_s) )
      allocate(rhs_NS_tor(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      allocate(rhs_NS_pol(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      rhs_NS_tor_old = 0.0d0
      rhs_NS_pol_old = 0.0d0
      initialised = .true.
      call poissoninv_init( pinv )
      ! Initialise the Crank-Nicholson matrices.
      call cnt_init(inv_lhs_NS_tor, h, 1.0d0,error)
      if(error.gt.0) return
      call cnp_init(inv_lhs_NS_pol, h, 1.0d0,error)
      if(error.gt.0) return
   end subroutine

   !-----------------------------------------------------------------------------
   !> This routine prepares the matrices for the poisson inversion necessary for
   !! solving for the poloidal flow.\n
   !! As opposed to inv_lhs_NS_tor and inv_lhs_TE, pinv will multiply an
   !! array in direct space and return an array in real space! \n
   !! pinv(j,i,l) = ( -l(l+1)/r^2 + d^2/dr^2 + 2/r*d/dr ) delta(j,i)
   !!
   !! Calculates the LU-decomposed radial laplace operator with fixed boundary conditions.
   !!
   !! Uses Numerical Recipes routine ludcmp.
   subroutine poissoninv_init(pinv)
      implicit none

      double precision, intent(out):: pinv(Nr,Nr,0:Nt_s)
      double precision:: inv(Nr,Nr)
      integer:: indx(Nr)
      integer:: i,j,l
      double precision:: dr(Nr), ddr(Nr), aux(Nr)

      pinv = 0.0d0

      !-- pinv(j,i,l) = ( -l(l+1)/r^2 + d^2/dr^2 + 2/r*d/dr ) delta(j,i)
      do i=1, Nr
         ! radarr = delta(j,i)
         ! forall (l=0:Nt_s, j=1:Nr) pinv(j,i,l) = -llp1(l)*radarr(j)/rcoll2(j)
         forall (l=0:Nt_s) pinv(i,i,l) = -llp1(l)/rcoll2(i)
         aux    = 0.d0
         aux(i) = 1.d0
         call radial_dr_ddr_1D_r2r(aux,dr,ddr)
         ! deriv  = d/dr delta(j,i)
         ! radarr = d^2/dr^2 delta(j,i)
         forall(l=0:Nt_s, j=1:Nr)
            pinv(j,i,l) = pinv(j,i,l) + ddr(j) + 2.0d0*dr(j)/rcoll(j)
         endforall
      enddo

      ! This implements boundary conditions: fixed values at boundaries
      ! note: i=1 corresponds to r_o, i=Nr corresponds to r_i
      pinv(1, 1:Nr,0:Nt_s) = 0.d0
      pinv(Nr,1:Nr,0:Nt_s) = 0.d0
      pinv(1, 1, 0:Nt_s) = 1.d0
      pinv(Nr,Nr,0:Nt_s) = 1.d0

      do l=0, Nt_s
         inv = pinv(:,:,l)
         call ludcmp(inv, Nr, Nr, indx)
         !-- invert poisson matrix:
         call matinv(inv, indx, Nr, Nr)
         pinv(:,:,l) = inv
      enddo
   end subroutine poissoninv_init

   !---------------------------------------------------------------
   !> Computes matrices for the Crank-Nicholson scheme on the Toroidal flow.
   !! @par dinv is the resulting matrix.
   !! @par dt is the timestep currently being used
   !! @par fact is a factor for the time step. In the current implementation this is 1.0.
   !! @param error Error code to be forwarded.
   subroutine cnt_init(dinv,dt,fac,error)
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
         return
      endif

      dto2fac = dt/(2.0d0*fac)
      do l=0, Nt_s
         aux => dinv(:,:,l)
         ! dinv contains the matrices describing the effect of Dl
         do j=1,Nr ! mode
            do i=1,Nr ! radius
               aux(i,j) = - poly_ddr(i,j) &
                          - (2-llp1(l))*poly(i,j)/rcoll2(i) &
                          - 4*poly_dr(i,j)/rcoll(i)
               call drs_apply_hypDiff(aux(i,j), l)
               aux(i,j) = aux(i,j)*dto2fac + poly(i,j)
            enddo
         enddo

         ! implements boundary conditions: derivative=0 at boundaries
         ! note: i=1 corresponds to r_o, i=Nr corresponds to r_i
         if(flowBC(1).eq.FreeSlip) then
            ! d/dr u = 0
            do j=1, Nr_s
               aux(Nr,j) = poly_dr(Nr,j)
            enddo
         elseif(flowBC(1).eq.NoSlip) then
            ! u = 0
            aux(Nr,1:Nr) = poly(Nr,1:Nr)
         endif
         if(flowBC(2).eq.FreeSlip) then
            ! d/dr u = 0
            do j=1, Nr_s
               aux(1, j) = poly_dr(1, j)
            enddo
         elseif(flowBC(2).eq.NoSlip) then
            ! u = 0
            aux(1, 1:Nr) = poly(1, 1:Nr)
         endif

         ! Dealiase
         if (Nr_s+1 .le. Nr) aux(:, Nr_s+1:Nr) = 0.d0

         call ludcmp(aux,Nr,Nr,indx)

         !-- invert all cn-matrixes for tor. flow:
         ! Invert each block independently.
         call matinv(aux, indx, Nr, Nr)
      enddo
   end subroutine cnt_init

   !---------------------------------------------------------------
   !> initializes matrix for Crank-Nicholson on the toroidal flow at l=1,m=0
   !! to correct the angular momentum
   !! @param error Error code to be forwarded.
   subroutine cnam(dinv,indx,h,error)
      implicit none
      integer, intent(out):: error
      double precision:: dinv(Nr,Nr),pinvw(Nr,Nr)
      double precision:: h, ho2
      double precision:: weight(Nr),xgau(Nr),x,aux
      integer:: indx(Nr)
      integer:: i,j,l,n
      error=0
      if(.not.initialised) then
         error=1
         return
      endif

      ho2 = 0.5d0*h
      l=1
      do j=1,Nr
         do i=1,Nr
            pinvw(i,j) = ( -poly_ddr(i,j) + (llp1(l)-2)*poly(i,j)/rcoll2(i))*ho2 + poly(i,j)
         enddo
      enddo

      ! implements boundary conditions: derivative=0 at boundaries
      ! note: i=1 corresponds to x=1, i=Nr corresponds to x=-1
      if(flowBC(1).eq.FreeSlip) then
         pinvw(Nr,1:Nr_s) = poly_dr(Nr,1:Nr_s)
      elseif(flowBC(1).eq.NoSlip) then
         pinvw(Nr,1:Nr_s) = 2.0d0*poly(Nr,1:Nr_s)
      endif
      if(flowBC(2).eq.FreeSlip) then
         pinvw(1, 1:Nr_s) = poly_dr(1, 1:Nr_s)
      elseif(flowBC(2).eq.NoSlip) then
         pinvw(1, 1:Nr_s) = 2.0d0*poly(1, 1:Nr_s)
      endif
      ! force angular momentum to zero
      if (l.eq.1) then
         pinvw(1, 1) = 0.d0
         pinvw(Nr,1) = 0.d0
      endif

      call gauleg(-1.D0,1.D0,xgau,weight,Nr)
      do j=1,Nr
         n=j-1
         aux=0.d0
         do i=1, Nr
            x   = xgau(i)
            aux = aux + weight(i)*cos(n*acos(x))*( 0.5*(x+1) + eta/(1-eta) )**4
         enddo
         pinvw(Nr/2,j) = aux
      enddo
      do j=Nr_s+1,Nr
         pinvw(Nr/2,j) = 0.d0
      enddo

      ! Dealiase
      do j=Nr_s+1,Nr
         pinvw(1, j) = 0.d0
         pinvw(Nr,j) = 0.d0
      enddo

      call ludcmp(pinvw,Nr,Nr,indx)

      ! dinv contains the LU decomposed matrix ready for use by Crank-Nicholson
      dinv(1:Nr,1:Nr) = pinvw(1:Nr,1:Nr)
   end subroutine cnam

   !---------------------------------------------------------------
   !> Computes matrices for the Crank-Nicholson scheme on the Poloidal flow.
   !! @par dinv is the resulting matrix.
   !! @par dt is the timestep currently being used
   !! @par fact is a factor for the time step. In the current implementation
   !! this is 1.0.
   !! @param error Error code to be forwarded.
   subroutine cnp_init(dinv, dt, fac, error)
      implicit none
      integer, intent(out):: error
      double precision, target, intent(out):: dinv(Nr,Nr,0:Nt_s)
      double precision, intent(in):: dt, fac
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
      do l=0, Nt_s
         aux => dinv(1:Nr,1:Nr,l)
         do j=1,Nr ! mode
            do i=1,Nr ! radius
               aux(i,j) = - poly_ddr(i,j) - 2*poly_dr(i,j)/rcoll(i) + llp1(l)*poly(i,j)/rcoll2(i)
               call drs_apply_hypDiff(aux(i,j), l)
               aux(i,j) = aux(i,j)*dto2fac + poly(i,j)
            enddo
         enddo
         ! implement boundary conditions: fixed values at boundaries
         ! Rigid boundary conditions are equivalent to settong the
         ! Poloidal flow to zero there.
         aux(1, 1:Nr) = poly(1, 1:Nr)
         aux(Nr,1:Nr) = poly(Nr,1:Nr)

         ! Dealiase
         if (Nr_s+1 .le. Nr) aux(:, Nr_s+1:Nr) = 0.d0

         call ludcmp(aux,Nr,Nr,indx)

         !-- invert all cn-matrixes for pol. flow:
         call matinv(aux, indx, Nr, Nr)
      enddo
   end subroutine cnp_init

   !---------------------------------------------------------------
   !> The non-linear terms of the Navier-Stokes equation
   subroutine NavierStokes(rhs_NS_tor,rhs_NS_pol)
      implicit none
      double precision, intent(out):: rhs_NS_tor(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision, intent(out):: rhs_NS_pol(0:Nt_s, blk_ps_size(mpi_rank), Nr)
      double precision:: force_r(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: force_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: force_p(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: force_aux_r(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: force_aux_t(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: force_aux_p(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: Ra !< The Rayleigh number speciffic to each type of Temperature/Composition profile
      integer:: ps_end, ps_start, t_end, t_start
      integer:: i,l,j,m

      rhs_NS_tor = 0.0d0
      rhs_NS_pol = 0.0d0
      ps_start = blk_ps_start(mpi_rank)
      ps_end   = blk_ps_size(mpi_rank)
      t_start  = blk_t_start(mpi_rank)
      t_end    = blk_t_size(mpi_rank)-1

      ! Only add Reynolds stresses if we are not computing the linear onset
      if((drs_calc_type.ne.LinearThermalOnset.and.&
          drs_calc_type.ne.LinearCompositionalOnset)) then
         force_r = rot_flow_r_t
         force_t = rot_flow_t_t
         force_p = rot_flow_p_t
      else
         force_r = 0.0d0
         force_t = 0.0d0
         force_p = 0.0d0
      endif
      !-- Add the Coriolis force
      forall(j=1:Np, l=0:t_end, i=1:Nr)
         force_r(l,j,i) = force_r(l,j,i) + Ta*costheta(blk_t_start(mpi_rank)+l)
         force_t(l,j,i) = force_t(l,j,i) - Ta*sintheta(blk_t_start(mpi_rank)+l)
      endforall
      !-- Omega' x u = (rotu + tau*z/|z|) x u    (theta,phi,r)
      force_aux_r = force_t*flow_p_t - force_p*flow_t_t
      force_aux_t = force_p*flow_r_t - force_r*flow_p_t
      force_aux_p = force_r*flow_t_t - force_t*flow_r_t
      force_r = force_aux_r
      force_t = force_aux_t
      force_p = force_aux_p

      !-- Lorentz force
      !--  (Omega' x u)-(rotB x B)    (theta,phi,r)
      if(field_present) then
         force_r = force_r - (rot_field_t_t*field_p_t - rot_field_p_t*field_t_t)
         force_t = force_t - (rot_field_p_t*field_r_t - rot_field_r_t*field_p_t)
         force_p = force_p - (rot_field_r_t*field_t_t - rot_field_t_t*field_r_t)
      endif

      ! Prepare to compute poloidal and toroidal terms
      ! Last two arguments are inverted relatively to the intended usage due to the
      ! different definitions of pol/tor for field and flow.
      call vectorField2PolTor_common(force_r,force_t,force_p,rhs_NS_tor,rhs_NS_pol)
      call PolTor_common2PolTor_flow(rhs_NS_pol, rhs_NS_tor)

      if(temp_present) then
         ! Add buoyancy due to temperature anomalies
         if(tempProf.eq.Conduction) then
            Ra = Ra_t*(1.0d0-eta)
         elseif(tempProf.eq.InternalHeating) then
            Ra = Ra_t
         else
            Ra = Ra_t
         endif
         do i=1, Nr
            do j=1,ps_end
               m = blk_ts_start(j)
               forall(l=m:Nt_s, l.ne.0) rhs_NS_pol(l,j,i) = rhs_NS_pol(l,j,i) - Ra*temp(l,j,i)
            enddo
         enddo
      endif
#ifdef COMP
      if (comp_present) then
         ! Add buoyancy due to composition anomalies
         if(compProf.eq.Diffusive) then
            Ra = Ra_c*(1.0d0-eta)
         elseif(compProf.eq.WellMixed) then
            Ra = Ra_c
         else
            Ra = Ra_c
         endif
         do i=1, Nr
            do j=1,ps_end
               m = blk_ts_start(j)
               forall(l=m:Nt_s, l.ne.0) rhs_NS_pol(l,j,i) = rhs_NS_pol(l,j,i) - Ra*comp(l,j,i)
            enddo
         enddo
      !   debug_call save_lmr_quantity(Ra*comp, 'comp_lmr')
      endif
#endif
      !debug_call drs_abort(666)
   end subroutine NavierStokes

   !----------------------------------------------------------------------
   !>
   subroutine just_apply_flow_BC()
      implicit none
      double precision:: flow_tor_new(Nr)
      double precision:: flow_pol_new(Nr)
      integer:: l, j

      jl_do(j,l)
         flow_pol_new(1:Nr) = flow_pol(l, j, 1:Nr)
         call apply_flow_pol_BC(flow_pol_new)
         flow_pol(l,j,1:Nr) = flow_pol_new(1:Nr)
      jl_enddo
      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)

      jl_do(j,l)
         flow_tor_new(1:Nr) = flow_tor(l, j, 1:Nr)
         call apply_flow_tor_BC(flow_tor_new)
         flow_tor(l,j,1:Nr) = flow_tor_new(1:Nr)
      jl_enddo
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)
   end subroutine

   !----------------------------------------------------------------------
   !> Updates the flow using a hybrid Crank-Nicholson/Addams-Bashford implicit
   !! integration scheme.
   !! @param h is the time step;
   !! @param h1 and @param h2 are the Addams-Bashford weights
   !! @param error Error code to be forwarded.
   subroutine update_flow(h, h1, h2,error)
      implicit none
      double precision, intent(in):: h, h1, h2
      integer, intent(out):: error
      double precision:: ho2
      double precision:: pw(Nr)
      double precision:: flow_tor_new(Nr)
      double precision:: flow_pol_new(Nr)

      integer:: l, j, i, blk_size

      error=0
      if(.not.initialised) then
         error=1
         return
      endif
      if(.not.flow_evolves) return

      blk_size = blk_ps_size(mpi_rank)
      ho2 = h/2.0d0
      call update_Green_functions()
      jl_do(j,l)
         pw(1:Nr) = rhs_NS_pol(l, j, 1:Nr)
         rhs_NS_pol(l, j, 1:Nr) = matmul(pinv(1:Nr,1:Nr,l),pw(1:Nr))
      jl_enddo

      jl_do(j,l)
         do i=1, Nr
            flow_pol_new(i) = flow_pol(l,j,i) + flow_pol_lap(l,j,i)*ho2 + rhs_NS_pol(l,j,i)*h1 - rhs_NS_pol_old(l, j, i)*h2
         enddo
         call apply_flow_pol_BC(flow_pol_new)
         ! multiply with inverted CN-matrices:
         ! This leaves all the fields in lmn space
         flow_pol(l,j,1:Nr) = matmul(inv_lhs_NS_pol(1:Nr, 1:Nr, l), flow_pol_new(1:Nr))
         rhs_NS_pol_old(l,j,1:Nr) = rhs_NS_pol(l,j,1:Nr)
      jl_enddo
      if((drs_calc_type.eq.LinearThermalOnset) .or. &
         (drs_calc_type.eq.LinearCompositionalOnset) ) then
         flow_pol(:,1,:) = 0.0d0
      endif
      !-- all fields are given in (lmn) space
      !-- new derivatives are in (lmr).
      ! All fields are dealised in radial_dr_ddr_3D_n2r.
      call radial_dr_ddr_3D_n2r(flow_pol, flow_pol_dr, flow_pol_ddr)
      call applyGreen()

      ! The toroidal flow
      jl_do(j,l)
         do i=1, Nr
            flow_tor_new(i) = flow_tor(l, j, i) + flow_tor_lap(l,j,i)*ho2 + rhs_NS_tor(l, j, i)*h1 - rhs_NS_tor_old(l, j, i)*h2
         enddo
         call apply_flow_tor_BC(flow_tor_new)
         ! multiply with inverted CN-matrices:
         ! This leaves all the fields in lmn space
         flow_tor(l,j,1:Nr) = matmul(inv_lhs_NS_tor(1:Nr, 1:Nr, l), flow_tor_new(1:Nr))
      jl_enddo
      rhs_NS_tor_old = rhs_NS_tor
      if((drs_calc_type.eq.LinearThermalOnset) .or. &
         (drs_calc_type.eq.LinearCompositionalOnset) ) then
         flow_tor(:,1,:) = 0.0d0
      endif
      call radial_dr_ddr_3D_n2r(flow_tor, flow_tor_dr, flow_tor_ddr)
      call update_flow_pol_lap()
      call update_flow_tor_lap()
   end subroutine update_flow

   !----------------------------------------------------------------------
   !> Encapsulates dealing with the Green's functions
   subroutine update_Green_functions()
      !-- get first Green's function
      call mk_green(1.0d0, 0.0d0, green1, greenD1, greenS1)
      !-- get second Green's function
      call mk_green(0.0d0, 1.0d0, green2, greenD2, greenS2)
   end subroutine

   !----------------------------------------------------------------------
   !> Generates the appropriate Green's function for a given boundary value.
   subroutine mk_green(ib, ob, green, greenD, greenS)
      implicit none
      double precision, intent(in):: ib !< Boundary value for the inner boundary.
      double precision, intent(in):: ob !< Boundary value for the outer boundary.
      double precision, intent(out):: green(Nr, 0:Nt_s)
      double precision, intent(out):: greenD(Nr, 0:Nt_s)
      double precision, intent(out):: greenS(Nr, 0:Nt_s)
      double precision:: pw(Nr), y1(Nr), pw_ddr(Nr), pw_dr(Nr)
      integer:: l

      do l=0, Nt_s
         y1(:)  = 0.0d0
         y1(1)  = ib
         y1(Nr) = ob
         pw(1:Nr) = matmul(pinv(1:Nr,1:Nr,l), y1(1:Nr))

         pw(1)  = 0.0d0
         pw(Nr) = 0.0d0
         y1(:)  = pw(:)
         pw(1:Nr) = matmul(inv_lhs_NS_pol(1:Nr,1:Nr,l), y1(1:Nr))

         call radial_dr_ddr_1D_n2r(pw, pw_dr, pw_ddr)

         green(1:Nr,l)  = pw(1:Nr)
         greenD(1:Nr,l) = pw_dr(1:Nr)
         greenS(1:Nr,l) = pw_ddr(1:Nr)
      enddo
   end subroutine mk_green

   !----------------------------------------------------------------------
   !> Apply Green's functions
   subroutine applyGreen()
      implicit none
      integer:: j,l
      double precision:: d1i, d1o, d2i, d2o, di, do, a1, a2, det

      jl_do(j,l)
         if(flowBC(1).eq.FreeSlip) then
            d1i = greenS1(Nr, l)
            d2i = greenS2(Nr, l)
            di  = flow_pol_ddr(l, j, Nr)
         else
            d1i = greenD1(Nr, l)
            d2i = greenD2(Nr, l)
            di  = flow_pol_dr(l, j, Nr)
         endif
         if(flowBC(2).eq.FreeSlip) then
            d1o = greenS1(1, l)
            d2o = greenS2(1, l)
            do  = flow_pol_ddr(l, j, 1)
         else
            d1o = greenD1(1, l)
            d2o = greenD2(1, l)
            do  = flow_pol_dr(l, j, 1)
         endif

         det = d1i*d2o - d1o*d2i
         a1  = (-d2o*di + d2i*do)/det
         a2  = ( d1o*di - d1i*do)/det
         !-- copy new fields and radial flow_tor_dr_new. to old arrays (we are in (lmr) now):
         flow_pol(l,j,:)     = flow_pol(l,j,:)     + a1*green1(:,l)  + a2*green2(:,l)
         flow_pol_dr(l,j,:)  = flow_pol_dr(l,j,:)  + a1*greenD1(:,l) + a2*greenD2(:,l)
         flow_pol_ddr(l,j,:) = flow_pol_ddr(l,j,:) + a1*greenS1(:,l) + a2*greenS2(:,l)
         !-- time step done for this (lm).
      jl_enddo
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
