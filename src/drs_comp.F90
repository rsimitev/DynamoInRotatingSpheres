! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Module dealing with the composition.
module drs_comp
#include "drsDefs.F90"
   use drs_mpi
   use drs_dims
   use drs_params
   use drs_legendre
   use drs_transforms
   use drs_radial
   use drs_hypDiff
   implicit none
   save

   double precision, allocatable:: comp(:,:,:)
   double precision, allocatable:: comp_dr(:,:,:)
   double precision, allocatable:: comp_ddr(:,:,:)
   double precision, allocatable:: comp_avg(:,:,:)
   double precision, allocatable:: comp_dr_avg(:)
   double precision, allocatable:: comp_profile(:)
   double precision, allocatable:: comp_profile_dr(:)
   double precision, allocatable:: comp_t(:,:,:)

   logical, private:: initialised=.false.
   logical, private:: comp_allocated=.false.
contains

   !---------------------------------------------------------------
   !> Allocates the variables required for computations envolving composition.
   subroutine drs_comp_allocation()
      implicit none
      if(comp_allocated) return
      allocate(  comp_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate(    comp(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate( comp_dr(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(comp_ddr(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(comp_avg(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(comp_dr_avg(Nr))
      allocate(comp_profile(Nr))
      allocate(comp_profile_dr(Nr))
      comp_allocated = .true.
   end subroutine drs_comp_allocation

   !---------------------------------------------------------------
   !> Initialises the composition boundary conditions, derivatives and profiles.
   subroutine drs_comp_init(error)
      implicit none
      integer, intent(out):: error
      integer::i
      error=0
      if(initialised) then
         error=-1
         return
      endif
      call drs_bcast(compBC, 2)
      call drs_bcast(compProf)
      call radial_dr_ddr_3D_r2r(comp, comp_dr, comp_ddr)
      select case(compProf)
         ! Well mixed core
         case( WellMixed )
            ! Composition is uniform accross the shell.
            comp_profile    = 1.0d0
            comp_profile_dr = 0.0d0
         case( Diffusive )
            ! Pure diffusion profile considering fixed values at the surface
            forall(i=1:Nr)
               comp_profile(i)    = ( eta/rcoll(i) - (1-eta) )/(1-eta)**2/Pc
               comp_profile_dr(i) = -eta/rcoll2(i)/(1-eta)**2/Pc
            endforall
         case( InternalSources )
            ! Internal composition sources (don't ask!)
            forall(i=1:Nr)
               comp_profile(i)    = (1.0d0 - 0.5d0*rcoll2(i))/Pc
               comp_profile_dr(i) = -rcoll(i)/Pc
            endforall
         case default
            error = 100
            initialised=.false.
            return
      endselect
      initialised=.true.
   end subroutine

   !---------------------------------------------------------------
   !> Outputs a human readable name for the composition profiles.
   !! @since 1.6.1
   function compProfName()
      implicit none
      character(len=16):: compProfName
      select case(compProf)
         ! Well mixed core
         case(WellMixed)
            compProfName = 'well mixed'
         case(Diffusive)
            compProfName = 'diffusive'
         case(InternalSources)
            compProfName = 'internally gen.'
         case default
            compProfName = 'diffusive'
      endselect
   end function

   !---------------------------------------------------------------
   !> Resets the composition and its derivatives to 0.
   subroutine drs_comp_reset()
      implicit none
      comp     = 0.0d0
      comp_dr  = 0.0d0
      comp_ddr = 0.0d0
   end subroutine

   !---------------------------------------------------------------
   !> Computes the laplacian of the composition.
   double precision pure function comp_lap(l,j,i)
      implicit none
      integer, intent(in)::l,j,i
      comp_lap = comp_ddr(l, j, i) + 2.0/rcoll(i)*comp_dr(l, j, i) - llp1(l)/rcoll2(i)*comp(l, j, i)
      call drs_apply_hypDiff(comp_lap, l)
   end function

   !---------------------------------------------------------------
   subroutine drs_comp_randomize(noise)
      implicit none
      double precision, intent(in):: noise
      double precision:: ranr
      integer:: ran
      integer:: i, j, jg, l, m
      if(noise.le.0.D0) return
      ran=1
      do j=2,min(6,Np_s)
         jg = blk_ps_start(mpi_rank) + j
         m  = m0*(jg/2)
            do l=max(1,m),min(4*m0,Nt_s)
               ran=mod(ran*106+1283,6075)
               ranr=1.0*ran/6075.0
               do i=2,Nr-1
                  if (noise.eq.100) then
                     comp(l,j,i) = comp(l,j,i) + ranr*comp(l,j,Nr/2)/(l+1)**2/(m+1)**4*sin(pi*(rcoll(i)-eta/(1-eta)))
                  else
                     comp(l,j,i) = comp(l,j,i) + noise*ranr/(l+1)**2/(j/2)**4*sin(pi*(rcoll(i)-eta/(1-eta)))
                  endif
            enddo
         enddo
      enddo
   end subroutine

   !---------------------------------------------------------------
   !>  These lines take care of boundary conditions
   !!  If the value at a boundary is bc different from 0, insert bc/2
   !!  because of the factor 2 between Chebychev and radial differentiation
   !!  If the boundary condition is not on the derivative but the value of the
   !!  function itself (as is the case for the temperature), no factor of 2 is needed
   !!  if the boundary value bc is different from 0.
   subroutine apply_comp_BC(pencil)
      implicit none
      !> A pencil with the forces in lmr space.
      double precision, intent(inout):: pencil(Nr)
      pencil(1)  = 0.0d0
      pencil(Nr) = 0.0d0
   end subroutine

   !---------------------------------------------------------------
   !> Computes the composition in real space from the composition in spectral
   !! space.
   subroutine calc_comp(comp_spec,comp_real)
      implicit none
      !> Composition in real space.
      double precision, intent(out):: comp_real(0:blk_t_size(mpi_rank)-1,Np,Nr)
      !> Composition in spectral space.
      double precision, intent(in):: comp_spec(0:Nt_s,1:blk_ps_size(mpi_rank),Nr)
      integer:: i
      double precision:: aux1(0:Nt_s,blk_ps_size(mpi_rank))
      double precision:: aux2(0:blk_t_size(mpi_rank)-1,Np)

      do i=1,Nr
         aux1(0:Nt_s,1:blk_ps_size(mpi_rank)) = comp_spec(0:Nt_s,1:blk_ps_size(mpi_rank),i)
         call ylmb(aux1,aux2)
         comp_real(0:blk_t_size(mpi_rank)-1,1:Np,i) = aux2(0:blk_t_size(mpi_rank)-1,1:Np)
      enddo
   end subroutine
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
