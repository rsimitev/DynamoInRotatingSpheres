! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Temperature related operations
module drs_temp
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

   double precision, allocatable:: temp(:,:,:)
   double precision, allocatable:: temp_dr(:,:,:)
   double precision, allocatable:: temp_ddr(:,:,:)
   double precision, allocatable:: temp_lap(:,:,:)
   double precision, allocatable:: temp_avg(:,:,:)
   double precision, allocatable:: temp_dr_avg(:)
   double precision, allocatable:: temp_profile(:)
   double precision, allocatable:: temp_profile_dr(:)
   double precision, allocatable:: temp_t(:,:,:)

   logical, private:: initialised=.false.
   logical, private:: temp_allocated=.false.
contains
   !------------------------------------------------------------------
   !> Allocates the temperature related variables.
   subroutine drs_temp_allocation()
      implicit none
      if(temp_allocated) return
      allocate(    temp(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(  temp_t(0:(blk_t_size(mpi_rank)-1), Np, Nr))
      allocate( temp_dr(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(temp_ddr(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(temp_lap(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(temp_avg(0:Nt_s,1:blk_ps_size(mpi_rank),Nr))
      allocate(temp_dr_avg(Nr))
      allocate(temp_profile(Nr))
      allocate(temp_profile_dr(Nr))
      temp_allocated = .true.
   end subroutine drs_temp_allocation

   !------------------------------------------------------------------
   !> Precomputes the adimensional radial temperature profile.
   subroutine drs_temp_init(error)
      implicit none
      integer, intent(out):: error
      integer:: i
      error=0
      if(initialised) then
         error=-1
         return
      endif
      ! If we have not been previously allocated, allocate and
      ! initialise to 0
      if(.not.temp_allocated) then
         call drs_temp_allocation()
         temp = 0.0d0
         error=-2 ! The flow was allocated and initialised here.
      endif
      call drs_bcast(tempBC, 2)
      call drs_bcast(tempProf)
      call radial_dr_ddr_3D_r2r(temp,temp_dr,temp_ddr)
      call update_temp_lap()
      if(tempProf.eq.Conduction) then
         ! No heat sources
         forall(i=1:Nr)
            temp_profile(i)    = ( eta/rcoll(i) - (1-eta) )/(1-eta)**2/Pt
            temp_profile_dr(i) = -eta/rcoll2(i)/(1-eta)**2/Pt
         endforall
      elseif(tempProf.eq.InternalHeating) then
         ! Internal heat sources:
         forall(i=1:Nr)
            temp_profile(i)    = (1.0d0 - 0.5d0*rcoll2(i))/Pt
            temp_profile_dr(i) = -rcoll(i)/Pt
         endforall
      endif
   end subroutine

   !------------------------------------------------------------------
   !> Outputs a human readable name for the temperature profiles.
   !! @since 1.6.1
   function tempProfName()
      implicit none
      character(len=16):: tempProfName
      select case(tempProf)
         case(Conduction)
            tempProfName = 'conductive'
         case(InternalHeating)
            tempProfName = 'internal heating'
         case default
            tempProfName = 'internal heating'
      endselect
   end function

   !------------------------------------------------------------------
   !> Resets the temperature and its derivatives to zero.
   subroutine drs_temp_reset()
      implicit none
      temp  = 0.0d0
      temp_dr = 0.0d0
      temp_ddr = 0.0d0
   end subroutine

   !------------------------------------------------------------------
   !> Recomputes and caches the laplacian of the temperature.
   subroutine update_temp_lap()
      implicit none
      integer::l,j,i
      do i=1, Nr
         jl_do(j,l)
            temp_lap(l,j,i) = temp_ddr(l,j,i) + 2.0/rcoll(i)*temp_dr(l,j,i) - llp1(l)/rcoll2(i)*temp(l,j,i)
            call drs_apply_hypDiff(temp_lap(l, j, i), l)
         jl_enddo
      enddo
   end subroutine

   !------------------------------------------------------------------
   subroutine drs_temp_randomize(noise)
      implicit none
      double precision, intent(in):: noise
      double precision:: ranr
      integer:: ran
      integer:: i, j, jg, l, m
      if(noise.le.0.D0) return
      ran=1
      do j=2,min(6,Np_s)
         jg = blk_ps_start(mpi_rank) + j - 1
         m  = m0*(jg/2)
            do l=max(1,m),min(4*m0,Nt_s)
               ran=mod(ran*106+1283,6075)
               ranr=1.0*ran/6075.0
               do i=2,Nr-1
                  if (int(noise).eq.100) then
                     temp(l,j,i) = temp(l,j,i) + temp(l,j,Nr/2)*ranr/(l+1)**2/(m+1)**4*sin(pi*(rcoll(i)-eta/(1-eta)))
                  else
                     temp(l,j,i) = temp(l,j,i) + noise*         ranr/(l+1)**2/(j/2)**4*sin(pi*(rcoll(i)-eta/(1-eta)))
                  endif
            enddo
         enddo
      enddo
   end subroutine

   !------------------------------------------------------------------
   !> Boundary conditions for the temperature.
   !! For fixed temperature, the value of the anomaly is set to be zero
   !! \todo make this value depend on l and m when we can specify the
   !! full 2D anomaly at the boundaries.
   subroutine apply_temp_BC_RHS(pencil)
      implicit none
      !> A pencil of temperature data in lmr space.
      double precision, intent(inout):: pencil(Nr)
      ! Assume the temperature deviation is zero at the boundary
      ! Inner boundary
      if(tempBC(1)==FixedTemperature) pencil(Nr) = 0.0d0
      ! Outer boundary
      if(tempBC(2)==FixedTemperature) pencil(1)  = 0.0d0
      ! Assume the flux due to the temperature deviation is zero at the boundary
      ! Inner boundary
      if(tempBC(1)==FixedHeatFlux) pencil(Nr) = 0.0d0
      ! Outer boundary
      if(tempBC(2)==FixedHeatFlux) pencil(1)  = 0.0d0
   end subroutine

   !------------------------------------------------------------------
   !> Computes the temperature anomaly in real space.
   subroutine calc_temp(temp_t)
      implicit none
      double precision, intent(out):: temp_t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      integer:: i
      double precision:: aux1(0:Nt_s,blk_ps_size(mpi_rank))
      double precision:: aux2(0:blk_t_size(mpi_rank)-1,Np)

      do i=1,Nr
         aux1(0:Nt_s,1:blk_ps_size(mpi_rank)) = temp(0:Nt_s,1:blk_ps_size(mpi_rank),i)
         call ylmb(aux1,aux2)
         temp_t(0:blk_t_size(mpi_rank)-1,1:Np,i) = aux2(0:blk_t_size(mpi_rank)-1,1:Np)
      enddo
   end subroutine
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
