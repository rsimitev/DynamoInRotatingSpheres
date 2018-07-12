! Copyright  L. Silva (lacsilva@gmail.com), 2015
#include "drsDefs.f90"
! -----------
!> Computes quantities relevant to the 1st benchmark exercise based on a given state.
! -----------
program test_calc_flow
   use drs_params
   use drs_dims
   use drs_mpi
   use drs_fftw3
   use drs_legendre
   use drs_radial
   use drs_flow
   implicit none

   double precision:: vol, Res
   double precision, allocatable:: vec_r(:,:,:), Bvec_r(:,:,:)
   double precision, allocatable:: vec_t(:,:,:), Bvec_t(:,:,:)
   double precision, allocatable:: vec_p(:,:,:), Bvec_p(:,:,:)
   double precision, allocatable:: scal(:,:,:)
   double precision, allocatable:: ur(:)
   double precision:: dur, up, T, Bt, T_prof
   double precision:: dp, dp1, dp2, rr
   integer:: error, rcut, pcut, i, mmax


   call init(error)
   if(error.ne.0) call drs_abort(error)
   ! At this point we need to initialise the flow

   ! Initialize derivatives
   call drs_flow_init(flow_tor_dr,flow_tor_ddr,flow_pol_dr,flow_pol_ddr)
   ! Compute the kinetic energy density
   call calc_flow(vec_r, vec_t, vec_p)
   ! Test vector components at individual points
   Ekin = energy(vec_r, vec_t, vec_p)
   ! test the kinetic energy


contains


   subroutine selectEquatorMidShell(field, line, rcut)
      implicit none
      double precision, intent(in)::  field(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(out):: line(Np)
      double precision:: halfrad, dr1, dr2, dct1, dct2, dct
      integer, intent(out):: rcut
      integer:: tcut, tstart, i
      halfrad = 0.50d0*(1.0d0+eta)/(1.0d0-eta)

      ! Radius starts with ro and decreases with i
      rcut = 1
      do i=1, Nr
         if (rcoll(i).ge.halfrad) then
            rcut = i ! Just above mid-shell. rcoll(rcut+1) is just below mid-shell
         else
            exit
         endif
      enddo
      ! costheta varies between -1 and 1 and increases with i, therefore
      ! theta decreases with i
      tcut = 0
      do i=0, Nt
         if (costheta(i).le.0.0d0) then
            tcut = i ! Just below the equator. costheta(tcut+1) is just above equator
         else
            exit
         endif
      enddo
      tstart = blk_t_start(mpi_rank)

      line(:) = field(tcut-tstart,:,rcut)

      if(rcoll(rcut).gt.halfrad) then !No point at mid-shell
         dr1 = rcoll(rcut)-halfrad
         dr2 = halfrad-rcoll(rcut+1)

         if (abs(costheta(tcut)).lt.1.0d-50) then !There is a point at the equator
            line(:) = (line(:)*dr2 + field(tcut-tstart,:,rcut+1)*dr1)/drcoll(rcut)
         else !There is no point at the equator
            dct1 = -costheta(tcut)
            dct2 = costheta(tcut+1)
            dct  = dct1 + dct2
            line(:) = (line(:)*dr2*dct2 + &
                  field(tcut-tstart,:,rcut+1)*dr1*dct2 + &
                  field(tcut-tstart+1,:,rcut+1)*dr1*dct1 + &
                  field(tcut-tstart+1,:,rcut)*dr2*dct1&
                  )/(drcoll(rcut)*dct)
         endif
      else
         if (abs(costheta(tcut)).lt.1.0d-50) then !There is no point at the equator
            dct1 = -costheta(tcut)
            dct2 = costheta(tcut+1)
            dct  = dct1 + dct2
            line(:) = (line(:)*dct2 + field(tcut-tstart+1,:,rcut)*dct1)/dct
         endif
      endif
   end subroutine selectEquatorMidShell

   double precision function volume(eta)
      implicit none
      double precision, intent(in):: eta

      volume = 4.0d0*pi/3.0d0*(1.0d0-eta**3)/(1.0d0-eta)**3
   end function volume

   subroutine init(error)
      implicit none
      integer, intent(inout):: error
      character(len=100):: deflate, inflate
      character:: first
      error = 0
      deflate = 'compressdata '
      inflate = 'uncompressdata '

      ! Initialise everything
      drs_want_hypDiff = .FALSE.

      ! Set the resolution to use in the calculations.
      m0 = 1
      Nr = 5
      Nt = 13
      Np = 26
      eta = 0.5
      Nr_s = 5
      Nt_s = 13
      Np_s = 26
      lform = 1
      drift = 0
      comment = ''
      drs_calc_type = 3

      call drs_mpi_init()
      if(mpi_size.ne.1) then
         spew 'This program should be ran on a single cpu'
         call drs_abort(1)
      endif
      call select_calc_type(error)
      ! Start the initializations
      call drs_params_init()
      call drs_dims_init(error)

      call mpi_dims_init(Nt, Np_s, m0, error)
      call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)

      call drs_legendre_allocation()
      call drs_flow_allocation()
      call drs_probes_allocation()

      call drs_legendre_init()
      call drs_radial_init(eta)


      allocate(vec_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(vec_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(vec_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(ur(Np))

   end subroutine
end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
