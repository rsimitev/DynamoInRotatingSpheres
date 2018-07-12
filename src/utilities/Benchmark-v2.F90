! -----------
!> Computes quantities relevant to the 2nd benchmark exercise based on a given state.
! -----------
program Benchmarkv2
#include "drsDefs.F90"
   use drs_time
   use drs_params
   use drs_dims
   use drs_mpi
   use drs_fftw3
   use drs_legendre
   use drs_radial
   use drs_hypDiff
   use drs_flow
   use drs_field
   use drs_temp
   use drs_io_state
   use drs_probes
   use drs_real_space
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


   error = 0
   Ekin = -1.0d-100
   Eb = -1.0d-100
   up = -1.0d-100
   T  = -1.0d-100
   Bt = -1.0d-100
   call init(error)
   rr = 0.50d0*(1.0d0+eta)/(1.0d0-eta)

   ! Compute the adimensional volume of the shell
   vol  = volume(eta)
   mmax = (Np_s - 1)/2
   Res  = Nr**(1.0d0/3.0d0)*(Nt_s*Np_s - mmax**2 + mmax + 1)**(1.0d0/3.0d0)
   dp   = 2.0d0*pi/Np

   ! Compute the magnetic energy density
   if(field_evolves) then
      call calc_field(Bvec_r, Bvec_t, Bvec_p)
      Eb = energy(Bvec_r, Bvec_t, Bvec_p)
   endif

   ! Compute the kinetic energy density
   call calc_flow(vec_r, vec_t, vec_p)
   Ekin = energy(vec_r, vec_t, vec_p)

   !Compute the temperature in real space
   call calc_temp(scal)

   ! Construct the flow at the equator and mid-shell
   ! TODO deal with the case mod(mpi_size,2)==0
   !if((blk_t_start(mpi_rank) .le. tcut) .and. (blk_t_start(mpi_rank)+blk_t_size(mpi_rank) .gt. tcut)) then
      call selectEquatorMidShell(vec_r, ur, rcut)
      ! Select the appropriate point for the benchmark quantities
      open(unit=30, file='Bench.'//trim(io_calc_file_in)//'.ur.dat', status='UNKNOWN')
      do i=1, Np
         Write(30,*) (i-1)*dp, ur(i)
      enddo
      close(30)
      pcut = -1
      do i=1, Np
        dur = -1.0d0
        ! The point for pcut must be a zero (or close to) for ur
        if(i.ne.Np) then
           if(ur(i)*ur(i+1).lt.0.0d0) dur = ur(i+1)-ur(i)
        else
           if(ur(Np)*ur(1).lt.0.0d0) dur = ur(1)-ur(Np)
        endif
        ! It also must have dur/dp > 0
        if (dur.gt.0.0d0) then
           pcut = i ! In fact, this point is somewhere between i and i+1
           dp1  = abs(ur(pcut))*dp/dur
           dp2  = dp - dp1
           exit
        endif
      enddo
      ! Compute the quantities
      if(pcut.ne.-1) then
         ! Reuse ur to store u_p
         call selectEquatorMidShell(vec_p, ur, rcut)
         open(unit=30, file='Bench.'//trim(io_calc_file_in)//'.up.dat', status='UNKNOWN')
         do i=1, Np
            Write(30,*) (i-1)*dp, ur(i)*Pm
         enddo
         close(30)
         up = (ur(pcut)*dp2 + ur(pcut+1)*dp1)/dp

         call selectEquatorMidShell(scal, ur, rcut)
         call cacheTemperatureProfile(rcut, T_prof)
         open(unit=30, file='Bench.'//trim(io_calc_file_in)//'.T.dat', status='UNKNOWN')
         do i=1, Np
            Write(30,*) (i-1)*dp, ur(i) + T_prof
         enddo
         close(30)
         T = (ur(pcut)*dp2 + ur(pcut+1)*dp1)/dp
         ! This is the adimensional temperature for fixed value BC's
         T = T + T_prof

         if(field_evolves) then
            call selectEquatorMidShell(Bvec_t, ur, rcut)
            open(unit=30, file='Bench.'//trim(io_calc_file_in)//'.Bt.dat', status='UNKNOWN')
            do i=1, Np
               Write(30,*) (i-1)*dp, -ur(i)*sqrt(Pm/Ta)
            enddo
            close(30)
            Bt = (ur(pcut)*dp2 + ur(pcut+1)*dp1)/dp
         endif
      endif

   !endif
   !call drs_maximize(up)
   !call drs_maximize(T)
   !call drs_maximize(Bt)

   spew 'Resolution'
   spew Res
   spew 'Position: pcut, phi(pcut), phi'
   spew pcut, pcut*dp, pcut*dp + dp1
   spew 'Time'
   spew time/Pm
   spew 'Ek, Eb, T, up, Bt'
   Write(*,'(5F15.5)') Ekin*Pm**2, Eb*Pm**2, T, up*Pm, -Bt*sqrt(Pm/Ta)

contains
   !> Cache the temperature from the profile
   subroutine cacheTemperatureProfile(rcut, T)
      implicit none
      integer, intent(in):: rcut
      double precision, intent(out):: T
      double precision:: halfrad, dr1, dr2
      halfrad = 0.50d0*(1.0d0+eta)/(1.0d0-eta)
      T = temp_profile(rcut)

      if(rcoll(rcut).gt.halfrad) then !No point at mid-shell
         dr1 = rcoll(rcut)-halfrad
         dr2 = halfrad-rcoll(rcut+1)

         T = (T*dr2 + temp_profile(rcut+1)*dr1)/drcoll(rcut)
      endif
      T = T + 1.0d0/Pt

   end subroutine cacheTemperatureProfile


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


      !TODO replace with a call to readparv2
      ! Read computation name
      open(700,file='benchmarkv1.in',status="OLD")
      do
         read(700,'(A50)', iostat=error) io_calc_file_in
         io_calc_file_in = adjustl(io_calc_file_in)
         first = io_calc_file_in(1:1)
         if((first.ne.'#').and.(first.ne.'*')) exit
      enddo
      close(700)
      ! Read parameters in
      open(900, file=trim(io_calc_file_in)//'.par', status='OLD')
      call drs_read_state_par(900) ! This is needed before initialising the modules
      close(900)

      ! Initialise everything
      drs_want_hypDiff = .FALSE.

      ! Set the resolution to use in the calculations.
      m0 = m0i
      Nr = Nri
      Nt = Nti
      Np = Npi
      eta = etai
      Nr_s = Nri_s
      Nt_s = Nti_s
      Np_s = Npi_s
      lform = lformi
      drift = drifti
      comment = commenti
      drs_calc_type = drs_calc_typei
      ! Update the simulation parameters
      Pt   = Pti
      Ta   = Tai
      Ra_t = Ra_ti
      Pm   = Pmi

      call drs_mpi_init(error)
      if(mpi_size.ne.1) then
         spew 'This program should be ran on a single cpu'
         call drs_abort(1)
      endif
      call select_calc_type(error)
      ! Start the initializations
      call drs_time_init()
      call drs_params_init()
      call drs_dims_init(error)

      call mpi_dims_init(Nt, Np_s, m0, error)
      call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)

      call drs_legendre_allocation()
      call drs_flow_allocation()
      call drs_field_allocation()
      call drs_temp_allocation()
      call drs_probes_allocation()

      call drs_legendre_init()
      call drs_radial_init(eta)
      call drs_hypDiff_init(Nt)
      call drs_probes_init(time)

      call system (trim(inflate)//' '//trim(io_calc_file_in))
      call drs_load_state(error)
      if(error>0) call drs_abort(error)
      call system (trim(deflate)//' '//trim(io_calc_file_in))

      ! Initialize derivativs
      call drs_flow_init(error)
      if(field_evolves) call drs_field_init(error)
      call drs_temp_init(error)

      allocate(vec_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(vec_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(vec_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(Bvec_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(Bvec_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(Bvec_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(scal(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(ur(Np))

   end subroutine
end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
