! -----------
!> Computes the plots require for the Yokoi paper
! -----------
program YokoiPlots
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

   double precision, allocatable:: u_r(:,:,:), rotu_r(:,:,:), B_r(:,:,:), rotB_r(:,:,:)
   double precision, allocatable:: u_t(:,:,:), rotu_t(:,:,:), B_t(:,:,:), rotB_t(:,:,:)
   double precision, allocatable:: u_p(:,:,:), rotu_p(:,:,:), B_p(:,:,:), rotB_p(:,:,:)
   double precision, allocatable:: uS_r(:,:,:), rotuS_r(:,:,:), BS_r(:,:,:), rotBS_r(:,:,:)
   double precision, allocatable:: uS_t(:,:,:), rotuS_t(:,:,:), BS_t(:,:,:), rotBS_t(:,:,:)
   double precision, allocatable:: uS_p(:,:,:), rotuS_p(:,:,:), BS_p(:,:,:), rotBS_p(:,:,:)
   double precision, allocatable:: aux(:,:,:)
   double precision, allocatable:: EMF(:,:,:)
   double precision, allocatable:: alphaBz(:,:,:)
   double precision, allocatable:: betaJz(:,:,:)
   double precision, allocatable:: gammaOz(:,:,:)
   double precision:: c, c2
   integer:: i, l, j
   integer:: error, nstates
   character:: first

   call init(error)

   call computeAndSaveAverage(nstates)

   !--main loop
   main_loop:do
      read (700,'(A50)', iostat=error) io_calc_file_in
      if (error.ne.0) exit
      io_calc_file_in = adjustl(io_calc_file_in)
      first = io_calc_file_in(1:1)
      if((first=='#').or.(first=='*')) cycle

      ! Read a new state
      call system (inflate_state//trim(io_calc_file_in))
      call drs_load_state(error)
      if(error>0) call drs_abort(error)
      call system (deflate_state//trim(io_calc_file_in))

      ! Remove the steady part
      flow_pol = flow_pol - flow_pol_avg
      flow_tor = flow_tor - flow_tor_avg
      if(field_evolves) then
         field_pol = field_pol - field_pol_avg
         field_tor = field_tor - field_tor_avg
      endif

      ! Initialize derivatives and quantities in real space
      call drs_flow_init(error)
      call calc_flow(u_r, u_t, u_p)
      call calc_rot_flow(rotu_r, rotu_t, rotu_p)
      if(field_evolves) then
         call drs_field_init(error)
         call calc_field(B_r, B_t, B_p)
         call calc_rot_field(rotB_r, rotB_t, rotB_p)
      endif

      ! Compute the z-alligned turbulent electromotive force
      call computeEMF(EMF)

      ! Compute alpha
      call computeAlpha(alphaBz)

      ! Compute beta
      call computeBeta(BetaJz)

      ! Compute gamma
      call computeGamma(gammaOz)

   enddo main_loop
   close(700)

   EMF  = EMF/nstates

   call saveDXmeridional(EMF, 0.0d0, 'EMF')

   ! Compute B_z
   do l=0, Nt
      c  = costheta(l)
      c2 = -sqrt(1.0d0 - c**2)
      do j=1, Np
         do i=1, Nr
            aux(l,j,i) = BS_r(l,j,i)*c + BS_t(l,j,i)*c2
         enddo
      enddo
   enddo
   ! Compute  alpha*Bz:
   alphaBz = alphaBz*aux/nstates
   call saveDXmeridional(alphaBz, 0.0d0, 'alphaBz')

   ! Compute J_z:
   do l=0,Nt
      c  = costheta(l)
      c2 = -sqrt(1.0d0-c**2)
      do j=1,Np
         do i=1,Nr
            aux(l,j,i) = rotBS_r(l,j,i)*c + rotBS_t(l,j,i)*c2
         enddo
      enddo
   enddo
   ! Compute -beta*Jz:
   betaJz  = -betaJz*aux/nstates
   call saveDXmeridional(betaJz, 0.0d0, 'betaJz')

   ! Compute Omega_z
   do l=0,Nt
      c  = costheta(l)
      c2 =-sin(acos(costheta(l)))
      do j=1,Np
         do i=1,Nr
            aux(l,j,i) = rotuS_r(l,j,i)*c + rotuS_t(l,j,i)*c2
         enddo
      enddo
   enddo
   ! Compute gamma*Omega_z:
   gammaOz = gammaOz*aux/nstates
   call saveDXmeridional(gammaOz, 0.0d0, 'gammaOz')


contains
   subroutine init(error)
      implicit none
      integer, intent(inout):: error
      character:: first
      error = 0
      open(700,file='Yokoi.in',status="OLD")
      do
         read(700,'(A)', iostat=error) io_calc_file_in
         io_calc_file_in = adjustl(io_calc_file_in)
         first = io_calc_file_in(1:1)
         if((first.ne.'#').and.(first.ne.'*')) exit
      enddo
      rewind(700)
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
      call check_dims(error)
      if (error.ne.0) call drs_abort(error)

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

      allocate(u_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(u_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(u_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(B_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(B_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(B_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(uS_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(uS_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(uS_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(BS_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(BS_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(BS_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotu_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotu_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotu_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotB_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotB_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotB_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotuS_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotuS_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotuS_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotBS_r(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotBS_t(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(rotBS_p(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(aux(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(EMF(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(alphaBz(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(betaJz(0:(blk_t_size(mpi_rank)-1),Np,Nr))
      allocate(gammaOz(0:(blk_t_size(mpi_rank)-1),Np,Nr))

   end subroutine init

   subroutine computeAndSaveAverage(nstates)
      implicit none
      integer, intent(out):: nstates
      io_calc_file_out = "Steady.0"
      nstates = 0
      !--main loop
      ! Compute the time average
      do
         read (700,'(A50)', iostat=error) io_calc_file_in
         if (error.ne.0) exit
         io_calc_file_in = adjustl(io_calc_file_in)
         first = io_calc_file_in(1:1)
         if((first=='#').or.(first=='*')) cycle
         nstates = nstates + 1

         call system (inflate_state//trim(io_calc_file_in))
         call drs_load_state(error)
         if(error>0) call drs_abort(error)
         call system (deflate_state//trim(io_calc_file_in))

         flow_pol_avg = flow_pol_avg + flow_pol
         flow_tor_avg = flow_tor_avg + flow_tor
         if(field_evolves) then
            field_pol_avg = field_pol_avg + field_pol
            field_tor_avg = field_tor_avg + field_tor
         endif
      enddo
      rewind(700)
      flow_pol = flow_pol_avg/nstates
      flow_tor = flow_tor_avg/nstates
      if(field_evolves) then
         field_pol = field_pol_avg/nstates
         field_tor = field_tor_avg/nstates
      endif
      ! Compute the zonal average
      flow_pol(:,2:Np_s,:) = 0.0d0
      flow_tor(:,2:Np_s,:) = 0.0d0
      if(field_evolves) then
         field_pol(:,2:Np_s,:) = 0.0d0
         field_tor(:,2:Np_s,:) = 0.0d0
      endif
      ! Save the steady zonally averaged state
      call drs_save_state()
      call system(deflate_state//trim(io_calc_file_out))
      ! And store it for future use.
      flow_pol_avg = flow_pol
      flow_tor_avg = flow_tor
      if(field_evolves) then
         field_pol_avg = field_pol
         field_tor_avg = field_tor
      endif
      ! Initialize derivatives and steady quantities in real space
      call drs_flow_init(error)
      call calc_flow(uS_r, uS_t, uS_p)
      call calc_rot_flow(rotuS_r, rotuS_t, rotuS_p)
      if(field_evolves) then
         call drs_field_init(error)
         call calc_field(BS_r, BS_t, BS_p)
         call calc_rot_field(rotBS_r, rotBS_t, rotBS_p)
      endif
   end subroutine computeAndSaveAverage

   subroutine computeEMF(EMF)
      implicit none
      integer:: i,l,j
      double precision:: c, c2
      double precision, intent(inout):: EMF(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: aux(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      ! <u' x B'>_z
      do j=1, Np
         do l=0, Nt
            c  = costheta(l)
            c2 = -sqrt(1.0d0 - c**2)
            do i=1, Nr
               aux(l,j,i)  = &
                  (u_t(l,j,i)*B_p(l,j,i) - u_p(l,j,i)*B_t(l,j,i) )*c + &
                  (u_p(l,j,i)*B_r(l,j,i) - u_r(l,j,i)*B_p(l,j,i) )*c2
            enddo
         enddo
      enddo
      !  Take the phi average
      do l=0, Nt
         do i=1, Nr
            aux(l,1,i) = sum(aux(l,1:Np,i))/Np
         enddo
      enddo
      ! Average in time
      do j=1,Np
            EMF(:,j,:) = EMF(:,j,:) + aux(:,1,:)
      enddo
   end subroutine

   subroutine computeAlpha(alpha)
      implicit none
      double precision, intent(inout):: alpha(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: aux(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      integer:: i,l

      aux = ( B_r*rotB_r + B_t*rotB_t + B_p*rotB_p ) - &
            ( u_r*rotu_r + u_t*rotu_t + u_p*rotu_p )

      ! Take the phi average
      do i=1,Nr
         do l=0,Nt
            aux(l,1,i) = sum(aux(l,1:Np,i))/Np
         enddo
      enddo
      ! Average in time
      do j=1,Np
         alpha(:,j,:) = alpha(:,j,:) + aux(:,1,:)
      enddo
   end subroutine

   subroutine computeBeta(beta)
      implicit none
      double precision, intent(inout):: beta(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: aux(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      integer:: i,l

      aux = ( ( B_r**2 + B_t**2 + B_p**2 ) + &
              ( u_r**2 + u_t**2 + u_p**2 ) )

      ! Take the phi average
      do l=0,Nt
         do i=1,Nr
            aux(l,1,i) = sum(aux(l,1:Np,i))/Np
         enddo
      enddo
      ! Average in time
      do j=1,Np
         beta(:,j,:) = beta(:,j,:) + aux(:,1,:)
      enddo
   end subroutine

   subroutine computeGamma(gamma)
      implicit none
      double precision, intent(inout):: gamma(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision:: aux(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      integer:: i,l

      aux = u_r*B_r + u_t*B_t + u_p*B_p

      ! Take the phi average
      do i=1,Nr
         do l=0,Nt
            aux(l,1,i) = sum(aux(l,1:Np,i))/Np
         enddo
      enddo
      ! Average in time
      do j=1,Np
         gamma(:,j,:) = gamma(:,j,:) + aux(:,1,:)
      enddo
   end subroutine

   subroutine saveDXmeridional(field, phi, filename)
      implicit none
      double precision, intent(in):: field(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: phi
      character(len=*), intent(in):: filename
      double precision:: render_out(0:(blk_t_size(mpi_rank)-1),Nr)
      double precision:: phi_n, dphi, dphi1, dphi2
      double precision:: mcut_min, mcut_max, mcut_rms
      integer:: i,l,icut

      dphi  = 2*pi/m0/Np
      icut  = int(phi/dphi) + 1
      phi_n = icut*dphi
      dphi1 = phi - phi_n
      dphi2 = dphi - dphi1

      ! Interpolate in phi
      render_out = (field(:,icut,:)*dphi2 + field(:,icut+1,:)*dphi1)/dphi

      ! find minimum and maximum
      mcut_min = minval(render_out)
      mcut_max = maxval(render_out)
      mcut_rms = sum(render_out(:,:)*dOmega(:,:))/(pi*(rcoll(1)-rcoll(Nr)))
      mcut_rms = sqrt(mcut_rms)

      open(2,file=trim(filename)//'_rms.txt',STATUS="UNKNOWN")
      Write(2,*) mcut_rms
      close(unit=2)

      open(unit=300, file=filename//'.general', status='UNKNOWN')

      ! Write the DX header
      write(300,*) 'file = ', filename//'.dat'
      write(300,"(' grid = ',I5,' x ',I5)") Nr, Nt+3
      write(300,*) 'format = ascii'
      write(300,*) 'interleaving = field'
      write(300,*) 'majority = column'
      write(300,*) 'field = locations, field0'
      write(300,*) 'structure = 2-vector, scalar'
      write(300,*) 'type = float, float'
      write(300,*) 'end'
      close(300)
      ! And now the content
      open(unit=300, file=filename//'.dat', status='UNKNOWN')
      do i=1,Nr
         write(300,'(3E15.5)') 0.0d0, rcoll(i), render_out(Nt,i)
      enddo
      do l=Nt,0,-1
         do i=1,Nr
            write(300,'(3E15.5)') rcoll(i)*sqrt(1-costheta(l)**2), rcoll(i)*costheta(l), render_out(l,i)
         enddo
      enddo
      do i=1,Nr
         write(300,'(3E15.5)') 0.0d0, -rcoll(i), render_out(0,i)
      enddo

      close(300)

   end subroutine saveDXmeridional

end program

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
