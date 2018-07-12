! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Deals with input and output of state files and derived quantities.
module drs_io
#include "drsDefs.F90"
   use drs_params
   use drs_io_state
   use drs_io_units
   use drs_dims
   use drs_time
   use drs_fftw3 ! remesh
   use drs_mpi
   use drs_flow
   use drs_field
   use drs_temp
#ifdef COMP
   use drs_comp
#endif
   use drs_radial ! rcoll et al.
   use drs_real_space !evaluate_real_space
   use drs_probes
   use drs_error_codes
   implicit none
   save
   private
   character(len=13), parameter:: deflate = 'compressdata '
   character(len=15), parameter:: inflate = 'uncompressdata '

   public:: drs_open_output, drs_close_output
   public:: dump_state
   public:: deflate, inflate
   public:: save_l_spec, save_m_spec, save_n_spec, save_cfl

contains

   !---------------------------------------------------
   !> Checks that the symmetries in the input files are compatible
   !! with the symmetries we require from this computation.
   subroutine check_input_par(error)
      implicit none
      integer, intent(out):: error
      if(m0.ne.m0i) then
         if(m0.lt.1 .or. m0.eq.7 .or. m0.eq.11 .or. m0.eq.13 .or. m0.ge.17) then
            write(*,*) 'Invalid symmetry: ',m0
            error = ERR_INVALID_M_SYMMETRY
         endif
         if(mod(m0i,m0).ne.0 .and. mod(m0,m0i).ne.0) then
            spew 'm0i is',m0i,' m0 is',m0
            spew 'This is not very sensible.'
            error = ERR_INCOMPATIBLE_SYMMETRY
         endif
      endif
      Write(*,*) 'Input files are described by the following parameters:'
      !-- dims of input file must match dims of calc. (imperative for radial dir.):
      write(*,'(A11,8A7)') 'Dimensions:','Nr','Nt','Np','Nr_s','Nt_s','Np_s','lsymm','m0'
      write(*,'(A11,8I7)') ' ', Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi, m0i
#ifdef COMP
      write(*,'(8A11)') 'Parameters:','eta','Pt','Pc','Tau','Ra_t','Ra_c','Pm'
      write(*,'(A11,7D11.4)') ' ',etai, Pti, Pci, Tai, Ra_ti, Ra_ci, Pmi
#else
      write(*,'(6A11)') 'Parameters:','eta','Pt','Tau','Ra_t','Pm'
      write(*,'(A11,5D11.4)') ' ',etai, Pti, Tai, Ra_ti, Pmi
#endif
   end subroutine check_input_par

   !----------------------------------------------------------------
   !> Opens units for regularly probed quantities to be saved.
   subroutine drs_open_output()
      implicit none
      !-- open some output files:
      if(mpi_rank.eq.0) then
         open (unit_cfl,file=trim(io_calc_file_out)//'.cfl',status='unknown')
         if(flow_evolves) then
            open (unit_u_mid,file=trim(io_calc_file_out)//'.u_mid',status='unknown')
            open (unit_ek,file=trim(io_calc_file_out)//'.ek',status='unknown')
            open (unit_koeu,file=trim(io_calc_file_out)//'.koeu',status='unknown')
            open (unit_dissu,file=trim(io_calc_file_out)//'.dissu',status='unknown')
            open (unit_am,file=trim(io_calc_file_out)//'.am',status='unknown')
         endif
         if(temp_evolves) then
            open (unit_nu,file=trim(io_calc_file_out)//'.nu',status='unknown')
         endif
         if(field_evolves) then
            open (unit_eb,file=trim(io_calc_file_out)//'.eb',status='unknown')
            open (unit_koeb,file=trim(io_calc_file_out)//'.koeb',status='unknown')
            open (unit_dissB,file=trim(io_calc_file_out)//'.dissB',status='unknown')
         endif
         if(drs_calc_type.eq.KinematicDynamo) then
            open (unit_evp,file=trim(io_calc_file_out)//'.evp',status='unknown')
            open (unit_evt,file=trim(io_calc_file_out)//'.evt',status='unknown')
         endif
      endif
      call wait_for_everyone()
   end subroutine drs_open_output
   
   !----------------------------------------------------------------
   !> Closes the open units for regularly probed quantities.
   subroutine drs_close_output()
      implicit none
      if(mpi_rank.eq.0) then
         close (unit_cfl) !cfl
         if(flow_evolves) then
            close (unit_ur) ! u_r @Nr/2
            close (unit_ek) ! Ek
            close (unit_koeu) ! flow coeffs
            close (unit_dissu) ! flow dissipation
            close (unit_am) ! Angular momentum
         endif
         if(temp_evolves) then
            close (unit_nu) ! Nu
         endif
         if(field_evolves) then
            close (unit_eb) ! Eb
            close (unit_koeb) ! field coefs
            close (unit_dissB) ! magnetic dissipation
         endif
         if(drs_calc_type.eq.KinematicDynamo) then
            close (unit_evp)
            close (unit_evt)
         endif
      endif  !  root
      call wait_for_everyone()
   end subroutine drs_close_output

   !****************************************************************
   !> Saves the final state and writes the final quantities to file
   subroutine dump_state()
      !-- 30.09.96 M.A. partial parallelization.
      implicit none
      !--   workarray  for the averaged fields:
 
      ! Start by saving the state
      call drs_save_state()
      !>\todo  Describe the format of the following files.
      if(mpi_rank.eq.0) then
         if(temp_evolves) then
            ! Advection 
            call save_advection()
            ! Temperatures
            call save_average_temperature()
         endif
         if(flow_evolves) then
            ! Mean radial flow
            call save_mean_radial_flow()
            ! Mean horizontal flow
            call save_mean_horizontal_flow()
         endif
      endif  !  mpi_rank.eq.0

      ! Save final values for energies, etc.
      call evaluate_real_space()
      call save_stuff(nsample)
      !-- Output (normalized) power spectra
      call save_l_spec()
      call save_m_spec()
      call save_n_spec()

      if(flow_evolves) then
    call save_zonal_averages()
      endif

      spew 'state dumped after ',steps,'steps, nsample = ',nsample
      spew 'file = ',trim(io_calc_file_out)
   end subroutine dump_state

   !****************************************************************
   !>
   subroutine save_advection()
      implicit none
      integer:: i
      double precision:: alpha, dt

      dt = time - time_start

      open (unit_adv,file=trim(io_calc_file_out)//'.adv',status='unknown')
      do i=1,Nr
    !-- write averaged fields:
    alpha = rcoll2(i)*(temp_dr_avg(i) - Pt*adv_avg(i))/dt/(-eta/(1-eta)**2) + 1
    write (unit_adv,'(4(D17.9,X))') rcoll(i), adv_avg(i)/dt, temp_dr_avg(i)/dt, alpha
      enddo
      close (unit_adv)
   end subroutine
   
   !****************************************************************
   !>
   subroutine save_average_temperature()
      implicit none
      integer:: i
      double precision::dt
      dt = time - time_start

      open (unit_t,  file=trim(io_calc_file_out)//'.t',  status='unknown')
      do i=1,Nr
    write (unit_t,'(4(D17.9,X))')   rcoll(i), t2_avg(i)/dt,  temp_avg(0,1,i)/dt, temp_profile(i)
      enddo
      close (unit_t)
   end subroutine
   
   !****************************************************************
   !>
   subroutine save_mean_radial_flow()
      implicit none
      integer:: i
      double precision::dt
      dt = time - time_start

      open (unit_ur, file=trim(io_calc_file_out)//'.ur', status='unknown')
      do i=1,Nr
    !-- write averaged fields:
    write (unit_ur,'(2(D17.9,X))')  rcoll(i), ur_avg(i)/dt
      enddo
      close (unit_ur)
   end subroutine

   !****************************************************************
   !>
   subroutine save_mean_horizontal_flow()
      implicit none
      integer:: i
      double precision::dt
      dt = time - time_start

      open (unit_uaz,file=trim(io_calc_file_out)//'.uaz',status='unknown')
      do i=1,Nr
    write (unit_uaz,'(5(D17.9,X))') rcoll(i), ut_avg(i)/dt,  ut2(i), up_avg(i)/dt, up2(i)
      enddo
      close (unit_uaz)
   end subroutine
   
   !****************************************************************
   !> Computes the average zonal flow.
   !! This routine has the side effect of replacing the flow by the 
   !! average flow. It should be used with care.
   subroutine save_zonal_averages()
      implicit none
      double precision:: u_zonal(0:Nt,Nr)
      integer:: i, l
      !-- get the derivative flow_pol_dr of average poloidal flow:
      flow_pol = flow_pol_avg
      call radial_dr_ddr_3D_r2r(flow_pol, flow_pol_dr, flow_pol_ddr)

      flow_tor = flow_tor_avg
      call radial_dr_ddr_3D_r2r(flow_tor, flow_tor_dr, flow_tor_ddr)

      call calc_flow(flow_r_t, flow_t_t, flow_p_t)
      u_zonal = 0.0d0
      do i=1,Nr
    ! Each cpu takes care of its own theta section
    do l=0, blk_t_size(mpi_rank)-1
       u_zonal(l+blk_t_start(mpi_rank),i) = sum(flow_p_t(l,1:Np,i))/Np
    enddo
    ! Merge all the theta sections into the final array.
    call sum_over_all_cpus(u_zonal(:,i))
      enddo
      ! Write everything down
      if(mpi_rank.eq.0) then
    open (unit_uzon, file=trim(io_calc_file_out)//'.uzon', status='unknown')
    do i=1,Nr
       do l=0,Nt
          write (unit_uzon,*) real(u_zonal(l,i))
       enddo
    enddo
    close (unit_uzon)
      endif
   end subroutine

   !-----------------------------------------------------------------------------
   !>  Saves the normalized power spectra with respect to l
   subroutine save_l_spec()
      implicit none
      integer:: l
      double precision:: tspec(0:Nt_s)
#ifdef COMP
      double precision:: cspec(0:Nt_s)
#endif
      double precision:: uspec(0:Nt_s),bspec(0:Nt_s)

      if(temp_evolves)  call l_spec_of_scalar_field(temp,tspec)
#ifdef COMP
      if(comp_evolves)  call l_spec_of_scalar_field(comp,cspec)
#endif
      if(flow_evolves)  call calc_flow_lspec(uspec)
      if(field_evolves) call calc_field_lspec(bspec)

      if(mpi_rank.eq.0) then
         if(temp_evolves.or.flow_evolves) then
            open (unit_lspec,file=trim(io_calc_file_out)//'.lspec',status='unknown')
            if(field_evolves) then
#ifdef COMP
               if(comp_evolves) then
                  do l=1,Nt_s
                     write (unit_lspec,'(I3,4(D17.9,X))') l,tspec(l),cspec(l),uspec(l),bspec(l)
                  enddo
               else
#endif
                  do l=1,Nt_s
                     write (unit_lspec,'(I3,3(D17.9,X))') l,tspec(l),uspec(l),bspec(l)
                  enddo
#ifdef COMP
               endif
#endif
            else
#ifdef COMP
               if(comp_evolves) then
                  do l=1,Nt_s
                     write (unit_lspec,'(I3,3(D17.9,X))') l,tspec(l),cspec(l),uspec(l)
                  enddo
               else
#endif
                  do l=1,Nt_s
                     write (unit_lspec,'(I3,2(D17.9,X))') l,tspec(l),uspec(l)
                  enddo
#ifdef COMP
               endif
#endif
            endif               ! if(field_evolves)
            close (unit_lspec)
         endif
      endif                     !  mpi_rank.eq.0
   end subroutine save_l_spec

   !-----------------------------------------------------------------------------
   !>  Saves the normalized power spectra with respect to m.
   subroutine save_m_spec()
      implicit none

      double precision:: tspec(m0*Np_s+1)
#ifdef COMP
      double precision:: cspec(m0*Np_s+1)
#endif
      double precision:: uspec(m0*Np_s+1)
      double precision:: Bspec(m0*Np_s+1)
      integer:: m, mmax, jg

      if(temp_evolves)  call m_spec_of_scalar_field(temp,tspec)
      if(flow_evolves)  call calc_flow_mspec(uspec)
      if(field_evolves) call calc_field_mspec(Bspec)
#ifdef COMP
      if(comp_evolves)  call m_spec_of_scalar_field(comp,cspec)
#endif

      mmax  = m0*(Np_s/2)
      if(mpi_rank.eq.0) then
         if(temp_evolves.or.flow_evolves) then
            open (unit_mspec,file=trim(io_calc_file_out)//'.mspec',status='unknown')
            if(field_evolves) then
#ifdef COMP
               if(comp_evolves) then
                  do m=0,mmax,m0
                     jg = 2*(m/m0)+1
                     write (unit_mspec,'(I6,4(D17.9,X))') m,tspec(jg),cspec(jg),uspec(jg),Bspec(jg)
                  enddo
               else
#endif
                  do m=0,mmax,m0
                     jg = 2*(m/m0)+1
                     write (unit_mspec,'(I6,3(D17.9,X))') m,tspec(jg),uspec(jg),Bspec(jg)
                  enddo
#ifdef COMP
               endif
#endif
            else
#ifdef COMP
               if(comp_evolves) then
                  do m=0,mmax,m0
                     jg = 2*(m/m0)+1
                     write (unit_mspec,'(I6,3(D17.9,X))') m,tspec(jg),cspec(jg),uspec(jg)
                  enddo
               else
#endif
                  do m=0,mmax,m0
                     jg = 2*(m/m0)+1
                     write (unit_mspec,'(I6,2(D17.9,X))') m,tspec(jg),uspec(jg)
                  enddo
#ifdef COMP
               endif
#endif
            endif
            close (unit_mspec)
         endif                  !  root
      endif
   end subroutine save_m_spec

   !-----------------------------------------------------------------------------
   !>  Saves the normalized power spectra of all quantities with respect to n
   subroutine save_n_spec()
      implicit none
      integer:: i
      double precision:: tspec(Nr_s), uspec(Nr_s), bspec(Nr_s)
#ifdef COMP
      double precision:: cspec(Nr_s)
#endif

      if(flow_evolves)  call calc_flow_nspec(uspec)
      if(field_evolves) call calc_field_nspec(Bspec)
      if(temp_evolves)  call n_spec_of_scalar_field(temp, tspec)
#ifdef COMP
      if(comp_evolves)  call n_spec_of_scalar_field(comp, cspec)
#endif

      ! Output
      if(mpi_rank.eq.0) then
         if(temp_evolves.or.flow_evolves) then
            open (unit_nspec,file=trim(io_calc_file_out)//'.nspec',status='unknown')
            if(field_evolves) then
#ifdef COMP
               if(comp_evolves) then
                  do i=1,Nr_s
                     write (unit_nspec,'(I4,4(D17.9,X))') i,tspec(i),cspec(i),uspec(i),bspec(i)
                  enddo
               else
#endif
                  do i=1,Nr_s
                     write (unit_nspec,'(I4,3(D17.9,X))') i,tspec(i),uspec(i),bspec(i)
                  enddo
#ifdef COMP
               endif
#endif
            else
#ifdef COMP
               if(comp_evolves) then
                  do i=1,Nr_s
                     write (unit_nspec,'(I4,3(D17.9,X))') i,tspec(i),cspec(i),uspec(i)
                  enddo
               else
#endif
                  do i=1,Nr_s
                     write (unit_nspec,'(I4,2(D17.9,X))') i,tspec(i),uspec(i)
                  enddo
#ifdef COMP
               endif  ! if(comp_evolves)
#endif
            endif ! if(field_evolves)
            close (unit_nspec)
         endif                  ! mpi_rank.eq.0
      endif
   end subroutine save_n_spec

   !-----------------------------------------------------------------------------
   !> Save the cfl numbers to file
   subroutine save_cfl()
      implicit none
      !-- save CFL-numbers:
      if(mpi_rank.eq.0) write (unit_cfl,'(5(D17.9,X))') time, cfl(1:4)
   end subroutine
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
