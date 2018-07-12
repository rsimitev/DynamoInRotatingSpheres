! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2014
!> Deals with input and output of state files and .par files.
module drs_io_state
#include "drsDefs.F90"
   use drs_params
   use drs_dims
   use drs_time
   use drs_legendre
   use drs_fftw3
   use drs_mpi
   use drs_flow
   use drs_field
   use drs_temp
#ifdef COMP
   use drs_comp
#endif
   use drs_error_codes
   implicit none
   save
   private
   integer:: input_data_file_ver
   integer, parameter:: unit_par=60
   !> These are the output_data_file_ver numbers we are going to write
#ifdef COMP
   integer, public, parameter:: output_data_file_ver = 10209
#else
   integer, public, parameter:: output_data_file_ver = 10205
#endif
   !> These are the output_data_file_ver numbers we know about.
   !! They correspond to the versions of the par files we can read.
   !! Add 100 for the new normalisation of the SH.
   integer, parameter:: DATA_FILE_VER_1 = 10101
   integer, parameter:: DATA_FILE_VER_2 = 10102
   integer, parameter:: DATA_FILE_VER_3 = 10103
   integer, parameter:: DATA_FILE_VER_4 = 10104
   integer, parameter:: DATA_FILE_VER_5 = 10105
   integer, parameter:: DATA_FILE_VER_6 = 10106
   integer, parameter:: DATA_FILE_VER_7 = 10107
   integer, parameter:: DATA_FILE_VER_9 = 10109
   ! MPI
   integer, target:: usr_dimsi(8)
   ! Public variables
   integer, public:: lformi, drs_calc_typei, tempBCi, flowBCi, magBCi
   integer, public:: compBCi
   integer, public:: Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi, m0i
   integer, public:: stepmaxi, sample_ratei, transienti
   double precision, public:: etai,Pti,Tai,Ra_ti,Pmi,hi
   double precision, public:: Ra_ci !< Compositional Rayleigh of input
   double precision, public:: Pci   !< Compositional Prandtl of input
   double precision, public:: drifti
   character(len=60), public:: commenti

   character(len=13), parameter, public:: deflate_state = 'compressdata '
   character(len=15), parameter, public:: inflate_state = 'uncompressdata '
   public:: drs_load_state, drs_save_state, drs_read_state_par, drs_write_state_par

contains

   !---------------------------------------------------------------
   !>     reads the parameterfile 'file'.par
   subroutine drs_read_state_par(unit_in)
      implicit none
      integer, intent(in):: unit_in
      integer:: error, iaux
      error = 0
      read(unit_in,*)
      ! Read the input_data_file_ver number so we can make an
      ! informed decision on how to read the rest of this file.
      read(unit_in,*) input_data_file_ver
      backspace(unit_in)
      select case(input_data_file_ver)
         ! output_data_file_verC+100 are the new normalizations
         case(DATA_FILE_VER_1, DATA_FILE_VER_1+100)
            read(unit_in,*) iaux, lformi, drs_calc_typei, tempBCi, flowBCi, magBCi
            read(unit_in,*)
            read(unit_in,*) etai,Pti,Tai,Ra_ti,Pmi
            Ra_ci   = 0.0d0
            Pci     = 0.0d0
            compBCi = 0
            read(unit_in,*)
            read(unit_in,*) Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi
            m0i = 1
            read(unit_in,*)
            read(unit_in,*) hi, stepmaxi, transienti, sample_ratei
            stepstart = 0
            time   = 0.0d0
            drifti = 0.0d0
         case(DATA_FILE_VER_2, DATA_FILE_VER_2+100)
            read(unit_in,*) iaux, lformi, drs_calc_typei, tempBCi, flowBCi, magBCi
            read(unit_in,*)
            read(unit_in,*) etai,Pti,Tai,Ra_ti,Pmi
            Ra_ci   = 0.0d0
            Pci     = 0.0d0
            compBCi = 0
            read(unit_in,*)
            read(unit_in,*) Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi
            m0i = 1
            read(unit_in,*)
            read(unit_in,*) hi, stepmaxi, transienti, sample_ratei, stepstart
            time   = hi*stepstart
            drifti = 0.0d0
         case(DATA_FILE_VER_3, DATA_FILE_VER_3+100)
            read(unit_in,*) iaux, lformi, drs_calc_typei, tempBCi, flowBCi, magBCi
            read(unit_in,*)
            read(unit_in,*) etai,Pti,Tai,Ra_ti,Pmi
            Ra_ci   = 0.0d0
            Pci     = 0.0d0
            compBCi = 0
            read(unit_in,*)
            read(unit_in,*) Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi
            m0i = 1
            read(unit_in,*)
            read(unit_in,*) hi, stepmaxi, transienti, sample_ratei, stepstart, time, drifti
         case(DATA_FILE_VER_4,DATA_FILE_VER_4+100,DATA_FILE_VER_5,DATA_FILE_VER_5+100)
            read(unit_in,*) iaux, lformi, drs_calc_typei, tempBCi, flowBCi, magBCi
            read(unit_in,*)
            read(unit_in,*) etai,Pti,Tai,Ra_ti,Pmi
            Ra_ci   = 0.0d0
            Pci     = 0.0d0
            compBCi = 0
            read(unit_in,*)
            read(unit_in,*) Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi, m0i
            read(unit_in,*)
            read(unit_in,*) hi, stepmaxi, transienti, sample_ratei, stepstart, time, drifti
         case(DATA_FILE_VER_7,DATA_FILE_VER_7+100,DATA_FILE_VER_9,DATA_FILE_VER_9+100)
            read(unit_in,*) iaux,lformi,drs_calc_typei,tempBCi,compBCi,flowBCi,magBCi
            read(unit_in,*)
            read(unit_in,*) etai,Pti,Pci,Tai,Ra_ti,Ra_ci,Pmi
            read(unit_in,*)
            read(unit_in,*) Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi, m0i
            read(unit_in,*)
            read(unit_in,*) hi, stepmaxi, transienti, sample_ratei, stepstart, time, drifti
         case default
            write(*,*) 'Unknown output_data_file_ver number in input file: ',input_data_file_ver
            stop 1
      end select

      read(unit_in,*) commenti
      Tai = sqrt(Tai)
   end subroutine drs_read_state_par

   !---------------------------------------------------------------
   !> Writes the parameter file 'file'.par
   subroutine drs_write_state_par(unit_out)
      implicit none
      integer, intent(in):: unit_out
      double precision:: TaOut
      TaOut=Ta*Ta
#ifdef COMP
      write(unit_out,'(A)') '| output_data_file_ver | lformat | drs_calc_type | tempBC | compBC | flowBC | magBC |'
      write(unit_out,'(7I8)') output_data_file_ver,lform,drs_calc_type,tempBC(1),compBC(1),flowBC(1),magBC(1)
#else
      write(unit_out,'(A)') '| output_data_file_ver | lformat | drs_calc_type | tempBC | flowBC | magBC |'
      write(unit_out,'(6I8)') output_data_file_ver,lform,drs_calc_type,tempBC(1),flowBC(1),magBC(1)
#endif

#ifdef COMP
      write(unit_out,'(A)') '|  eta  | Therm. Prandtl | Compo. Prandtl | Taylor | Therm. Rayleigh | Compo. Rayleigh |  Magn. Prandtl  |'
      write(unit_out,'(7D14.5)') eta, Pt, Pc, TaOut, Ra_t, Ra_c, Pm
#else
      write(unit_out,'(A)') '|  eta   | Prandtl |   Taylor   |  Rayleigh   |  magPr  |'
      write(unit_out,'(5D14.5)') eta, Pt, TaOut, Ra_t, Pm
#endif
      write(unit_out,'(A)') '| Nr | Nt | Np | Nr_s | Nt_s |Np_s | lsymm | m0 |'
      write(unit_out,'(8I8)') Nr, Nt, Np, Nr_s, Nt_s, Np_s, lsymm, m0
      write(unit_out,'(A)') '|delta_t|nsteps|transient|sample_rate|step|time|drift|'
      write(unit_out,'(D11.3,I7,I7,I6,I8,D14.5,D11.3)') &
      h, stepmax, transient, sample_rate_sim_time, steps, time, drift
      write(unit_out,*) '''',comment,''''

   end subroutine drs_write_state_par


   !---------------------------------------------------------------
   !> Reads a state performing interpolation as needed.
   !! The state is stored in the files with name
   !! given by \a io_calc_file_in and are described by the file with extension
   !! .par.
   subroutine drs_load_state(error)
      implicit none
      integer, intent(out):: error
      logical:: par_file_exists

      error=0
      if(mpi_rank.eq.0) then
         inquire(file=trim(io_calc_file_in)//'.par',exist=par_file_exists)
         if(par_file_exists) then
            !-- read parameter of input data:
            open(unit_par,file=trim(io_calc_file_in)//'.par',form='formatted',status='old')
            call drs_read_state_par(unit_par)
            close(unit_par)
            call check_input_par(error)

            steps = stepstart
            usr_dimsi = (/ Nri, Nti, Npi, Nri_s, Nti_s, Npi_s, lsymmi, m0i/)
         else
            spew 'No parameter input file. ABORTING!'
            error = ERR_NO_PAR_FILE
         endif
      endif
      call drs_bcast(error)
      if(error>0) return

      call drs_bcast(usr_dimsi, size(usr_dimsi, 1))
      Nri    = usr_dimsi(1)
      Nti    = usr_dimsi(2)
      Npi    = usr_dimsi(3)
      Nri_s  = usr_dimsi(4)
      Nti_s  = usr_dimsi(5)
      Npi_s  = usr_dimsi(6)
      lsymmi = usr_dimsi(7)
      m0i    = usr_dimsi(8)
      call drs_bcast(lformi)

      ! This is the model time of the input state
      time_start = time

      ! Allocate the flow, just in case.
      ! This will do nothing if the flow was already allocated.
      call drs_flow_allocation()
      !-- read input data:
      call read_3Dfieldp(trim(io_calc_file_in)//'.pol',  flow_pol, error)
      if (error==1) then
         spew 'No poloidal flow input data. ABORTING!'
         error = ERR_NO_FLOW_POL_FILE
         return
      endif
      if(mpi_rank.eq.0) flow_pol(0,1,1:Nr) = 0.d0

      call read_3Dfieldp(trim(io_calc_file_in)//'.tor',  flow_tor, error)
      if (error==1) then
         spew 'No toroidal flow input data. ABORTING!'
         error = ERR_NO_FLOW_TOR_FILE
         return
      endif
      if(mpi_rank.eq.0) flow_tor(0,1,1:Nr) = 0.d0

      if(temp_present) then
         ! Allocate the temperature, just in case.
         ! This will do nothing if the temperature was already allocated.
         call drs_temp_allocation()
         call read_3Dfieldp(trim(io_calc_file_in)//'.temp', temp, error)
         if (error==1) then
            spew 'No temperature input data. ABORTING!'
            error =  ERR_NO_TEMPERATURE_FILE
            return
         endif
      endif

      if (field_present) then
         ! Allocate the field, just in case.
         ! This will do nothing if the field was already allocated.
         call drs_field_allocation()
         call read_3Dfieldp(trim(io_calc_file_in)//'.Bp', field_pol, error)
         if (error==1) then
            spew 'No poloidal field input data. ABORTING!'
            error = ERR_NO_MAG_POL_FILE
            return
         endif
         if(mpi_rank.eq.0) field_pol(0,1,1:Nr) = 0.d0

         call read_3Dfieldp(trim(io_calc_file_in)//'.Bt', field_tor, error)
         if (error==1) then
            spew 'No toroidal field input data. ABORTING!'
            error = ERR_NO_MAG_TOR_FILE
            return
         endif
         if(mpi_rank.eq.0) field_tor(0,1,1:Nr) = 0.d0
      endif
#ifdef COMP
      if(comp_present) then
         ! Allocate the composition, just in case.
         ! This will do nothing if the composition was already allocated.
         call drs_comp_allocation()
         call read_3Dfieldp(trim(io_calc_file_in)//'.thetac', comp, error)
         if (error==1) then
            spew 'No composition input data. ABORTING!'
            error = ERR_NO_COMPOSITION_FILE
            return
         endif
      endif
#endif
   end subroutine drs_load_state

   !---------------------------------------------------
   !> Subroutine to deal with a change of axial symmetry.\n
   !! Applyed when m0 on input is not the same as the m0 for the present
   !! calculation.
   subroutine axial_symmetry_transform(m0i, m0, N, t)
      implicit none
      integer, intent(in):: m0, m0i, N
      double precision, intent(inout):: t(Nri, 0:Nti_s, N)
      double precision:: t1(Nri, 0:Nti_s, N)
      integer:: j, m
      integer:: ji, joff

      !-- if the input has another wavenumber than m0, we have to push
      !-- the coefficients to their correct positions:
      if(m0.eq.m0i) return
      t1 = 0.0d0
      if(m0.gt.m0i) then
         !-- increasing m0: picking up the right data
         do j=1, N
            m    = (j/2)*m0 ! Wavenumbers I want contiguous in j
            joff = mod(j,2) ! Real or imaginary
            ji   = 2*(m/m0i) + joff ! Index of the wavenumbers I want in the old matrix
            if(ji.le.N) then
               t1(1:Nri,0:Nti_s,j) = t(1:Nri,0:Nti_s,ji)
            endif
         enddo
      elseif(m0.lt.m0i) then
         !-- reducing m0: expanding the data
         do ji=1, N
            m    = (ji/2)*m0i ! Wavenumbers that were contiguous in the old matrix
            joff = mod(ji,2)  ! Real or imaginary
            j    = 2*(m/m0) + joff ! Index of the wavenumbers I want in the new matrix
            if(j.le.N) then
               t1(1:Nri,0:Nti_s,j) = t(1:Nri,0:Nti_s,ji)
            endif
         enddo
      endif  ! m0i.gt.m0

      !-- copy data back to the input fields
      t = 0.0d0
      do j=1,N
         m = m0*(j/2)
         t(1:Nri,m:Nti_s,j) = t1(1:Nri,m:Nti_s,j)
      enddo

      !-- transformation of wavenumber m0 done.
      !-- the mode m=mmax has to be zero if N is even (only real part)
      if((mod(N,2).eq.0)) then
         t(1:Nri,0:Nti_s,N) = 0.D0
      endif  !  mod(N,2).eq.0
   end subroutine axial_symmetry_transform

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

   !****************************************************************
   ! subroutine write_3Dfieldp(file,t,buf,usr_dims,mpi_dims,lform)
   !
   !> write one field t to disk.
   !
   ! usr_dims(*),mpi_dims(*) are global and local dimensions
   ! lform             format of datafiles: 0=unformatted, 1=formatted
   !! \todo Implement netcdf output.
   !****************************************************************
   subroutine write_3Dfieldp(file,t)
      implicit none
      character*(*), intent(in):: file
      double precision, intent(in):: t(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      double precision, allocatable:: cache(:,:,:)
      integer:: nadl
      integer:: j,l,i

      nadl = blk_ps_size(mpi_rank)

      if(mpi_rank.eq.0) then !-- write own data:
         allocate(cache(0:Nt_s, Np_s, Nr))
         cache = 0.0d0
      else
         allocate(cache(0:Nt_s, blk_ps_size(mpi_rank), Nr))
      endif
      cache(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr) = t(0:Nt_s,1:blk_ps_size(mpi_rank),1:Nr)
      call gather_from_m(cache, Nt_s, Nr)
      if(mpi_rank.eq.0) then !-- write data:
         if(lform.eq.0) then
            open(420,FILE=trim(file),form='unformatted',status='unknown')
         else
            open(420,FILE=trim(file),status='unknown')
         endif
         do j=1, Np_s
            do l=m0*(j/2), Nt_s
               do i=1, Nr
                  if(lform.eq.0) then
                     write(420) cache(l,j,i)
                  else
                     write (420,'(D19.12)') cache(l,j,i)
                  endif
               enddo
            enddo
         enddo
         close(420)
      endif
      if(allocated(cache)) deallocate(cache)
   end subroutine write_3Dfieldp


   !****************************************************************
   ! subroutine read_3Dfieldp(file,t,buf,usr_dims,mpi_dims,lform)
   !
   ! read one field t from disk in lmr space.
   !
   ! usr_dims(*),mpi_dims(*) are global and local dimensions
   ! lform             format of datafiles: 0=unformatted, 1=formatted
   !****************************************************************
   subroutine read_3Dfieldp(file,t,error)
      implicit none

      character*(*), intent(in):: file
      double precision, intent(out):: t(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      integer, intent(out):: error
      double precision, allocatable:: cache(:,:,:)
      double precision, allocatable:: oldNorm(:,:), newNorm(:,:)
      double precision:: aux1(Nri), aux2(Nr)
      integer:: i,j,l, N, m
      integer:: blk_size,blk_start
      logical:: file_exists

      error = 0
      blk_size  = blk_ps_size(mpi_rank)
      blk_start = blk_ps_start(mpi_rank)

      ! Check whether the file we are about to open exists
      if(mpi_rank.eq.0) then
         inquire(file=trim(file),exist=file_exists)
      endif
      ! Tell everyone and return in case it doesn't.
      call drs_bcast(file_exists)
      if(.not.file_exists) then
         error = 1
         return
      endif

      N = max(Npi_s,Np_s)
      if(mpi_rank.eq.0) then
         allocate(cache(Nri, 0:Nti_s, N))
      else
         allocate(cache(Nri, 0:Nti_s, blk_ps_size(mpi_rank)))
      endif
      cache = 0.0d0
      t     = 0.0d0
      if(mpi_rank.eq.0) then
         !-- read own data:
         if(lformi.eq.0) then
            open (410,FILE=trim(file),form='unformatted',status='OLD')
         else
            open (410,FILE=trim(file),form='formatted',status='OLD')
         endif
         !-- Npi_s is dim of input file, blk_ps_size(mpi_rank) is dim for this calculation!
         ! If the input state contains unnormalized coefficients
         ! apply a norm transformation
         if (input_data_file_ver < 10200) then
            allocate( oldNorm(0:Nti_s+2, 0:Nti_s+2) )
            allocate( newNorm(0:Nti_s+2, 0:Nti_s+2) )
            call initNormalization(UnNormalized, Nti_s, oldNorm)
            call initNormalization(Normalised, Nti_s, newNorm)
         endif
         do j=1, Npi_s
            m = m0i*(j/2)
            do l = m, Nti_s
               do i=1, Nri
                  if(lformi.eq.0) then
                     read (410) cache(i,l,j)
                  else
                     read (410,*) cache(i,l,j)
                  endif
                  ! If the input state contains unnormalized coefficients
                  ! apply a norm transformation
                  if (input_data_file_ver < 10200) then
                     cache(i,l,j) = cache(i,l,j)*(oldNorm(l,m)/newNorm(l,m))
                  endif

               enddo
            enddo
         enddo
         if (input_data_file_ver < 10200) then
            deallocate(newNorm)
            deallocate(oldNorm)
         endif
         close(410)
         call axial_symmetry_transform(m0i, m0, N, cache)
      endif
      call distribute_in_m(cache, Nti_s, Nri)
      t = 0.0d0

      do j=1, blk_ps_size(mpi_rank)
         do l = m0*(j/2), min(Nti_s, Nt_s)
            aux1(1:Nri) = cache(1:Nri,l,j)
            call remesh(Nri, aux1, Nr, aux2)
            t(l,j,1:Nr) = aux2(1:Nr)
         enddo
      enddo
      deallocate(cache)
   end subroutine read_3Dfieldp

   !> Saves the present state to file.
   !! At this point all files are saved with the file name
   !! given by \a io_calc_file_out
   subroutine drs_save_state()
      implicit none
      !-- write parameterfile .par:
      if(mpi_rank.eq.0) then
         open(unit_par, file=trim(io_calc_file_out)//'.par',  form='formatted',   status='unknown')
         call drs_write_state_par(unit_par)
         close(unit_par)
      endif
      !-- write the fields flow_pol,flow_tor,temp,field_pol,field_tor:
      if(flow_evolves) then
         call write_3Dfieldp(trim(io_calc_file_out)//'.pol',  flow_pol)
         call write_3Dfieldp(trim(io_calc_file_out)//'.tor',  flow_tor)
      endif
      if(temp_evolves) then
         call write_3Dfieldp(trim(io_calc_file_out)//'.temp', temp)
      endif
#ifdef COMP
      if(comp_evolves) then
         call write_3Dfieldp(trim(io_calc_file_out)//'.thetac', comp)
      endif
#endif
      if(field_evolves) then
         call write_3Dfieldp(trim(io_calc_file_out)//'.Bp', field_pol)
         call write_3Dfieldp(trim(io_calc_file_out)//'.Bt', field_tor)
      endif
   end subroutine drs_save_state
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
