! Copyright  L. Silva (lacsilva@gmail.com), 2015
#include "drsDefs.f90"
program test_vectorField2Divergence
   use drs_dims
   use drs_params
   use drs_mpi
   use drs_fftw3
   use drs_radial
   use drs_legendre
   use drs_transforms
   implicit none
   integer:: l, j, i
   double precision, allocatable:: ur(:,:,:)
   double precision, allocatable:: ut(:,:,:)
   double precision, allocatable:: up(:,:,:)
   double precision, allocatable:: div(:,:,:)
   logical:: pass

   call init()
   ut=0.0d0
   up=0.0d0
   ur=0.0d0
   div=0.0d0
   do i=1, Nr
      ur(:,:,i) = rcoll(i)
   enddo
   call vectorField2Divergence(ur, ut, up, div)
   do i=1, Nr
      jl_do(j,l) 
         if((l==0).and.(j==1)) then
            if(dabs(div(l,j,i)-dble(3))>1.0d-10) then
              pass=.false.
              Write(*,*) l,j,i, div(l,j,i)
            endif
         elseif(dabs(div(l,j,i))>1.0d-10) then
            pass=.false.
            Write(*,*) l,j,i, div(l,j,i)
         endif
      jl_enddo
   enddo
   if (pass) Write(*,*) 'Pass.'
contains
   subroutine init()
      implicit none
      integer:: error
      error=0
      call drs_mpi_init(error)
      if(error.ne.0) stop
      Nr = 11
      Nt = 13
      Np = 27
      Nr_s = 11
      Nt_s = 13
      Np_s = 27
      lsymm = 1
      m0 = 1
      eta = 0.5
      pass = .true.
      call drs_params_init()
      call drs_dims_init(error)
      call mpi_dims_init(Nt, Np_s, m0, error)
      ! Initialize the dft routines.
      call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)

      call drs_legendre_allocation()
      call drs_legendre_init()
      call drs_radial_init(eta)
      allocate(ur(0:blk_t_size(mpi_rank)-1,Np,Nr))
      allocate(ut(0:blk_t_size(mpi_rank)-1,Np,Nr))
      allocate(up(0:blk_t_size(mpi_rank)-1,Np,Nr))
      allocate(div(0:blk_t_size(mpi_rank)-1,Np,Nr))
   end subroutine init      
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
