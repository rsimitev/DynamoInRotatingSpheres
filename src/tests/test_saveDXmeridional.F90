! Copyright  L. Silva (lacsilva@gmail.com), 2015
#include "drsDefs.f90"
! ----------------
!> 
! ----------------
program test_saveDXMer
   use drs_time
   use drs_params
   use drs_dims
   use drs_mpi
   use drs_fftw3
   use drs_legendre
   use drs_radial
   implicit none
   integer:: error
   integer:: l, j, i
   double precision, allocatable:: testField(:,:,:)

   Nr = 30
   Nt = 180
   Np = 361
   Nr_s = 30
   Nt_s = 180
   Np_s = 361
   m0  = 1
   eta = 0.35
   call drs_mpi_init(error)
   if(error.ne.0) stop
   if(mpi_size.ne.1) then
    spew 'This program should be ran on a single cpu'
    call drs_abort(1)
   endif
   call drs_params_init()
   call drs_dims_init(error)
   call mpi_dims_init(Nt, Np_s, m0, error)

   call drs_fftw3_init(Nr, blk_t_size(mpi_rank), Np)
   call drs_legendre_allocation()
   call drs_legendre_init()
   call drs_radial_init(eta)

   allocate(testField(0:(blk_t_size(mpi_rank)-1),Np,Nr))

   do j=1, Np
      do l=0, Nt
         do i=1, Nr
            testField(l,j,i) = dble(100**3*i + 100*l + j)
         enddo
      enddo
   enddo

   call saveDXmeridional(testField, pi*0.5d0, 'test_saveDXmeridional')
contains
   subroutine saveDXmeridional(field, phi, filename)
      implicit none
      double precision, intent(in):: field(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: phi
      character(len=*), intent(in):: filename
      double precision:: render_out(0:(blk_t_size(mpi_rank)-1),Nr)
      double precision:: phi_n, dphi, dphi1, dphi2
      double precision:: mcut_min,mcut_max
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

      open(unit=300, file=filename//'.general', status='UNKNOWN')

      ! Write the DX header
      write(300,*) 'file = ', filename//'.dat'
      write(300,"(' grid = ',I5,' x ',I5)") Nr, Nt+1
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
      do l=Nt,0,-1
         do i=1,Nr
            write(300,'(3E15.5)') rcoll(i)*sqrt(1-costheta(l)**2), rcoll(i)*costheta(l), render_out(l,i)
         enddo
      enddo

      close(300)

   end subroutine
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
