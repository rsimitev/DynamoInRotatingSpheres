module drs_io_DX
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
   double precision:: cut_phi !< the azimuth to use on meridional cuts
   double precision:: cut_z !< the azimuth to use on equator parallell cuts
   double precision:: where_to_cut = 0.0d0
   integer:: cut_type !< the type of cut or render to save

   interface save2DX
      module procedure save2DXscalar, save2DXvector
   end interface
contains

   !> Saves the contents of a scalar field  to file.
   !! @par is a field in real (tpr) space.
   !! @par filename is the base name for the output files.
   subroutine save2DXscalar(field, filename)
      implicit none
      double precision, intent(in):: field(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      character(len=*), intent(in):: filename
      select case(cut_type)
         case(1)
            call saveDXvolume(field, filename)
         case(2)
            cut_phi = where_to_cut
            call saveDXmeridional(field, filename)
      end select
   end subroutine save2DXscalar

   !> Saves the contents of a vector field to file given its three components. 
   !! @par XX, YY, ZZ are the real space components of the vector field
   !! @par filename is the base name for the output files.
   subroutine save2DXvector(XX, YY, ZZ, filename)
      implicit none
      double precision, intent(in):: XX(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: YY(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: ZZ(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      character(len=*), intent(in):: filename
      select case(cut_type)
         case(1)
            call saveDXvolume3DVec(XX,YY,ZZ, filename)
         case(2)
            cut_phi = where_to_cut
            call saveDXmeridional3DVec(XX,YY,ZZ, filename)
      end select
   end subroutine save2DXvector

   !> Writes a meridional slice of the field @par field
   !! into the files with basename @par filename.
   subroutine saveDXmeridional(field, filename)
      implicit none
      double precision, intent(in):: field(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      character(len=*), intent(in):: filename
      double precision:: render_out(0:(blk_t_size(mpi_rank)-1),Nr)
      double precision:: phi_n, dphi, dphi1, dphi2
      double precision:: mcut_min, mcut_max, mcut_rms
      integer:: i,l,icut

      dphi  = 2*pi/m0/Np
      icut  = int(cut_phi/dphi) + 1
      phi_n = icut*dphi
      dphi1 = cut_phi - phi_n
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

      open(unit=300, file=trim(filename)//'.general', status='UNKNOWN')

      ! Write the DX header
      ! write(300,*) 'file = ', trim(filename)//'.dat'
      write(300,"(' grid = ',I5,' x ',I5)") Nr, Nt+3
      write(300,*) 'format = ascii'
      write(300,*) 'interleaving = field'
      write(300,*) 'majority = column'
      write(300,*) 'field = locations, field0'
      write(300,*) 'structure = 2-vector, scalar'
      write(300,*) 'type = float, float'
      write(300,*) 'end'
      !close(300)
      ! And now the content
      ! open(unit=300, file=trim(filename)//'.dat', status='UNKNOWN')
      do i=1,Nr
         write(300,'(3F15.9)') 0.0d0, rcoll(i), render_out(Nt,i)
      enddo
      do l=Nt,0,-1
         do i=1,Nr
            write(300,'(3F15.9)') rcoll(i)*sqrt(1-costheta(l)**2), rcoll(i)*costheta(l), render_out(l,i)
         enddo
      enddo
      do i=1,Nr
         write(300,'(3F15.9)') 0.0d0, -rcoll(i), render_out(0,i)
      enddo

      close(300)
   end subroutine saveDXmeridional

   !> Writes a meridional slice of the field @par field
   !! into the files with basename @par filename.
   subroutine saveDXmeridional3DVec(field_x,field_y,field_z, filename)
      implicit none
      double precision, intent(in):: field_x(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: field_y(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: field_z(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      character(len=*), intent(in):: filename
      double precision:: render_x(0:(blk_t_size(mpi_rank)-1),Nr)
      double precision:: render_y(0:(blk_t_size(mpi_rank)-1),Nr)
      double precision:: render_z(0:(blk_t_size(mpi_rank)-1),Nr)
      double precision:: phi_n, dphi, dphi1, dphi2
      integer:: i,l,icut

      dphi  = 2*pi/m0/Np
      icut  = int(cut_phi/dphi) + 1
      phi_n = icut*dphi
      dphi1 = cut_phi - phi_n
      dphi2 = dphi - dphi1

      ! Interpolate in phi
      render_x = (field_x(:,icut,:)*dphi2 + field_x(:,icut+1,:)*dphi1)/dphi
      render_y = (field_y(:,icut,:)*dphi2 + field_y(:,icut+1,:)*dphi1)/dphi
      render_z = (field_z(:,icut,:)*dphi2 + field_z(:,icut+1,:)*dphi1)/dphi

      open(unit=300, file=trim(filename)//'-vec-merid.general', status='UNKNOWN')

      ! Write the DX header
      ! write(300,*) 'file = ', trim(filename)//'-vec-merid.dat'
      write(300,"(' grid = ',I5,' x ',I5)") Nr, Nt+3
      write(300,*) 'format = ascii'
      write(300,*) 'interleaving = field'
      write(300,*) 'majority = column'
      write(300,*) 'field = locations, field0'
      write(300,*) 'structure = 2-vector, 3-vector'
      write(300,*) 'type = float, float'
      write(300,*) 'end'
      ! close(300)
      ! And now the content
      ! open(unit=300, file=trim(filename)//'-vec-merid.dat', status='UNKNOWN')
      do i=1,Nr
         write(300,'(5F15.9)') 0.0d0, rcoll(i), render_x(Nt,i), render_y(Nt,i), render_z(Nt,i)
      enddo
      do l=Nt,0,-1
         do i=1,Nr
            write(300,'(5F15.9)') rcoll(i)*sqrt(1-costheta(l)**2), rcoll(i)*costheta(l), render_x(l,i), render_y(l,i), render_z(l,i)
         enddo
      enddo
      do i=1,Nr
         write(300,'(5F15.9)') 0.0d0, -rcoll(i), render_x(0,i), render_y(0,i), render_z(0,i)
      enddo

      close(300)
   end subroutine saveDXmeridional3DVec

   !> Writes a volume rendeer of the field @par field
   !! into the files with basename @par filename.
   subroutine saveDXvolume(field, filename)
      implicit none
      double precision, intent(in):: field(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      character(len=*), intent(in):: filename
      double precision:: x, y, z, sinp, cosp, dphi
      double precision:: vol_min, vol_max
      integer:: i,l, j

      dphi  = 2*pi/m0/Np

      ! find minimum and maximum
      vol_min = minval(field)
      vol_max = maxval(field)

      open(unit=300, file=trim(filename)//'-vol.general', status='UNKNOWN')

      ! Write the DX header
      ! write(300,*) 'file = ',trim(filename)//'-vol.dat'
      write(300,"(' grid = ',I5,' x ',I5,' x ',I5)") Nr, Nt+3, Np+1
      write(300,*) 'format = ascii'
      write(300,*) 'interleaving = field'
      write(300,*) 'majority = column'
      write(300,*) 'field = locations, field0'
      write(300,*) 'structure = 3-vector, scalar'
      write(300,*) 'type = float, float'
      write(300,*) 'end'
      ! close(300)
      ! And now the content
      ! open(unit=300, file=trim(filename)//'-vol.dat', status='UNKNOWN')
      do j=1, Np
         sinp = sin((j-1)*dphi)
         cosp = cos((j-1)*dphi)
         ! The North pole
         do i=1,Nr
            !write(300,'(3F15.9)') x,     y,     z       , interpolated value
            write(300,'(4F15.9)')  0.0d0, 0.0d0, rcoll(i), sum(field(Nt,:,i))/Np
         enddo
         do l=Nt,0,-1
            do i=1,Nr
               x = rcoll(i)*sqrt(1-costheta(l)**2)*cosp
               y = rcoll(i)*sqrt(1-costheta(l)**2)*sinp
               z = rcoll(i)*costheta(l)
               write(300,'(4F15.9)') x, y, z, field(l,j,i)
            enddo
         enddo
         ! The South pole
         do i=1,Nr
            write(300,'(4F15.9)') 0.0d0, 0.0d0,  -rcoll(i), sum(field(0,:,i))/Np
         enddo
      enddo
      ! Duplicate Greenwich to get continuity in openDX
      sinp = 0.0d0
      cosp = 1.0d0
      ! The North pole
      do i=1,Nr
         !write(300,'(3F15.9)') x,     y,     z       , interpolated value
         write(300,'(4F15.9)')  0.0d0, 0.0d0, rcoll(i), sum(field(Nt,:,i))/Np
      enddo
      do l=Nt,0,-1
         do i=1,Nr
            x = rcoll(i)*sqrt(1-costheta(l)**2)*cosp
            y = rcoll(i)*sqrt(1-costheta(l)**2)*sinp
            z = rcoll(i)*costheta(l)
            write(300,'(4F15.9)') x, y, z, field(l,1,i)
         enddo
      enddo
      ! The South pole
      do i=1,Nr
         write(300,'(4F15.9)') 0.0d0, 0.0d0,  -rcoll(i), sum(field(0,:,i))/Np
      enddo

      close(300)
   end subroutine saveDXvolume

   !> Writes a volume rendeer of the field @par field
   !! into the files with basename @par filename.
   subroutine saveDXvolume_v2(field, filename)
      implicit none
      double precision, intent(in):: field(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      character(len=*), intent(in):: filename
      double precision:: x, y, z, sinp, cosp, dphi
      double precision:: vol_min, vol_max
      integer:: i,l, j

      dphi  = 2*pi/m0/Np

      ! find minimum and maximum
      vol_min = minval(field)
      vol_max = maxval(field)

      open(unit=300, file=trim(filename)//'-vol.general', status='UNKNOWN')

      ! Write the DX header
      ! write(300,*) 'file = ',trim(filename)//'-vol.dat'
      write(300,*) 'points = ', Nr*((Nt+1)*Np + 2)
      write(300,*) 'format = ascii'
      write(300,*) 'interleaving = field'
      write(300,*) 'majority = column'
      write(300,*) 'field = locations, field0'
      write(300,*) 'structure = 3-vector, scalar'
      write(300,*) 'type = float, float'
      write(300,*) 'end'
      ! close(300)
      ! And now the content
      ! open(unit=300, file=trim(filename)//'-vol.dat', status='UNKNOWN')
      ! The North pole
      do i=1,Nr
         !write(300,'(3F15.9)') x,     y,     z       , interpolated value
         write(300,'(4F15.9)')  0.0d0, 0.0d0, rcoll(i), sum(field(Nt,:,i))/Np
      enddo
      do j=1, Np
         sinp = sin((j-1)*dphi)
         cosp = cos((j-1)*dphi)
         do l=Nt,0,-1
            do i=1,Nr
               x = rcoll(i)*sqrt(1-costheta(l)**2)*cosp
               y = rcoll(i)*sqrt(1-costheta(l)**2)*sinp
               z = rcoll(i)*costheta(l)
               write(300,'(4F15.9)') x, y, z, field(l,j,i)
            enddo
         enddo
      enddo
      ! The South pole
      do i=1,Nr
         !write(300,'(3F15.9)') x,     y,     z       , interpolated value
         write(300,'(4F15.9)') 0.0d0, 0.0d0,  -rcoll(i), sum(field(0,:,i))/Np
      enddo
      close(300)
   end subroutine saveDXvolume_v2

   !> Writes a volume rendeer of the vector field components @par XX, YY and ZZ
   !! into the files with basename @par filename.
   subroutine saveDXvolume3DVec(XX, YY, ZZ, filename)
      implicit none
      double precision, intent(in):: XX(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: YY(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      double precision, intent(in):: ZZ(0:(blk_t_size(mpi_rank)-1),Np,Nr)
      character(len=*), intent(in):: filename
      double precision:: x, y, z, sinp, cosp, dphi
      integer:: i,l, j

      dphi  = 2*pi/m0/Np

      open(unit=300, file=trim(filename)//'-vec-vol.general', status='UNKNOWN')

      ! Write the DX header
      ! write(300,*) 'file = ',trim(filename)//'-vec-vol.dat'
      write(300,"(' grid = ',I5,' x ',I5,' x ',I5)") Nr, Nt+3, Np+1
      write(300,*) 'format = ascii'
      write(300,*) 'interleaving = field'
      write(300,*) 'majority = column'
      write(300,*) 'field = locations, field0'
      write(300,*) 'structure = 3-vector, 3-vector'
      write(300,*) 'type = float, float'
      write(300,*) 'end'
      ! close(300)
      ! And now the content
      ! open(unit=300, file=trim(filename)//'-vec-vol.dat', status='UNKNOWN')
      do j=1, Np
         sinp = sin((j-1)*dphi)
         cosp = cos((j-1)*dphi)
         ! The North pole
         do i=1,Nr
            !write(300,'(3F15.9)') x,     y,     z       , interpolated value
            write(300,'(6F15.9)')  0.0d0, 0.0d0, rcoll(i), sum(XX(Nt,:,i))/Np, sum(YY(Nt,:,i))/Np, sum(ZZ(Nt,:,i))/Np
         enddo
         do l=Nt,0,-1
            do i=1,Nr
               x = rcoll(i)*sqrt(1-costheta(l)**2)*cosp
               y = rcoll(i)*sqrt(1-costheta(l)**2)*sinp
               z = rcoll(i)*costheta(l)
               write(300,'(6F15.9)') x, y, z, XX(l,j,i), YY(l,j,i), ZZ(l,j,i)
            enddo
         enddo
         ! The South pole
         do i=1,Nr
            write(300,'(6F15.9)') 0.0d0, 0.0d0,  -rcoll(i), sum(XX(0,:,i))/Np, sum(YY(0,:,i))/Np, sum(ZZ(0,:,i))/Np
         enddo
      enddo
      ! Duplicate Greenwich to get continuity in openDX
      sinp = 0.0d0
      cosp = 1.0d0
      ! The North pole
      do i=1,Nr
         !write(300,'(3F15.9)') x,     y,     z       , interpolated value
         write(300,'(6F15.9)')  0.0d0, 0.0d0, rcoll(i), sum(XX(Nt,:,i))/Np, sum(YY(Nt,:,i))/Np, sum(ZZ(Nt,:,i))/Np
      enddo
      do l=Nt,0,-1
         do i=1,Nr
            x = rcoll(i)*sqrt(1-costheta(l)**2)*cosp
            y = rcoll(i)*sqrt(1-costheta(l)**2)*sinp
            z = rcoll(i)*costheta(l)
            write(300,'(6F15.9)') x, y, z, XX(l,1,i), YY(l,1,i), ZZ(l,1,i)
         enddo
      enddo
      ! The South pole
      do i=1,Nr
         write(300,'(6F15.9)') 0.0d0, 0.0d0,  -rcoll(i), sum(XX(0,:,i))/Np, sum(YY(0,:,i))/Np, sum(ZZ(0,:,i))/Np
      enddo

      close(300)
   end subroutine saveDXvolume3DVec

end module drs_io_DX

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
