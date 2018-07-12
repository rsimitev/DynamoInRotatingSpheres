! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Module with helper subroutines for debug
#include "drsDefs.F90"
module drs_debug
   use drs_dims
   use drs_mpi

contains
   ! Saves a quantity in lmr space to disk
   subroutine save_lmr_quantity(t, tname)
      implicit none
      !> Fragment of the array to be saved held by CPU with mpi_rank.
      double precision, intent(in):: t(0:Nt_s,blk_ps_size(mpi_rank),Nr)
      !> The name of the file that will be saved to disk.
      character(*),intent(in):: tname
      integer:: rank, l, j, i
      do rank=0, mpi_size-1
         if (mpi_rank==rank) then
            if (rank==0) then
               open(unit=472,file=tname,status='UNKNOWN')
            else
               open(unit=472,file=tname,status='UNKNOWN', access='APPEND')
            endif
            jl_do(j,l)
               do i=1, Nr
                  Write(472,*) l, j+blk_ps_start(rank)-1, i, t(l,j,i)
               enddo
            jl_enddo
            close(unit=472)
         endif
         call wait_for_everyone()
      enddo
   end subroutine save_lmr_quantity

   subroutine save_tpr_quantity(t, tname)
      implicit none
      double precision, intent(in):: t(0:blk_t_size(mpi_rank)-1,Np,Nr)
      character(*),intent(in):: tname
      integer:: rank, l, j, i
      do rank=0, mpi_size-1
         if (mpi_rank==rank) then
            if (rank==0) then
               open(unit=472,file=tname,status='UNKNOWN')
            else
               open(unit=472,file=tname,status='UNKNOWN', access='APPEND')
            endif
            do l=0, blk_t_size(rank)-1
               do j=1, Np
                  do i=1, Nr
                     Write(472,*) l+blk_t_start(rank), j, i, t(l,j,i)
                  enddo
               enddo
            enddo
            close(unit=472)
         endif
         call wait_for_everyone()
      enddo
   end subroutine save_tpr_quantity
end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
