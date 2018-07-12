! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Provides initialisation and variables to be used with the mpi implementation
module drs_mpi
#include "drsDefs.F90"
   use drs_error_codes
   implicit none
   save
#ifdef MPI
#include <mpif.h>
   integer:: mpi_err
   integer:: mpi_stat(MPI_STATUS_SIZE)
   double precision, allocatable:: mpibuffer(:)
#endif
   integer:: mpi_size !< How many CPU's are in use.
   integer:: mpi_rank !< The rank of the present CPU.
   !-- arrays that hold local dims of all procs:
   integer, target, allocatable:: blk_ps_start(:) !< Start index of the blocks in m for each CPU.
   integer, target, allocatable:: blk_ps_size(:)  !< Size of the blocks in m for each CPU.
   integer, target, allocatable:: blk_t_start(:)  !< Start index of the blocks in theta dor each CPU.
   integer, target, allocatable:: blk_t_size(:)   !< Size of the blocks in theta for each CPU.
   !> Stores the index of the first nonzero l value in the block.
   integer, target, allocatable:: blk_ts_start(:)
   integer, pointer:: mm(:) !< A convinience shorthand for blk_ts_start.
   integer:: blk_t_max_size !< Maximum size of theta block per cpu.
   integer:: blk_ps_max_size !< Maximum size of phi block per cpu.

   logical, private:: initialised=.false.
   logical, private:: dims_initialised=.false.
   ! Interfaces and overloads
   !> Encapsulates sums of several types and ranks.
   interface sum_over_all_cpus
      module procedure sum_over_all_cpus_scal, sum_over_all_cpus_vect
   end interface
   !> Encapsulates minimization of several types and ranks.
   interface drs_minimize
      module procedure drs_minimize_dble, drs_minimize_dble_scal
   end interface
   !> Encapsulates maximization of several types and ranks.
   interface drs_maximize
      module procedure drs_maximize_dble, drs_maximize_dble_scal
   end interface
   !> Encapsulates broadcast of several types and ranks.
   interface drs_bcast
      module procedure drs_bcast_int, drs_bcast_dble, &
                       drs_bcast_int_scal, drs_bcast_dble_scal, &
                       drs_bcast_logical_scal
   end interface
contains

   !------------------------------------------------------------------
   !> Gets initial values for mpi_size and mpi_rank. Allocates block indices accordingly.
   subroutine drs_mpi_init(error)
      implicit none
      integer, intent(out):: error
      error=0
      if(initialised) then
         error=-1
         return
      endif
#ifdef MPI
      call MPI_INIT(mpi_err)
      if(mpi_err.ne.0) then
         write(*,*) 'error in MPI_INIT().'
         error = mpi_err
         return
      endif
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpi_err)
      if(mpi_err.ne.0) then
         write(*,*) 'error in MPI_COMM_SIZE().'
         error = mpi_err
         return
      endif
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_err)
      if(mpi_err.ne.0) then
         write(*,*) 'error in MPI_COMM_RANK().'
         error = mpi_err
         return
      endif
#else
      mpi_size = 1
      mpi_rank = 0
#endif
      allocate(blk_ps_size(0:mpi_size-1),blk_ps_start(0:mpi_size-1))
      allocate(blk_t_size(0:mpi_size-1), blk_t_start(0:mpi_size-1))

      if(mpi_size.le.0) then
         spew 'Illegal number of processors: ',mpi_size
         call drs_abort(ERR_ILLEGAL_NPROCS)
      endif
      initialised = .true.
   end subroutine drs_mpi_init

   !------------------------------------------------------------------
   !> Initializes mpi variables and sizes.
   subroutine mpi_dims_init(Nt, Np_s, m0, error)
      implicit none
      integer, intent(out):: error
      integer, intent(in):: Nt, Np_s, m0
      integer:: blk_size, blk_start
      integer:: rank

      error=0
      if(dims_initialised) then
         error=-1
         return
      endif
      if(.not.initialised) then
         error = ERR_MPI_NOT_INITIALISED
         return
      endif
      !TODO: These conditions do not strictly make sense. They need to be reviewed.
      if(Np_s .lt. 2*mpi_size-1) then
         spew 'drs_dims_init(): Too many processors for this m-resolution.',Np_s
         error = ERR_TOO_MANY_PROCS_FOR_M
         return
      endif
      if(Nt+1 .lt. mpi_size) then
         spew 'drs_dims_init(): Too many processors for this theta-resolution.',Nt
         error = ERR_TOO_MANY_PROCS_FOR_TH
         return
      endif
#ifdef MPI
      spew "# CPU's in use: ",mpi_size
      if(.not.allocated(mpibuffer)) allocate(mpibuffer(4*mpi_size*Nt*Np_s*m0))
      call MPI_BUFFER_ATTACH(mpibuffer, 4*mpi_size*Nt*Np_s*m0, mpi_err)
#endif
      !-- special case mpi_size=1:
      if(mpi_size.eq.1) then
         blk_ps_start(0)  = 1
         blk_ps_size(0)   = Np_s
         ! TODO Deal with the cases where Np_s is even (last coeff does not have imaginary part).
         !if(mod(blk_ps_size(0),2).eq.0) blk_ps_size(0) = blk_ps_size(0) - 1
         blk_ps_max_size  = blk_ps_size(0)
         blk_t_start(0)   = 0
         blk_t_size(0)    = Nt + 1
         blk_t_max_size   = blk_t_size(0)
         call blk_ts_start_init(m0)
         error = 0
         return
      endif

      blk_ps_max_size = 1
      blk_start       = 1
      ! Number of azimuthal harmonics per cpu
      ! Also coincides with the block size in m
      ! The first coeficient corresponds to the static term (treated by
      ! mpi_rank=0). We have Np_s-1 terms that will be distributed by all the CPU's.
      blk_size = ceiling(dble(Np_s-1)/dble(mpi_size))
      ! Block size must be an even number to contain all real and imag coefficients.
      ! If it's not, we increment blk_size by one so that each CPU has
      ! consistent set of coeficients to treat. As a consequence, the last CPU will
      ! have a bit less work to do.
      if(mod(blk_size,2).ne.0) blk_size = blk_size + 1

      blk_ps_start(0) = blk_start
      ! Add one more coeff for m=0 to the first CPU.
      blk_ps_size(0)  = blk_size + 1
      ! Start of the next block
      blk_start       = blk_ps_start(0) + blk_ps_size(0)
      blk_ps_max_size = blk_ps_size(0)
      !-- loop over processors except last:
      do rank=1, mpi_size-2
         blk_size = ((Np_s - 1) - blk_start + 1)/(mpi_size-rank)
         if(mod(blk_size,2).ne.0) blk_size = blk_size + 1
         blk_ps_max_size = max(blk_ps_max_size, blk_size)
         blk_ps_start(rank) = blk_start
         blk_ps_size(rank)  = blk_size
         blk_start = blk_start + blk_size
      enddo
      !-- last processor:
      blk_size = Np_s - 1 - blk_start + 1
      blk_ps_start(mpi_size - 1) = blk_start
      blk_ps_size(mpi_size - 1)  = blk_size
      blk_ps_max_size = max(blk_ps_max_size, blk_size)

      ! Number of meridional points per cpu
      blk_t_max_size = ceiling(dble(Nt+1)/dble(mpi_size))
      blk_size  = blk_t_max_size
      blk_start = 0
      do rank=0, mpi_size-1
         blk_size = (Nt - blk_start + 1)/(mpi_size-rank)
         blk_t_max_size = max(blk_t_max_size, blk_size)
         blk_t_start(rank) = blk_start
         blk_t_size(rank)  = blk_size
         blk_start = blk_start + blk_size
      enddo
      ! everyone precomputes its own blk_ts_start(:)
      call blk_ts_start_init(m0)
      do rank=0, mpi_size-1
         spew "CPU", rank,"t-blocks:", blk_t_start(rank),"->", blk_t_start(rank)+blk_t_size(rank)-1
      enddo
      do rank=0, mpi_size-1
         spew "CPU", rank,"ps-blocks:", blk_ps_start(rank),"->", blk_ps_start(rank)+blk_ps_size(rank)-1
      enddo
      error = 0
      dims_initialised = .true.
   end subroutine mpi_dims_init

   !------------------------------------------------------------------
   !> transposition: t distrib(phi) --> tt_t distrib(theta):
   !!
   !!~~~~~~~~
   !!  input: t(0:Nt,Np_s) in (theta,phi)       (distr. in phi)
   !!         blk_ps_size(),blk_ps_start() contain the local input dims in phi for all processors.
   !!         blk_t_size(),blk_t_start() contain the local output dims in theta.
   !!         blk_ps_max_size,blk_t_max_size: maximum blocksizes in phi, theta.
   !!         mpi_size,mpi_rank
   !!         ir,tag: radial index,message-tag
   !!
   !! output:  tt_t(0:(Ntl-1),mpi_size*blk_ps_max_size)     (distr. in theta)
   !!          tt_t(0:(blk_t_max_size-1),Np)       (dynamic physical dims!)
   !!~~~~~~~~
   !! 04.10.96 M.A. original version (blocking).
   !-----------------------------------------------------------------------
   subroutine transpos_phi2theta(input, Nt, output, Np)
      implicit none
      double precision, intent(in):: input(0:Nt, 1:blk_ps_size(mpi_rank))
      integer, intent(in):: Np, Nt
      double precision, intent(out):: output(0:blk_t_size(mpi_rank)-1, 1:Np)    !  dynamic physical dims.
      integer:: j, l
#ifdef MPI
      integer:: tag, rank, n
      integer:: t_start, t_end, ps_start, ps_end
#endif
      output = 0.0d0
#ifdef MPI
      to_rank: do n=1, mpi_size-1
         ! Receiving CPU is the next in rank until there are no more ranks.
         ! Then we restart from 0.
         rank = mod(mpi_rank + n, mpi_size)
         ! Cut the data into pieces in theta:
         t_start = blk_t_start(rank)
         t_end   = blk_t_start(rank) + blk_t_size(rank) - 1
         ! construct a tag that encodes the rank of the sender and of the
         ! receiver, that is 1000*sender+receiver (works for less than 1000
         ! cpu's.)
         tag = 1000*mpi_rank + rank
         ! Send all mpi_rank has in phi and as much as rank needs in theta
         call mpi_bsend( input(t_start:t_end, 1:blk_ps_size(mpi_rank)), &
                         blk_t_size(rank)*blk_ps_size(mpi_rank), &
                         MPI_DOUBLE_PRECISION, &
                         rank, &
                         tag, &
                         MPI_COMM_WORLD, &
                         mpi_err )
      enddo to_rank

      from_rank: do n=1, mpi_size-1
         rank = mpi_rank - n
         if(rank.lt.0) rank = mpi_size + rank

         !-- put the pieces in phi together:
         ps_start = blk_ps_start(rank)
         ps_end   = blk_ps_start(rank) + blk_ps_size(rank) - 1

         ! Construct a tag that encodes the rank of the sender and of the
         ! receiver, that is 1000*sender+receiver (works for less than 1000
         ! cpu's.)
         tag = 1000*rank + mpi_rank
         call mpi_recv( output(0:blk_t_size(mpi_rank)-1, ps_start:ps_end), &
                        blk_t_size(mpi_rank)*blk_ps_size(rank), &
                        MPI_DOUBLE_PRECISION, &
                        rank, &
                        tag, &
                        MPI_COMM_WORLD, &
                        mpi_stat,&
                        mpi_err )
      enddo from_rank
#endif
      !-- mpi_rank --> mpi_rank:
      forall( l=0:blk_t_size(mpi_rank) - 1, j=1:blk_ps_size(mpi_rank) )
            output(l, blk_ps_start(mpi_rank) + j - 1) = input(blk_t_start(mpi_rank) + l, j)
      endforall
   end subroutine transpos_phi2theta

   !-----------------------------------------------------------------------
   !> transposition: tt_t distrib(theta) --> t distrib(phi):
   !!
   !!~~~~~~~
   !!  input:  tt_t(0:Ntl,mpi_size*blk_ps_max_size)     (transposed)
   !!          tt_t(0:(blk_t_max_size-1),Np)   (dynamic physical dims!)
   !!
   !!         blk_ps_size(),blk_ps_start() contain the local output dims in phi for all pes.
   !!         blk_t_size(),blk_t_start() contain the local input dims in theta.
   !!         blk_ps_max_size,blk_t_max_size: maximum blocksizes in phi, theta.
   !!         mpi_size,mpi_rank
   !!         ir,tag: radial index,message-tag
   !!
   !! output: t(0:Nt,Np_s) in (theta,phi)    (distr. in phi)
   !!~~~~~~~
   !!
   !! 04.10.96 M.A. original version.
   subroutine transpos_theta2phi(input, Np_s, output, Nt)
      implicit none
      integer, intent(in):: Np_s, Nt
      double precision, intent(in):: input(0:(blk_t_size(mpi_rank)-1),Np_s)
      double precision, intent(out):: output(0:Nt,blk_ps_size(mpi_rank))
      integer:: j,l
#ifdef MPI
      integer:: rank, n, tag
      integer:: t_start, t_end, ps_start, ps_end
#endif
      output = 0.0d0
#ifdef MPI
      ! Each cpu contains all values in the phi/m direction.
      ! It sends its data to the correct cpu with the corect size.
      if (mpi_size.gt.1) then
         do n=1, mpi_size-1
            ! Receiving CPU is the next in rank until there are no more ranks.
            ! Then we restart from 0.
            rank = mod(mpi_rank + n, mpi_size)
            ! Cut the data into pieces in phi:
            ps_start = blk_ps_start(rank)
            ps_end   = blk_ps_start(rank) + blk_ps_size(rank) - 1

            ! Construct a tag that encodes the rank of the sender and of the
            ! receiver, that is 1000*sender+receiver (works for less than 1000
            ! cpu's.)
            tag = 1000*mpi_rank + rank
            ! Send all mpi_rank has in theta and as much as rank needs in phi
            call mpi_bsend( input(0:blk_t_size(mpi_rank)-1, ps_start:ps_end), &
                           blk_t_size(mpi_rank)*blk_ps_size(rank),&
                           MPI_DOUBLE_PRECISION,&
                           rank,&
                           tag,&
                           MPI_COMM_WORLD,&
                           mpi_err )
         enddo  !  rank

         do n=1, mpi_size-1
            rank = mpi_rank - n
            if(rank.lt.0) rank = mpi_size + rank

            !-- put the pieces in theta together:
            t_start = blk_t_start(rank)
            t_end   = blk_t_start(rank) + blk_t_size(rank) - 1

            ! Construct a tag that encodes the rank of the sender and of the
            ! receiver, that is 1000*sender+receiver (works for less than 1000
            ! cpu's.)
            tag = 1000*rank + mpi_rank
            call mpi_recv( output(t_start:t_end, 1:blk_ps_size(mpi_rank)), &
                           blk_t_size(rank)*blk_ps_size(mpi_rank),&
                           MPI_DOUBLE_PRECISION,&
                           rank,&
                           tag,&
                           MPI_COMM_WORLD,&
                           mpi_stat,&
                           mpi_err )

         enddo  !  rank
      endif
#endif
      !-- mpi_rank --> mpi_rank:
      forall( l=0:blk_t_size(mpi_rank) - 1, j=1:blk_ps_size(mpi_rank) )
            output(blk_t_start(mpi_rank) + l, j) = input(l, blk_ps_start(mpi_rank) + j - 1)
      endforall
   end subroutine transpos_theta2phi

   !------------------------------------------------------------------
   !> Performs a one-to-all communication of the contents of buffer. It is essentially a
   !! targeted version of mpi_scatter.
   subroutine distribute_in_m(buffer, Nt, Nr)
      implicit none
      integer, intent(in):: Nt, Nr
      double precision, intent(inout):: buffer(:,:,:)
#ifdef MPI
      integer:: rank, blk_start, blk_end
      if(mpi_rank.eq.0) then
         !-- get data for other pe's:
         do rank=1, mpi_size-1
            blk_start = blk_ps_start(rank)
            blk_end   = blk_start + blk_ps_size(rank) - 1
            call MPI_SEND(buffer(1:Nr, 1:Nt+1, blk_start:blk_end), &
                       (Nt+1)*blk_ps_size(rank)*Nr, &
                       MPI_DOUBLE_PRECISION, &
                       rank, &
                       100+rank, &
                       MPI_COMM_WORLD, &
                       mpi_err)
        enddo  !  rank=1,mpi_size-1
     else   !  end root
         !-- receive data from root:
         call MPI_RECV(buffer(1:Nr, 1:Nt+1, 1:blk_ps_size(mpi_rank)), &
                     (Nt+1)*blk_ps_size(mpi_rank)*Nr, &
                     MPI_DOUBLE_PRECISION, &
                     0, &
                     100+mpi_rank, &
                     MPI_COMM_WORLD, &
                     mpi_stat, &
                     mpi_err)
     endif  ! end slave
#endif
   end subroutine distribute_in_m

   !------------------------------------------------------------------
   !> Performs an all-to-one communication of the contents of buffer. It is essentially a
   !! targeted version of mpi_gather.
   subroutine gather_from_m(buffer, Nt, Nr)
      implicit none
      integer, intent(in):: Nt, Nr
      double precision, intent(inout):: buffer(:,:,:)
#ifdef MPI
      double precision,allocatable:: aux(:,:,:)
      integer:: rank, blk_start, blk_end
      if(mpi_rank.ne.0) then !-- send data to root:
         call MPI_SEND(buffer(:,:,:), &
                        (Nt+1)*blk_ps_size(mpi_rank)*Nr, &
                        MPI_DOUBLE_PRECISION, &
                        0, &
                        200 + mpi_rank, &
                        MPI_COMM_WORLD, &
                        mpi_stat, &
                        mpi_err)
      else !-- get data from other pe's:
         do rank=1, mpi_size-1
            allocate(aux(Nt+1, blk_ps_size(rank), Nr))
            blk_start = blk_ps_start(rank)
            blk_end   = blk_start + blk_ps_size(rank) - 1
            call MPI_RECV(aux(:,:,:), &
                          (Nt+1)*blk_ps_size(rank)*Nr, &
                          MPI_DOUBLE_PRECISION, &
                          rank, &
                          200 + rank, &
                          MPI_COMM_WORLD, &
                          mpi_stat, &
                          mpi_err)
            buffer(:,blk_start:blk_end,:) = aux(:,:,:)
            deallocate(aux)
         enddo  !  rank=1,mpi_size-1
      endif
#else
      buffer = buffer
#endif
   end subroutine gather_from_m

   !------------------------------------------------------------------
   !> Initialises blk_ts_start.
   !!
   subroutine blk_ts_start_init(m0)
      implicit none
      integer, intent(in):: m0
      integer:: j
      if(.not.allocated(blk_ts_start)) allocate(blk_ts_start(blk_ps_max_size))
      do j=1, blk_ps_max_size
         blk_ts_start(j) = m0*((blk_ps_start(mpi_rank) + j - 1)/2)
      enddo
      mm => blk_ts_start
   end subroutine

   !------------------------------------------------------------------
   !> Subroutine to encapsulate sums across all the cpu's.\n
   subroutine sum_over_all_cpus_scal(val)
      implicit none
      double precision, intent(inout)::val
#ifdef MPI
      double precision:: aux

      if(mpi_size.gt.1) then
         call mpi_allreduce(val,aux,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_err)
         val = aux
      endif
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   !> Subroutine to encapsulate mpi calls that sum arrays over all cpu's
   subroutine sum_over_all_cpus_vect(val)
      implicit none
      double precision, intent(inout)::val(:)
#ifdef MPI
      integer:: length
      double precision,allocatable:: aux(:)

      length=size(val,1)
      allocate(aux(length))

      if(mpi_size.gt.1) then
         call mpi_reduce(val,aux,length,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)
         val = aux
      endif
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   !> Encapsulate mpi barrier
   subroutine wait_for_everyone()
      implicit none
#ifdef MPI
      if(mpi_size.gt.1) then
         call mpi_barrier(MPI_COMM_WORLD,mpi_err)
      endif
#endif
   end subroutine wait_for_everyone

   !------------------------------------------------------------------
   !> Encapsulate mpi_reduce min
   subroutine drs_minimize_dble(array)
      implicit none
      double precision, intent(inout):: array(:)
#ifdef MPI
      integer::length
      double precision, allocatable:: aux(:)
      length=size(array,1)
      allocate(aux(length))
      if(mpi_size.gt.1) then
         call MPI_REDUCE(array(1:length), aux(1:length), length, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, mpi_err)
         array = aux
      endif
      deallocate(aux)
#else
      array = array
#endif
   end subroutine

   !------------------------------------------------------------------
   !> Encapsulate mpi_reduce min (scalars)
   subroutine drs_minimize_dble_scal(val)
      implicit none
      double precision, intent(inout):: val
#ifdef MPI
      double precision:: aux
      if(mpi_size.gt.1) then
         call mpi_reduce(val, aux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, mpi_err)
         val = aux
      endif
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   !> Encapsulate mpi_reduce max
   subroutine drs_maximize_dble(array)
      implicit none
      double precision, intent(inout):: array(:)
#ifdef MPI
      integer::length
      double precision,allocatable:: aux(:)
      length=size(array,1)
      allocate(aux(length))
      if(mpi_size.gt.1) then
         call mpi_reduce(array, aux, length, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)
         array = aux
      endif
      deallocate(aux)
#else
      array = array
#endif
   end subroutine

   !------------------------------------------------------------------
   !> Encapsulate mpi_reduce max (scalars)
   subroutine drs_maximize_dble_scal(val)
      implicit none
      double precision, intent(inout):: val
#ifdef MPI
      double precision:: aux
      if(mpi_size.gt.1) then
         call mpi_reduce(val, aux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, mpi_err)
         val = aux
      endif
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   !> Collects double precision values stored in \a val from the ranks \a rank
   !! on the root process.
   subroutine drs_gather_vars(rank, val)
      implicit none
      integer, intent(in)::rank(:) !< list of ranks to collect from.
      double precision, intent(inout)::val(:) !< Values
#ifdef MPI
      integer:: N, i
      if(mpi_size.gt.1) then
         N = size(rank)
         do i=1, N
            if(rank(i)==0) cycle ! No need to send from root to root.
            if(mpi_rank.eq.rank(i)) then
               call mpi_bsend(val(i), 1, MPI_DOUBLE_PRECISION, 0, 200+i, MPI_COMM_WORLD, mpi_err)
            endif
            if(mpi_rank.eq.0) then
               call mpi_recv(val(i), 1, MPI_DOUBLE_PRECISION, rank(i), 200+i, MPI_COMM_WORLD, mpi_stat, mpi_err)
            endif
         enddo
      endif
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   subroutine drs_bcast_dble(array, num)
      implicit none
      double precision, intent(inout):: array(:)
      integer, intent(in):: num
#ifdef MPI
      if(mpi_size.gt.1) then
         call mpi_bcast(array, num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
      endif
#else
      array = array
#endif
   end subroutine

   !------------------------------------------------------------------
   subroutine drs_bcast_int(array, num)
      implicit none
      integer, intent(inout):: array(:)
      integer, intent(in):: num
#ifdef MPI
      !-- broadcast input parameters:
      call mpi_bcast(array, num, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#else
      array = array
#endif
   end subroutine

   !------------------------------------------------------------------
   subroutine drs_bcast_dble_scal(val)
      implicit none
      double precision, intent(inout):: val
#ifdef MPI
      !-- broadcast input parameters:
      call mpi_bcast(val, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   subroutine drs_bcast_int_scal(val)
      implicit none
      integer, intent(inout):: val
#ifdef MPI
      !-- broadcast input parameters:
      call mpi_bcast(val, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   subroutine drs_bcast_logical_scal(val)
      implicit none
      logical, intent(inout):: val
#ifdef MPI
      !-- broadcast input parameters:
      call mpi_bcast(val, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)
#else
      val = val
#endif
   end subroutine

   !------------------------------------------------------------------
   subroutine drs_abort(error)
      implicit none
      integer, intent(in):: error
      Write(*,*) "CPU", mpi_rank, "exiting with error", error
#ifdef MPI
      call MPI_ABORT(MPI_COMM_WORLD, error, mpi_err)
#else
      stop
#endif
   end subroutine

   !------------------------------------------------------------------
   subroutine mpi_cleanup()
      implicit none
#ifdef MPI
      integer:: freed_size
      call MPI_BUFFER_DETACH(mpibuffer, freed_size, mpi_err)
      if(allocated(mpibuffer)) deallocate(mpibuffer)
      call MPI_FINALIZE(mpi_err)
#endif
   end subroutine

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
