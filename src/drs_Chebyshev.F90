! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
!> Module containing the implementation of the Chebyshev polynomials.
module drs_Chebyshev
   implicit none
   save
   !> It makes use of fftw3.
   include  'fftw3.f'
   private
   double precision, parameter:: pi  = 3.141592653589793d0
   integer*8:: plan_x
   integer:: Nx, Nx_s
   double precision, allocatable:: ct_buffer(:)
   double precision, allocatable, public:: Cheb_x(:)
   double precision, allocatable, public, target:: Chebyshev(:,:)
   double precision, allocatable, public, target:: Chebyshev_dx(:,:)
   !> First index is radial point, second index is mode index
   double precision, allocatable, public, target:: Chebyshev_ddx(:,:)
   public:: Chebyshev_init, Chebyshev_cleanup
   public:: Chebyshev_x2n, Chebyshev_n2x
   public:: Cheb_compute_dx_ddx_x2x
   public:: Cheb_compute_dx_ddx_n2x, Cheb_compute_dx_n2n
contains
   !-----------------------------------------------------------------------------
   !> Computes the Chebyshev polynomials of order up to \a N as a function of r.
   !! \since 1.6.5
   subroutine Chebyshev_init(N, N_s)
      implicit none
      !> Number of points the polynomials
      integer, intent(in):: N
      !> maximum order of the polynomials
      integer, intent(in):: N_s
      integer::i
      Nx   = N
      Nx_s = N_s
      allocate(Chebyshev(Nx, Nx), Chebyshev_dx(Nx, Nx), Chebyshev_ddx(Nx, Nx))
      allocate(Cheb_x(Nx), ct_buffer(Nx))
      Chebyshev = 0.0d0
      call dfftw_plan_r2r_1d(plan_x, Nx, ct_buffer, ct_buffer, FFTW_REDFT00, FFTW_MEASURE)
      do i=1, Nx
         Cheb_x(i) = dcos(pi*(i-1.d0)/(dble(Nx)-1.d0))
         Chebyshev(i,i) = 1.0d0
         call Cheb_compute_dx_ddx_n2x(Chebyshev(:,i),Chebyshev_dx(:,i),Chebyshev_ddx(:,i))
         call Chebyshev_n2x(Chebyshev(:,i))
      enddo
   end subroutine Chebyshev_init

   !-----------------------------------------------------------------------------
   !! \since 1.6.5
   subroutine Chebyshev_cleanup()
      implicit none
      call dfftw_destroy_plan(plan_x)
      deallocate(Cheb_x, ct_buffer, Chebyshev, Chebyshev_dx, Chebyshev_ddx)
   end subroutine Chebyshev_cleanup

   !-----------------------------------------------------------------------------
   !> Returns second radial derivative in d2fdx2, first derivative in dfdx
   !! Input f is supposed to be given in Chebychev space, derivatives are returned
   !! in direct space.
   !! \since 1.6.5
   subroutine Cheb_compute_dx_ddx_n2x(f, dfdx, d2fdx2)
      implicit none
      !> Chebyshev coefficients of the input function
      double precision, intent(in):: f(Nx)
      !> First derivative of f at points 1..Nx
      double precision, intent(out):: dfdx(Nx)
      !> Second derivative of f at points 1..Nx
      double precision, intent(out):: d2fdx2(Nx)
      ! Differentiate in spectral space
      call Cheb_compute_dx_n2n(f, dfdx)
      ! Again for the second radial derivative
      call Cheb_compute_dx_n2n(dfdx, d2fdx2)
      ! Dealiase
      if (Nx_s+1 .le. Nx) dfdx(Nx_s+1:Nx) = 0.d0
      if (Nx_s+1 .le. Nx) d2fdx2(Nx_s+1:Nx) = 0.d0
      ! Convert to real space
      call Chebyshev_n2x(dfdx)
      call Chebyshev_n2x(d2fdx2)
   end subroutine Cheb_compute_dx_ddx_n2x

   !-----------------------------------------------------------------------------
   !> Returns second radial derivative in d2fdx2, first derivative in dfdx
   !! Input f is supposed to be given in real space, derivatives are returned
   !! in real space.
   !! \since 1.6.5
   subroutine Cheb_compute_dx_ddx_x2x(f, dfdx, d2fdx2)
      implicit none
      !> Chebyshev coefficients of the input function
      double precision, intent(in):: f(Nx)
      !> First derivative of f at points 1..Nx
      double precision, intent(out):: dfdx(Nx)
      !> Second derivative of f at points 1..Nx
      double precision, intent(out):: d2fdx2(Nx)
      double precision:: aux(Nx)
      aux=f
      call Chebyshev_x2n(aux(1:Nx))
      ! Differentiate in spectral space
      call Cheb_compute_dx_n2n(aux, dfdx)
      ! Again for the second radial derivative
      call Cheb_compute_dx_n2n(dfdx, d2fdx2)
      ! Dealiase
      if (Nx_s+1 .le. Nx) dfdx(Nx_s+1:Nx) = 0.d0
      if (Nx_s+1 .le. Nx) d2fdx2(Nx_s+1:Nx) = 0.d0
      ! Convert to real space
      call Chebyshev_n2x(dfdx)
      call Chebyshev_n2x(d2fdx2)
   end subroutine Cheb_compute_dx_ddx_x2x

   !-----------------------------------------------------------------------------
   !> Computes the Chebyshev coefficients of the first derivative of \a f with
   !! respect to x.
   !! \since 1.6.5
   pure subroutine Cheb_compute_dx_n2n(f, dfdx)
      implicit none
      !> Chebyshev coefficients of the input function
      double precision, intent(in):: f(Nx)
      !> Chebyshev coefficients of the derivative.
      double precision, intent(out):: dfdx(Nx)
      integer:: i

      dfdx(Nx)   = 0.0d0
      dfdx(Nx-1) = 2*(Nx-1)*f(Nx)
      dfdx(Nx-2) = 2*(Nx-2)*f(Nx-1)
      do i=Nx-3,2,-1
         dfdx(i) = 2*i*f(i+1) + dfdx(i+2)
      enddo
      dfdx(1) = f(2) + dfdx(3)/2
   end subroutine Cheb_compute_dx_n2n

   !-----------------------------------------------------------------------------
   !> The forward real to spectral cosinus transform
   !! Since \f$ T_n(\cos(t)) = \cos(nt) \f$, the forward cosinus transform gives
   !! us the coefficients of order \f$n\f$ of the expansion of a scalar function
   !! \f$f(x)\f$ in terms of Chebyshev polynomials.
   !! \since 1.6.5
   subroutine Chebyshev_x2n(input)
      implicit none
      !>
      double precision, intent(inout):: input(Nx)
      call dfftw_execute_r2r(plan_x, input, input)
      input(:) = input(:)/(2*(Nx-1))
      input(2:Nx-1) = input(2:Nx-1)*2
      if (Nx_s+1 .le. Nx) input(Nx_s+1:Nx) = 0.d0
   end subroutine Chebyshev_x2n

   !-----------------------------------------------------------------------------
   !> The backward spectral to real cosinus transform
   !! Since \f$ T_n(\cos(t)) = \cos(nt) \f$, the backward cosinus transform gives
   !! us the value of a scalar function \f$f(x)\f$ in terms of Chebyshev polynomials.
   !! \since 1.6.5
   subroutine Chebyshev_n2x(input)
      implicit none
      double precision, intent(inout):: input(Nx)
      ! Dealiase
      if (Nx_s+1 .le. Nx) input(Nx_s+1:Nx) = 0.d0
      input(2:Nx-1) = input(2:Nx-1)/2
      call dfftw_execute_r2r(plan_x, input, input)
   end subroutine Chebyshev_n2x

end module
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
