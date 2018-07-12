! Copyright  M. Ardes (), 1996
! Copyright  L. Silva (lacsilva@gmail.com), 2013-2015
module drs_hypDiff
   implicit none
   private
   double precision, allocatable::hypDiff(:) ! (0:Nt)
   logical:: drs_want_hypDiff = .FALSE.
   interface drs_apply_hypDiff
      module procedure drs_apply_hypDiff_dble, drs_apply_hypDiff_cmplx
   end interface

   public::drs_want_hypDiff, drs_hypDiff_init, drs_apply_hypDiff
contains

   subroutine drs_hypDiff_init(Nt)
      implicit none
      integer, intent(in):: Nt
      integer:: l
      if (drs_want_hypDiff) then
         allocate(hypDiff(0:Nt))
         do l=0,Nt
            ! hypdiff(l)=1.0
            ! hypdiff(l)=1+((1.0*l)/Nt_s)**3
            hypdiff(l) = 1.0d0 + 0.075d0*(l**3)
         enddo
      endif
   end subroutine drs_hypDiff_init

   elemental subroutine drs_apply_hypDiff_dble(a,l)
      implicit none
      double precision, intent(inout):: a
      integer, intent(in):: l
      if (drs_want_hypDiff) a = a*hypDiff(l)
   end subroutine drs_apply_hypDiff_dble

   elemental subroutine drs_apply_hypDiff_cmplx(a,l)
      implicit none
      double complex, intent(inout):: a
      integer, intent(in):: l
      if (drs_want_hypDiff) a = a*hypDiff(l)
   end subroutine drs_apply_hypDiff_cmplx
end module drs_hypDiff
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
