! matrix inversion routine
! a is a LU decomposed matrix, indx the index array returned by ludcmp,
! n the dimension and np the physical size of the arrays. y is a work array
! copied from numerical recipes
subroutine matinv(a, indx, n, np)
   implicit none
   integer, intent(in):: n,np,indx(n)
   double precision, intent(inout):: a(np,np)
   double precision:: y(np,np)
   integer:: i

   y = a
   a = 0.0d0
   do i=1,n
      a(i,i)=1.
   enddo
   do i=1,n
      call lubksb(y,n,np,indx,a(1:np,i))
   enddo
end subroutine matinv
