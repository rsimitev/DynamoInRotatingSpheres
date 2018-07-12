subroutine lubksb(a, n, np, indx, b)
   implicit none
   integer, intent(in):: n, np
   integer, intent(in):: indx(n)
   double precision, intent(in):: a(np, np)
   double precision, intent(inout):: b(np)
   integer:: i, ii, ll
   double precision:: tot

   ii = 0
   do i=1, n
      ll  = indx(i)
      tot = b(ll)
      b(ll) = b(i)
      if (ii.ne.0)then
         tot = tot - sum(a(i, ii:i-1)*b(ii:i-1))
      else if (tot.ne.0.) then
         ii=i
      endif
      b(i) = tot
   enddo
   do i=n, 1, -1
      tot = b(i)
      if(i.lt.n)then
         tot = tot - sum(a(i, i+1:n)*b(i+1:n))
      endif
      b(i) = tot/a(i, i)
   enddo
end subroutine lubksb
