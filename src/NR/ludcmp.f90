!-- converted to double precision
subroutine ludcmp(a,n,np,indx)
   implicit none
   integer, intent(in):: n, np
   integer, intent(out):: indx(n)
   double precision, intent(inout):: a(np, np)
   double precision, parameter::tiny=1.0D-20
   double precision:: d, aamax, sum, dum
   double precision:: vv(n)
   integer:: i,j,k, imax
   d=1.d0
   do i=1,n
      aamax=0.d0
      do j=1,n
         if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
      enddo
      !if (aamax.eq.0.D0) stop 'Error: Singular matrix.'
      vv(i)=1.D0/aamax
   enddo
   do j=1,n
      if (j.gt.1) then
         do i=1,j-1
            sum = a(i,j)
            if (i.gt.1)then
               do k=1,i-1
                  sum = sum - a(i,k)*a(k,j)
               enddo
               a(i,j)=sum
            endif
         enddo
      endif
      aamax=0.D0
      do i=j,n
         sum=a(i,j)
         if (j.gt.1)then
            do k=1,j-1
               sum = sum - a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
         endif
         dum = vv(i)*abs(sum)
         if (dum.ge.aamax) then
            imax  = i
            aamax = dum
         endif
      enddo
      if (j.ne.imax)then
         do k=1,n
            dum       = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k)    = dum
         enddo
         d = -d
         vv(imax) = vv(j)
      endif
      indx(j)=imax
      if(j.ne.n)then
         if(a(j,j).eq.0.D0) a(j,j) = tiny
         dum=1.d0/a(j,j)
         do i=j+1,n
            a(i,j)=a(i,j)*dum
         enddo
      endif
   enddo
   if(a(n,n).eq.0.D0) a(n,n) = tiny
END subroutine
