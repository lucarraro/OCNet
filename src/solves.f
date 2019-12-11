c     Reinhard Furrer, 2019
      
c     The file provides a forward solve Lx=b, essentially one-to-one from
c     spam. b is vector only and thus gCC10 ready.



      subroutine spamforward (n,x,b,l,jl,il,strow)
      implicit none
      logical,  external :: eqZERO

      integer n, jl(*),il(n+1), strow
      double precision  x(n), b(n), l(*)

      integer j,k 
      double precision  t

c-----------------------------------------------------------------------
c   solves    L x = b ; L = lower triang. /  CSR format
c                        sequential forward elimination 
c-----------------------------------------------------------------------
c
c On entry:
c----------
c n      = integer. dimensions of problem.
c b      = real array containg the right side.
c
c l, jl, il,    = Lower triangular matrix stored in CSR format.
c
c    only rows strow:n are updated!
c
c On return:
c-----------
c       x  = The solution of  L x  = b.
c--------------------------------------------------------------------
c     Reinhard Furrer June 2008, April 2012, Sept 2016
      

      if (strow .eq. 1) then
c        if first diagonal element is zero, break
         if (eqZERO(l(1))) then
            k = 1
            goto 5
         endif
c     first row has one element, increase starting row
         x(1) = b(1) / l(1)
         strow = 2
      endif  

      do 3 k = strow, n
         t = b(k)
         do 1 j = il(k), il(k+1)-1
            if (jl(j) .lt. k) then
               t = t-l(j)*x(jl(j))
            else
               if (jl(j) .eq. k) then
                 if (eqZERO(l(j))) goto 5 
c     diagonal element is not zero, hence we divide and leave the loop
                   x(k) = t / l(j)
                   goto 3
                endif 
             endif
 1        continue
 3     continue      

      return
 5    n = -k
      return
      end
c-----------------------------------------------------------------------


      
      function eqZERO( a)
      implicit none

      logical eqZERO
      double precision a

      eqZERO = (abs( a) .LE. 0d0)
      return
      end function

