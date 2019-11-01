C Reinhard Furrer, fall 2019
C
C
C The following function is based on diagaddmat from spamown.f, which
c adds a diagonal matrix to a sparse matrix:  A = Diag + A 

C Here we return  I-A, thus no diagonal values, and negative indices.
C A is a square matrix, diagonal elements are not passed!
C Algorithm is in place !!
      subroutine idsuba (n, a, ja, ia) 
      implicit none
      integer n
      double precision a(*)  
      integer ja(*), ia(n+1) 

c-----------------------------------------------------------------------
c on entry:
c ---------
c n = dimension of A 
c a, ja, ia   = Matrix A in compressed sparse row format. Sorted.
c iw   = n vector of zeros.
c
c on return:
c----------
c updated matrix A
c
c Reinhard Furrer
c-----------------------------------------------------------------
      logical insert
      integer i,j, k, k1, k2, icount
      integer iw(n)
      
      
      icount = n     

      do i=1, n 
         iw(i) = 0
      enddo

c     
      do 5 i=n, 1, -1 
         k1 = ia(i)
         k2 = ia(i+1)-1 

         ia(i+1) = ia(i+1)+icount

         if ((i .gt. n) .or. (iw(i) .gt. 0)) then
c     iw(ii) equal to 0, means no diagonal element in a, we need to insert it
c     test is thus true.

c     no fill-in, only copying
            do 4 k = k2,k1,-1 
               ja(k+icount) = ja(k)
               a(k+icount) = -a(k)
 4          continue  
            iw(i) = -i
         else
            insert=.TRUE.
            if (k2.lt.k1) then
               ja(k2+icount) = i
               a(k2+icount) = 1.0
               iw(i) = k2+icount
               icount = icount-1
               insert = .FALSE.
               if (icount .eq. 0) return
            else
               do 6 k = k2,k1,-1
                  if (ja(k).gt. i) then
                     ja(k+icount) = ja(k)
                     a(k+icount) = -a(k)
                  else  if  (insert) then
                     ja(k+icount) = i
                     a(k+icount) = 1.0
                     iw(i) = k+icount
                     icount = icount-1
                     insert = .FALSE.
                     if (icount .eq. 0) return
                  endif
                  if (ja(k).lt. i) then
                     ja(k+icount) = ja(k)
                     a(k+icount) = -a(k)
                  endif
 6             continue
c     in case there is only one element, larger than i, we still need to 
c     add the diagonal element
               if  (insert) then
                   ja(k+icount) = i
                   a(k+icount) = 1.0
                   iw(i) = k+icount
                   icount = icount-1
                   insert = .FALSE.
                   if (icount .eq. 0) return
                endif
            endif
         endif 
 5    continue
      return
c-----------------------------------------------------------------------
c------------end-of-setdiagmat------------------------------------------
      end
      
      
