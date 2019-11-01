C     Reinhard Furrer, fall 2019
c
c     Routine identical to spam code....
c-----------------------------------------------------------------------
      subroutine sortrows(nrow,a,ja,ia)

      implicit none
      integer nrow
      integer ia(nrow+1),ja(*)
      double precision  a(*)
c     
c     sorts the rows according to column entries
c
c On entry:
c----------
c     nrow    -- row dimension of the matrix
c     a,ja,ia -- input matrix in sparse format
c
c On return:
c-----------
c     a,ja,ia -- cleaned matrix
c
c Notes: 
c-------
c     Reinhard Furrer 2006-09-13
c-----------------------------------------------------------------------
c     Local variables
      integer i,j,k,ko,ipos
      double precision  tmp
c
c

c     .. order the entries according to column indices
c     burble-sort is used
c
      do 190 i = 1, nrow
         do 160 ipos = ia(i), ia(i+1)-1
            do 150 j = ia(i+1)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  tmp = a(k)
                  a(k) = a(j)
                  a(j) = tmp
               endif
 150        continue
 160     continue
 190  continue
      return
c---- end of sortrows --------------------------------------------------
c-----------------------------------------------------------------------
      end
