C Reinhard Furrer, fall 2019
C
C
C Painful implementation of a simple shift of elements...
C     Here, shift of elements corresponds to flipping a (0,1) pair
C     in a sparse matrix A. 

C Rudimentary tests are performed, errors are passed back by `ncol`
C   ncol=0:     success
C   ncol=-999:  value to replace was not zero
C   ncol=nrow:  replacing value was already zero 
C   ncol=-99:   (i,j)=(k,l)
      
c A more flexible version `exchange` is available upon request.
c     In that version,  a shift of elements corresponds to flipping
C     (`flip.eq.1`) a  (0,value) pair in a sparse matrix A.
c     Alternatively, the value can be passed via `val`.
      
      subroutine exchange2(nrow, ncol, a, ja, ia, i,j, k,l)
   
      implicit none
      
      integer nrow, ncol
      double precision a(*), val  
      integer  ja(*), ia(nrow+1)
      integer  i,j, k,l

      integer kk,rr, lpos, ipos
c     lpos contains the number of elements before l (case 2)
c     temporary assign a value to ipos  [-Wmaybe-uninitialized]
      ipos = 1


C     case 1: position in vector 


      val = 1.0
C For simplicity, we distinguish different cases, shifting forward
C backward, or same line manipulations. To not loose elements we have
C to follow the shift.

      if ((i .lt. 1).or.(k .lt. 1).or. (j.lt.1) .or. (k .lt. 1)) then
         ncol=-1
         return
      endif 
      if ((i .gt. nrow).or.(k .gt. nrow).or. (j .gt. ncol) .or. 
     &                   (l .gt. ncol)) then
         ncol=-2
         return
      endif 
 

C*CASE 1*   Changes shifting down/forward/right...      
      if (k.lt.i) then
         lpos=ia(k+1)
         do kk=ia(k),ia(k+1)-1 
            if (ja(kk) .eq. l) goto 999
            if (ja(kk) .gt. l) then
              lpos=kk
              goto 15
            endif
         enddo   
 15      continue
         do kk=ia(i),ia(i+1)-1
            if (ja(kk) .eq. j) then
               ipos = kk
               ncol = 0
            endif 
         enddo
         if (ncol.ne.0) return
         
         do rr=k+1,i
            ia(rr) = ia(rr) + 1
         enddo
         if (ipos.eq.lpos) then
             ja(lpos) = l  
             a(lpos) = val
             return 
             
         else 
             do kk=ipos,lpos+1,-1
               ja(kk) = ja(kk-1)
               a(kk) = a(kk-1)
             enddo
             ja(lpos) = l  
             a(lpos) = val
             return
         endif    
            
         ncol = -10
         return
       endif

C*CASE 2*    Changes shifting down      
      if (k.gt.i) then
         lpos=ia(k)-1
         do kk=ia(k+1)-1,ia(k),-1
            if (ja(kk) .eq. l) goto 999
            if (ja(kk) .lt. l) then
              lpos=kk
              goto 35
            endif
         enddo
 35      continue
         do kk=ia(i),ia(i+1)-1
            if (ja(kk) .eq. j) then
               ipos = kk
               ncol = 0
            endif 
         enddo
         if (ncol.ne.0) return
         
         do rr=i+1,k
            ia(rr) = ia(rr) - 1
         enddo
         if (ipos.eq.lpos) then
             ja(lpos) = l  
             a(lpos) = val
             return 
             
         else 
             do kk=ipos,lpos-1
               ja(kk) = ja(kk+1)
               a(kk) = a(kk+1)
             enddo
             ja(lpos) = l  
             a(lpos) = val
             return
         endif    
            
         ncol = -20
         return
       endif
       
       
   
c*CASE 3*     Changes are in same row. Three cases again
      if (l .lt. j) then
c     Move elements right, from kk-1 -> kk
c     We first test if there are any moves necessary. for simplicity
c     we work with two sub cases. 
         do kk=ia(k+1)-1,ia(k),-1
            if (ja(kk) .eq. j) then
               ncol = 0
               
               if (kk .eq. ia(k)) then
                 ja(kk) = l
                 a(kk) = val
                 return 
               else if (ja(kk-1) .lt. l) then
                 ja(kk) = l
                 a(kk) = val
                 return 
               endif                 
                 
               
               
            endif 
            if ((ja(kk) .le. j) .and. ((ja(kk) .gt. l)))then
               ja(kk) = ja(kk-1)
               a(kk) = a(kk-1)
            endif   
            if (ja(kk) .eq. l) goto 999
            if (ja(kk) .lt. l) then
               ja(kk) = l
               a(kk) = val
               return 
            endif   
         enddo
         ncol = -31
         return
      endif
      if (l .gt. j ) then
c move elements left, from kk+1 -> kk      
         do kk=ia(k),ia(k+1)-1
            if (ja(kk) .eq. j) then
               ncol = 0
               if(kk .eq. ia(k+1)-1) then
                  ja(kk) = l
                  a(kk) = val
                  return 
               else if (ja(kk+1).gt.l) then        
                  ja(kk) = l
                  a(kk) = val
                  return 
               endif
            endif
            if ((ja(kk) .ge. j) .and. (ja(kk) .lt. l)) then
               ja(kk) = ja(kk+1)
               a(kk) = a(kk+1)
            endif   
            if (ja(kk) .eq. l) goto 999
            if (ja(kk) .gt. l) then
               ja(kk) = l
               a(kk) = val
               return 
            endif   
         enddo
         ncol = -32
         return
      endif
      
C   Final case, (i,j)=(k,l), error by definition!      
      ncol = -99
      return
      
 999  ncol = -999
      return 
      end
      
      
      
