C     Reinhard Furrer, fall 2019
C
C     Functions to permute vectors and matrices.
C     Structure of code is heavily inspired by the code provided by the
C     R package `spam`

C     We have the following functions:
C     rowperm, colperm
C     permandsolve, allinone


      
C     By convention we use the construct a, ja, ia for sparse
C     input matrices and oa, oja, oia for output.
      

      subroutine rowperm (nrow,a,ja,ia,oa,oja,oia,perm)
      implicit none
      integer nrow, ja(*), ia(nrow+1)
      integer oja(*), oia(nrow+1), perm(nrow)
      double precision a(*), oa(*)

      integer i, j, k, kk, ii

c     determine pointers for output matix. 
c     
      do 10 j = 1, nrow
         i = perm(j)
         oia(i+1) = ia(j+1) - ia(j)
 10   continue
c
c get pointers from lengths
c
      oia(1) = 1
      do 20 j = 1, nrow
         oia(j+1) = oia(j+1) + oia(j)
 20   continue

      do 30 ii = 1, nrow
         kk = oia(perm(ii)) 
         do 32 k = ia(ii), ia(ii+1)-1 
            oja(kk) = ja(k) 
            oa(kk) = a(k)
            kk = kk+1
 32      continue
 30   continue

      return
      end


      subroutine colperm (nrow, oa, oja, oia, perm) 
      implicit none
      integer nrow, oja(*), oia(nrow+1), perm(*)
      double precision oa(*) 

      integer k, nnz

      nnz = oia(nrow+1)-1

      do 10 k=1,nnz
         oja(k) = perm(oja(k)) 
 10   continue

      call sortrows(nrow, oa, oja, oia)

      return

      end

      
      
C Note that oa, oja are augmented to hold an additional diagoal!
      subroutine permandsolve(nrow, a, ja, ia, oa, oja, oia, pperm, 
     &         x, b, strow, pow, energy)
      
      implicit none
      integer nrow, ja(*), ia(nrow+1), oja(*), oia(nrow+1) 
      integer pperm(nrow), strow
      double precision a(*), oa(*) 
      double precision x(nrow), b(nrow), pow, energy

      integer k

 
      call rowperm (nrow, a, ja, ia,   oa, oja, oia, pperm)
      call colperm (nrow, oa, oja, oia, pperm) 


      call idsuba (nrow, oa, oja, oia) 
      
      call spamforward (nrow, x, b, oa, oja, oia, strow)
      
      energy = 0.0
      if (abs(pow-0.5) .lt. 1.111D-5) then
         do 101 k=1,nrow
            energy = energy + sqrt(x(k))
 101      continue
      else  if (abs(pow-1) .lt. 1.111D-5) then
         do 102 k=1,nrow
            energy = energy + x(k)
 102      continue
       else  
         do 103 k=1,nrow
            energy = energy + x(k)**pow
 103      continue
      endif
     
      return
      end
c end permandsolve


      
C Note that oa, oja are augmented to hold an additional diagoal!
      subroutine allinone(nrow, ncol, no, 
     &    DownNode, node,  down_new, Anode,
     &     a, ja, ia, oa, oja, oia,
     &     pperm, upperm,
     &     x, b, strow, pow, energy,
     &     noDAG)
      
      implicit none
      integer nrow, ncol, ja(*), ia(nrow+1), oja(*), oia(nrow+1) 
      integer pperm(nrow), upperm(nrow), strow
      double precision a(*), oa(*) 
      double precision x(nrow), b(nrow), pow, energy
      integer noDAG, no, DownNode(nrow), node, down_new, Anode
      integer inv_perm(nrow)
      
      integer k

      call exchange2(nrow, ncol, a, ja, ia, DownNode(node), node,
     &                         down_new, node)

      if (ncol .ne. 0 ) then
c     error in when exchanging the entries
         return
      endif


      call updatePerm( nrow, pperm, inv_perm, DownNode, node,
     &                         down_new, Anode, no, noDAG, upperm)
 
      noDAG = noDAG + 1
      if (noDAG .eq. 1) then
         call iinvperm(nrow, upperm, inv_perm)
         
         call rowperm (nrow, a, ja, ia,   oa, oja, oia, inv_perm)
         call colperm (nrow, oa, oja, oia, inv_perm) 
c         write(*,*)inv_perm(1),inv_perm(2),inv_perm(3),inv_perm(4),
c     &    inv_perm(5),inv_perm(6)
      
         call idsuba (nrow, oa, oja, oia) 
         
         call spamforward (nrow, x, b, oa, oja, oia, strow)
      
         energy = 0.0
         if (abs(pow-0.5) .lt. 1.111D-5) then
            do 101 k=1,nrow
               energy = energy + sqrt(x(k))
 101        continue
         else  if (abs(pow-1) .lt. 1.111D-5) then
            do 102 k=1,nrow
               energy = energy + x(k)
 102        continue
         else  
            do 103 k=1,nrow
               energy = energy + x(k)**pow
 103        continue
         endif
      endif
      return
      end
c end allinone


