c     Reinhard Furrer, fall 2019
c     Routines specific for OCNet. 
c
c     update function in Fortran. R skeleton at the end of the file.


      subroutine iinvperm(n, old, new)
      implicit none
      integer n, old(n), new(n)

      integer i
      do i=1,n
         new( old(i)) = i
      enddo
      return
      end
      
c     Outlet <- which(DownNode==0)
c        here we know how many zeros there are!
      subroutine ifindloc(n,array, value, no, pos)  
      implicit none
      integer n, array(n), value, no, pos(no)
      integer ii, jj
c     array==value, no positions are written in pos 
            
      jj = 1
      do ii=1,n
         if (array(ii) .eq. value) then
            pos(jj) = ii
            jj = jj+1
         endif
      enddo
c      if (jj-1 .ne. no) => error
      return
      end

c     Fortran implementation of `!(k %in% outlet)` 
      PURE logical function kNotInOutlet(k, no, outlet)
      implicit none
      integer,  intent(in) :: k, no, outlet(no)
      integer ii

      kNotInOutlet = .TRUE.
      do ii=1,no
         if (k .eq. outlet(ii)) then
            kNotInOutlet = .FALSE.
            return
         endif
      enddo
      return
      end
        

c     Main routine:
c     input
c      n
c     perm, inv_perm,
c     DownNode, down_new, Anode,
c     no, noDAG, upperm        
      subroutine updatePerm( n, perm, inv_perm, DownNode,
     &     node, down_new, Anode,  no, noDAG, upperm)

      implicit none
      integer n, no, Anode, noDAG,
     &     perm(n), inv_perm(n),DownNode(n),
     &     node, down_new, upperm(n)

      integer i,k, flag_move
      integer Outlet(no)
      logical kNotInOutlet
c     flag_move <- "right"  ==1

      DownNode(node) = down_new
      

      call iinvperm(n, perm, inv_perm)
c invert permutation: j = invPerm[k] <==> k = perm[j]. It represents the order of a node      
      flag_move = 1
c if 'right'(==1), node and its upstream branch need to move rightwards in the permutation vector
      call ifindloc( n, DownNode, 0, no, Outlet)
      
      if (inv_perm(node) .gt. inv_perm(down_new)) then
c if inv_perm[node] > inv_perm[down_new], node needs to be moved leftwards
         flag_move = -1
         k = down_new

C the following might be faster          
c         do while (((inv_perm(node) .gt. inv_perm(k)))  .and.
c     &       kNotInOutlet(k, no, Outlet))
         do while ( kNotInOutlet(k, no, Outlet)  .and.
     &       (inv_perm(node) .gt. inv_perm(k)))
            k = DownNode(k)
         enddo
         if ((inv_perm(node) .eq. inv_perm(k))) then
c change to address  Impure function might not be evaluated        
            if (kNotInOutlet(k, no, Outlet)) then
              noDAG = 1
              return
            endif
         endif
      endif

c      write(*,*) "outlet, move", Outlet,flag_move
      
      if (flag_move .eq. -1) then
         k = 1
         if (inv_perm(down_new) .gt.1) then
            do i=1,inv_perm(down_new)-1
               upperm(k) = perm(i) 
               k=k+1
            enddo
         endif
         do i=(inv_perm(node)-Anode+1),inv_perm(node)
            upperm(k) = perm(i)
            k=k+1
         enddo
         do i=inv_perm(down_new),(inv_perm(node)-Anode)
            upperm(k)=perm(i)
            k=k+1
         enddo        
         if (inv_perm(node) .lt. n) then
            do i=inv_perm(node)+1,n
               upperm(k) = perm(i)
               k=k+1
            enddo
         endif
         
      endif
      if (flag_move .eq. 1) then
         k = 1
         if (inv_perm(node)-Anode+1 .gt.1) then
            do i=1,inv_perm(node)-Anode
               upperm(k) = perm(i)
               k=k+1
            enddo
         endif
         if(inv_perm(node)+1 .lt. inv_perm(down_new) ) then
            do i=inv_perm(node)+1, inv_perm(down_new)-1
               upperm(k)=perm(i)
               k=k+1
            enddo
         endif
         do i=(inv_perm(node)-Anode+1), inv_perm(node)
            upperm(k) = perm(i)
            k=k+1
         enddo
         do i=inv_perm(down_new),n
            upperm(k) = perm(i)
            k=k+1
         enddo        
      endif

      if (k .ne. n+1) then
         noDAG = -1
c         write(*,*) "Error", k
      endif
      return
      end

      
      
c$$$    if (flag_move=="left"){ # if the node needs to be moved leftwards, split vector perm accordingly
c$$$      if (inv_perm[down_new]>1){
c$$$        p1 <- perm[1:(inv_perm[down_new]-1)]
c$$$      } else {p1 <- NULL}
c$$$      p2 <- perm[inv_perm[down_new]:(inv_perm[node]-A[node])]
c$$$      p3 <- perm[(inv_perm[node]-A[node]+1):inv_perm[node]]
c$$$      if (inv_perm[node] < n){
c$$$        p4 <- perm[(inv_perm[node]+1):n]
c$$$      } else {p4 <- NULL}
c$$$      
c$$$    } else if (flag_move=="right"){ # if the node needs to be moved rightwards, split vector perm accordingly
c$$$      if ((inv_perm[node]-A[node]+1) > 1){
c$$$        p1 <- perm[1:(inv_perm[node]-A[node])]
c$$$      } else { p1 <- NULL}
c$$$      p2 <- perm[(inv_perm[node]-A[node]+1):inv_perm[node]]
c$$$      if ((inv_perm[node]+1)<inv_perm[down_new]){
c$$$        p3 <- perm[(inv_perm[node]+1):(inv_perm[down_new]-1)]
c$$$      } else {p3 <- NULL}
c$$$      p4 <- perm[inv_perm[down_new]:n]
c$$$    }
c$$$    
c$$$    perm <- c(p1,p3,p2,p4) # recompose perm by swapping the places of segments p2 and p3
c$$$  }
c$$$  
