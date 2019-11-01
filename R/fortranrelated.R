invPerm <- function(p)
    .Call("inv_permutation", p,PACKAGE = "OCNet")
          


permuteAddSolve <- function (A, P, startrow, prevsol, pow) 
{
    
    nrow <- A@dimension[1]
    ncol <- A@dimension[2]
    nz <- A@rowpointers[nrow + 1] - 1
    
    P <- as.integer(P)
    P <- invPerm(P)
    z <- .Fortran("permandsolve",
                  nrow,
                  A@entries,A@colindices,A@rowpointers,
                  entries = vector("double",nz+nrow),  # augment for additional diagonal
                  colindices = vector("integer", nz+nrow),
                  rowpointers = vector("integer", nrow + 1),
                  P,     # permutation vector, same for rows and columns!
                  x=prevsol, # previous solution that will be updated 
                  b=rep(1.0, nrow),
                  strow=as.integer(startrow),
                  pow=pow,
                  energy=0.0,
                  NAOK = TRUE, PACKAGE = "OCNet")
    return(list(Anew=z$x, energy=z$energy))
}




allinoneF <- function(N_outlet, perm, Wt, DownNode, node,  down_new, Anode, ExpEnergy) {

    nrow <- Wt@dimension[1]
    nz <- Wt@rowpointers[nrow + 1] - 1
    pas <- .Fortran("allinone",
                    nrow, ncol=Wt@dimension[2], as.integer(N_outlet),                    
                    as.integer(DownNode), as.integer(node), as.integer(down_new), as.integer(Anode),
                    
                    entries=Wt@entries, 
                    colindices=Wt@colindices,
                    rowpointers=Wt@rowpointers,
                    oa = vector("double",nz+nrow),  # augment for additional diagonal
                    oja = vector("integer", nz+nrow),
                    oia = vector("integer", nrow + 1),
                    
                    perm=as.integer(perm),
                    upperm=integer(nrow),
                   
                    Anew=numeric(nrow),
                    b=rep(1.0, nrow),
                    strow=1L,
                    pow=as.double(ExpEnergy),
                    energy=0.0,
                    flag=0L )

   if (pas$ncol!=0) {
      if (pas$ncol==-2)   stop('Indicies cannot be larger than dimension of matrix')
      if (pas$ncol==-1)   stop('Indicies have to larger than zero')
      if (pas$ncol==-999) stop('Value to set was not zero')
      if (pas$ncol==-99)  stop('(i,j)=(k,l), no flipping necessary')
      if (pas$ncol==nrow) stop('Value was already zero. ')
      if (pas$ncol!=0)    warning(paste('Error code ',pas$ncol))
   }
   if (pas$flag==1) {
       Wt_new <- new("spam")
       slot(Wt_new, "entries", check = FALSE) <- pas$entries
       slot(Wt_new, "colindices", check = FALSE) <- pas$colindices
       slot(Wt_new, "rowpointers", check = FALSE) <- pas$rowpointers
       slot(Wt_new, "dimension", check = FALSE) <- Wt@dimension

       return(list(flag=pas$flag, Anew=pas$Anew, energy=pas$energy, perm=pas$upperm, Wt_new=Wt_new))
    } else {
        return(list(flag=pas$flag))
    }
}





