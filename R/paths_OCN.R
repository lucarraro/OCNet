paths_OCN <- function(OCN, level = c("RN","AG"), whichNodes = NULL,
                      includePaths = FALSE, 
                      includeDownstreamNode = FALSE, 
                      includeUnconnectedPaths = FALSE, 
                      displayUpdates = FALSE){
  
  if (length(OCN$RN$nNodes)==0){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to paths_OCN.')
  }
  
  if (is.null(whichNodes)){
    whichNodes <- vector("list")
    for (str in level){ whichNodes[[str]] <- 1:OCN[[str]]$nNodes }
  } else {level <- names(whichNodes)}
  
  for (str in level){
    if (displayUpdates){message(sprintf("%s downstream paths... \n",str), appendLF = FALSE)}
    
    wN <- whichNodes[[str]]
    if (!(OCN[[str]]$outlet %in% wN)) wN <- c(wN, OCN[[str]]$outlet)
    
    ll <- paths_cpp(OCN, whichNodes=wN, str=str, includePaths = includePaths,
                    includeDownstreamNode = includeDownstreamNode,
                    includeUnconnectedPaths = includeUnconnectedPaths)
    tmp <- list(i=ll$i, j=ll$j, values=ll$values)
    OCN[[str]][["downstreamPathLength"]] <- spam(tmp, OCN[[str]]$nNodes, OCN[[str]]$nNodes)
    
    if (includePaths){OCN[[str]][["downstreamPath"]] <- ll$downstreamPath}
    if (includeUnconnectedPaths){
      OCN[[str]][["downstreamLengthUnconnected"]] <- ll$downstreamLengthUnc}
  }
  invisible(OCN)
}

# paths_OCN <- function(OCN,
#                       level=c("RN","AG"),
#                       includePaths=FALSE,
#                       includeDownstreamNode=FALSE,
#                       includeUnconnectedPaths=FALSE,
#                       displayUpdates=FALSE){
#   
#   if (length(OCN$RN$nNodes)==0){
#     stop('Missing fields in OCN. You should run aggregate_OCN prior to paths_OCN.')
#   }
#   
#   iDP <- includeUnconnectedPaths
#   
#   if("RN" %in% level){
#     # RN
#     RN_DownstreamPathLength <- spam(0,OCN$RN$nNodes,OCN$RN$nNodes)
#     if (iDP){RN_DwnstrLength_unconnected <- matrix(0,OCN$RN$nNodes,OCN$RN$nNodes)}
#     
#     indices_down <-  matrix(0,OCN$RN$nNodes*max(1000,ceiling(OCN$RN$nNodes*0.1)),2)
#     values_down <-   numeric(OCN$RN$nNodes*max(1000,ceiling(OCN$RN$nNodes*0.1)))
#     counter_down <-  1
# 
#     if(includePaths){
#       RN_DownstreamPath <- vector("list",OCN$RN$nNodes)
#       for (i in 1:OCN$RN$nNodes){RN_DownstreamPath[[i]] <- vector("list",OCN$RN$nNodes)}
#     }
#     
#     for (i in 1:OCN$RN$nNodes){
#       
#       if(includePaths){RN_DownstreamPath[[i]][[i]] <- i}
#       if(includeDownstreamNode){RN_DownstreamPathLength[i, i] <- OCN$RN$leng[i]}
#       Path <- i
#       j <- i
#       while (j != OCN$RN$outlet){
#         j <- OCN$RN$downNode[j]
#         Path <- c(Path,j)
#         
#         if(includePaths){RN_DownstreamPath[[i]][[j]] <- Path}
#         
#         indices_down[counter_down, ] <- c(i, j)
#         if (includeDownstreamNode){
#           values_down[counter_down] <- sum(OCN$RN$leng[Path])
#         } else {
#           values_down[counter_down] <- sum(OCN$RN$leng[Path]) - OCN$RN$leng[j]  
#         }
#         counter_down <- counter_down + 1
#         
#         # Ups: contains all nodes upstream of j, barring those upstream of the penultimate path node and  j 
#         # -> all the nodes for which j is the intersection node
#         if (iDP){
#           Ups <- setdiff(OCN$RN$upstream[[j]], c(OCN$RN$upstream[[Path[length(Path) - 1]]],j) ) 
#           
#           for (u in Ups){
#             if (includeDownstreamNode){
#               RN_DwnstrLength_unconnected[i,u] <- sum(OCN$RN$leng[Path])  # count the intersection node
#             } else {
#               RN_DwnstrLength_unconnected[i,u] <- sum(OCN$RN$leng[Path]) - OCN$RN$leng[j] # don't count the intersection node
#             }
#           }
#         }
#       }
#       if (displayUpdates){
#         if ((i %% round(OCN$RN$nNodes*0.001))==0){
#           message(sprintf("RN downstream paths... %.1f%%\r",i/(1.001*OCN$RN$nNodes)*100), appendLF = FALSE)}}
#     }
#     indices_down <- indices_down[1:(counter_down-1), ]
#     values_down <- values_down[1:(counter_down-1)]
#     RN_DownstreamPathLength[indices_down] <- values_down
#     
# 
#     if (displayUpdates){
#       message("RN downstream paths... 100.0%\n", appendLF = FALSE)}
#   }
#   
#   if ("AG" %in% level){
#     # AG
#     AG_DownstreamPathLength <- spam(0,OCN$AG$nNodes,OCN$AG$nNodes)
#     if(iDP){AG_DwnstrLength_unconnected <- matrix(0,OCN$AG$nNodes,OCN$AG$nNodes)} 
#     #if(iDP){AG_DwnstrLength_unconnected <- spam(0,OCN$AG$nNodes,OCN$AG$nNodes)}
#     
#     if(includePaths){
#       AG_DownstreamPath <- vector("list",OCN$AG$nNodes)
#       for (i in 1:OCN$AG$nNodes){AG_DownstreamPath[[i]] <- vector("list",OCN$AG$nNodes)}
#     }
#     
#     indices_down <-  matrix(0,OCN$AG$nNodes*max(1000,ceiling(OCN$AG$nNodes*0.1)),2)
#     values_down <-  numeric(OCN$AG$nNodes*max(1000,ceiling(OCN$AG$nNodes*0.1)))
#     counter_down <-  1
# 
#     for (i in 1:OCN$AG$nNodes){
#       
#       if(includePaths){AG_DownstreamPath[[i]][[i]] <- i}
#       if(includeDownstreamNode){AG_DownstreamPathLength[i, i] <- OCN$AG$leng[i]}
#       
#       Path <- i
#       j <- i
#       while (j != OCN$AG$outlet){
#         j <- OCN$AG$downNode[j]
#         Path <- c(Path,j)
#         
#         if(includePaths){AG_DownstreamPath[[i]][[j]] <- Path}
#         
#         indices_down[counter_down, ] <- c(i, j)
#         if (includeDownstreamNode){
#           values_down[counter_down] <- sum(OCN$AG$leng[Path])
#         } else {
#           values_down[counter_down] <- sum(OCN$AG$leng[Path]) - OCN$AG$leng[j]  
#         }
#         counter_down <- counter_down + 1
#         
#         # Ups: contains all nodes upstream of j, barring those upstream of the penultimate path node and  j 
#         # -> all the nodes for which j is the intersection node
#         if(iDP){
#           Ups <- setdiff(OCN$AG$upstream[[j]], c(OCN$AG$upstream[[Path[length(Path) - 1]]],j) ) 
#           
#           for (u in Ups){
#             if (includeDownstreamNode){
#               AG_DwnstrLength_unconnected[i, u] <- sum(OCN$AG$leng[Path])  # count the intersection node
#             } else {
#               AG_DwnstrLength_unconnected[i, u] <- sum(OCN$AG$leng[Path]) - OCN$AG$leng[j] # don't count the intersection node
#             }
#           }
#         }
#       }
#       if (displayUpdates){
#         if ((i %% max(1,round(OCN$AG$nNodes*0.01)))==0){
#           message(sprintf("AG downstream paths... %.1f%%\r",i/(1.001*OCN$AG$nNodes)*100), appendLF = FALSE)}}
#     }
#     indices_down <- indices_down[1:(counter_down-1), ]
#     values_down <- values_down[1:(counter_down-1)]
#     AG_DownstreamPathLength[indices_down] <- values_down
# 
#     if (displayUpdates){
#       message("AG downstream paths... 100.0%\n", appendLF = FALSE)}
#   }
#   
#   ## copy variables into OCN
#   if ("RN" %in% level){
#     if (includePaths){ OCN$RN$downstreamPath <- RN_DownstreamPath}
#     OCN$RN$downstreamPathLength <- RN_DownstreamPathLength
#     if(iDP){OCN$RN$downstreamLengthUnconnected <- RN_DwnstrLength_unconnected}
#   }
#   if ("AG" %in% level){
#     if (includePaths){ OCN$AG$downstreamPath <- AG_DownstreamPath}
#     OCN$AG$downstreamPathLength <- AG_DownstreamPathLength
#     if(iDP){OCN$AG$downstreamLengthUnconnected <- AG_DwnstrLength_unconnected}
#   }
#   
#   invisible(OCN)
#   
# }
