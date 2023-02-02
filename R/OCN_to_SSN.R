OCN_to_SSN <- function (OCN, level, obsDesign, obsSites, predDesign, predSites,
                         path, randomAllocation = FALSE, importToR = FALSE) 
{

  # check input
  if (missing(obsDesign) & missing(obsSites) ) 
    stop("At least one between obsDesign and obsSites must be specified")
  if (!missing(obsDesign) & !missing(obsSites) ) 
    stop("Only one between obsDesign and obsSites must be specified")
  if (missing(predDesign) & missing(predSites)){predDesign <- noPoints}
  if (!missing(predDesign) & !missing(predSites))
    stop("predDesign and predSites cannot both be specified")
  if (missing(path)) 
    stop("Path cannot be missing")
  if (missing(OCN)) {
    stop("Input OCN cannot be missing")
  }
  if (!(level %in% c("FD","RN","AG"))) {
    stop("Invalid level")
  }
  if (length(path) != 1) 
    stop("Please enter a single path")
  info <- file.info(path)
  isdir <- info$isdir
  if (is.na(isdir)) {
    dir.create(path)
  } else if (isdir == FALSE) {
    stop("Unable to create directory")
  }
  if (!("RN" %in% names(OCN)) && (level %in% c("RN","AG"))){
    stop('Missing aggregation level in OCN. Run landscape_OCN and/or aggregate_OCN prior to OCN_to_SSN.')
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd)) # when function exits, reset original wd
  setwd(path)
  
  sub_OCN <- NULL
  eval(parse(text=(paste("sub_OCN <- OCN$",level,sep=""))))
  
  # recalculate coordinates if periodicBoundaries == TRUE
  if (OCN$periodicBoundaries==TRUE){
    X <- vector(mode="numeric", length=sub_OCN$nNodes)
    Y <- vector(mode="numeric", length=sub_OCN$nNodes)
    for (i in 1:sub_OCN$nNodes){
      node <- which(OCN$FD$X==sub_OCN$X[i] & OCN$FD$Y==sub_OCN$Y[i])
      X[i] <- OCN$FD$XDraw[node]
      Y[i] <- OCN$FD$YDraw[node]
    }
  } else {
    X <- sub_OCN$X
    Y <- sub_OCN$Y
  }
  
  #initialization
  n_networks <- length(unique(sub_OCN$toCM)) #length(sub_OCN$outlet)
  edges <- vector(mode = "list", length = n_networks)
  tree.graphs <- edges
  locations <- edges
  rids <- edges
  initial_points <- n_edges <- vector(mode = "numeric", 
                                      length = n_networks)
  edge_lengths <- edges
  cumulative_nedges <- 0
  max_x <- vector(mode = "numeric", length = n_networks)
  min_x <- vector(mode = "numeric", length = n_networks)

  list_OCN <- vector(mode = "list", length = n_networks)
  for (i in 1:n_networks){
    indices_sub <- which(sub_OCN$toCM==i)
    list_OCN[[i]]$W <- sub_OCN$W[indices_sub,indices_sub,drop=FALSE]
    list_OCN[[i]]$X <- X[indices_sub] 
    list_OCN[[i]]$Y <- Y[indices_sub] 
    list_OCN[[i]]$leng <- sub_OCN$leng[indices_sub] 
    list_OCN[[i]]$outlet <- which(rowSums(list_OCN[[i]]$W)==0)
    #list_OCN[[i]]$A <- sub_OCN$A[indices_sub] 
    list_OCN[[i]]$nNodes <- length(indices_sub) 
  }
  
  treeFunction <- function (OCN) # rewritten   
  {
    mm <- as.dgCMatrix.spam(t(OCN$W)) # transpose because trees are seen from the outlet towards the branches
    # add fake node at the outlet
    mm <- cbind(rbind(mm,rep(0,OCN$nNodes)),rep(0,OCN$nNodes + 1) )
    mm[OCN$nNodes + 1,OCN$outlet] <- 1
    
    g <- graph_from_adjacency_matrix(mm)

    locations_this_network <- matrix(c(c(OCN$X,OCN$X[OCN$outlet]),c(OCN$Y,OCN$Y[OCN$outlet])),ncol=2,nrow=OCN$nNodes+1)
    return( list(locations = locations_this_network, graph = g, initialPoint = OCN$nNodes + 1) )
  }
  
  # generate network
  for (i in 1:n_networks) {
    graph <- treeFunction(list_OCN[[i]])
    edges_this_network <- get.edgelist(graph$graph, names = FALSE)
    reordering <- order(edges_this_network[, 1])
    edges_this_network <- edges_this_network[reordering, 
                                             ]
    locations_this_network <- graph$locations[reordering]
    tree.graphs[[i]] <- graph.edgelist(edges_this_network)
    initial_points[i] <- graph$initialPoint
    locations_this_network <- graph$locations
    locations[[i]] <- locations_this_network
    edges[[i]] <- edges_this_network
    n_edges[i] = nrow(edges_this_network)
    
     tmp <- vector(mode = "numeric", length = n_edges[i])
     for (k in 1:n_edges[i]) {tmp[k] <- list_OCN[[i]]$leng[edges[[i]][k,2]]} # "reverse" length (as the OCN is reversed)
     edge_lengths[[i]] <- tmp
    
    rids[[i]] <- (1:n_edges[i]) + cumulative_nedges
    names(edge_lengths[[i]]) <- rids[[i]]
    cumulative_nedges <- cumulative_nedges + n_edges[i]
    min_x[i] <- min(locations_this_network)
    max_x[i] <- max(locations_this_network)
  }
  # cumulative_x <- max_x[1]
  # if (n_networks > 1) {
  #   for (i in 2:n_networks) {
  #     locations[[i]][, 1] <- locations[[i]][, 1] + cumulative_x - 
  #       min_x[i] + 0.1
  #     cumulative_x <- cumulative_x + max_x[i] - min_x[i] + 
  #       0.1
  #   }
  # }
  spatial_edges <- vector(mode = "list", length = sum(n_edges))
  cumulative_nedges <- 0
  for (netid in 1:n_networks) {
    locations_this_network <- locations[[netid]]
    edges_this_network <- edges[[netid]]
    for (edge.index in 1:n_edges[netid]) {
      edge <- edges_this_network[edge.index, ]
      first.location = locations_this_network[edge[1],   ]
      second.location = locations_this_network[edge[2],    ]
      spatial_edges[[edge.index + cumulative_nedges]] = Lines(list(Line(rbind(second.location, first.location))), ID=as.character(edge.index+cumulative_nedges))
    }
    cumulative_nedges = cumulative_nedges + n_edges[netid]
  }
  sl <- SpatialLines(spatial_edges)
  edge_updist <- vector(mode = "list", length = n_networks)
  line_data <- data.frame()
  shreve <- vector(mode = "list", length = n_networks)
  for (netid in 1:n_networks) {
    rids_this_network <- rids[[netid]]
    edges_this_network <- edges[[netid]]
    edge_updist_this_network = 0
    known_points = initial_points[netid]
    remaining_edges <- edges_this_network
    edge_lengths_this_network <- edge_lengths[[netid]]
    remaining_edge_lengths <- edge_lengths_this_network
    known_rids <- c()
    remaining_rids <- rids_this_network
    while (TRUE) {
      can_calculate <- (remaining_edges[, 1] %in% known_points) & 
        (!(remaining_edges[, 2] %in% known_points)) # edges departing from knowing nodes and not entering known nodes
      upstream_point_indicies <- remaining_edges[can_calculate, 
                                                 2] # find nodes upstream and downstream of these nodes
      downstream_point_indicies <- remaining_edges[can_calculate, 
                                                   1]
      edge_updist_this_network <- c(edge_updist_this_network, 
                                    edge_updist_this_network[match(downstream_point_indicies, 
                                                                   known_points)] + remaining_edge_lengths[can_calculate]) 
      known_points <- c(known_points, upstream_point_indicies) # all points upstream of known points are now known
      remaining_edges <- remaining_edges[!can_calculate, 
                                         , drop = FALSE] # remove edges that have just been assessed
      remaining_edge_lengths <- remaining_edge_lengths[!can_calculate]
      known_rids <- c(known_rids, remaining_rids[can_calculate])
      remaining_rids <- remaining_rids[!can_calculate]
      if (length(remaining_edges) == 0) 
        break
    }
    edge_updist_this_network <- edge_updist_this_network[-1][order(known_rids)]
    names(edge_updist_this_network) <- sort(known_rids)
    edge_updist[[netid]] <- edge_updist_this_network
    shreve_this_network <- vector(mode = "numeric", 
                                  length = n_edges[netid])
    is_initial <- !(edges_this_network[, 2] %in% edges_this_network[, 
                                                                    1])
    known_points <- edges_this_network[which(is_initial), 
                                       2]
    remaining_points <- edges_this_network[which(!is_initial), 
                                           2]
    shreve_values_points <- rep(1, length(known_points))
    shreve_values_rids <- c()
    known_rids <- c()
    remaining_rids <- rids_this_network - min(rids_this_network) + 
      1
    remaining_edges <- edges_this_network
    while (TRUE) {
      can_calculate <- (remaining_edges[, 2] %in% known_points) # now evaluate edges entering known nodes
      shreve_values_rids <- c(shreve_values_rids, shreve_values_points[match(remaining_edges[can_calculate, 
                                                                                             2], known_points)])
      remaining_edges <- remaining_edges[!can_calculate, 
                                         , drop = FALSE]
      known_rids <- c(known_rids, remaining_rids[can_calculate])
      remaining_rids <- remaining_rids[!can_calculate]
      if (length(remaining_edges) == 0) 
        break
      can_calculate <- !(remaining_points %in% remaining_edges[, 
                                                               1])
      can_calculate_function <- function(index) {
        relevant_edges <- which(edges_this_network[, 
                                                   1] == remaining_points[index])
        return(all(relevant_edges %in% known_rids))
      }
      can_calculate[can_calculate] <- sapply(which(can_calculate), 
                                             can_calculate_function)
      calculate_shreve <- function(index) {
        relevant_edges <- which(edges_this_network[, 
                                                   1] == remaining_points[index])
        relevant_known_rids <- match(relevant_edges, 
                                     known_rids)
        return(sum(shreve_values_rids[relevant_known_rids]))
      }
      new_shreve_values_points <- sapply(which(can_calculate), 
                                         calculate_shreve)
      shreve_values_points <- c(shreve_values_points, new_shreve_values_points)
      known_points <- c(known_points, remaining_points[can_calculate])
      remaining_points <- remaining_points[!can_calculate]
    }
    known_rids_ordering <- order(known_rids)
    shreve_values_rids <- shreve_values_rids[known_rids_ordering]
    additive_function_values <- shreve_values_rids/max(shreve_values_rids)
    additional_line_data <- data.frame(rid = rids_this_network, 
                                       netID = as.factor(rep(netid, n_edges[netid])), upDist = edge_updist_this_network, 
                                       shreve = shreve_values_rids, Length = edge_lengths_this_network, 
                                       addfunccol = additive_function_values)
    rownames(additional_line_data) <- rids_this_network
    line_data <- rbind(additional_line_data, line_data)
  }
  sldf <- SpatialLinesDataFrame(sl, data = line_data, match.ID = TRUE)
  writeOGR(sldf, ".", "edges", verbose = FALSE,  driver = "ESRI Shapefile")
  setwd(old_wd)
  binary_ids_tables <- list()
  for (netid in 1:n_networks) {
    remaining_edges <- edges[[netid]]
    remaining_points <- unique(remaining_edges)
    known_point_indicies <- c(initial_points[netid])
    known_binaryids <- c("")
    known_rids <- c()
    remaining_rids <- rids[[netid]]
    while (TRUE) {
      next_points <- (remaining_edges[, 2] %in% remaining_points) & (remaining_edges[, 1] %in% known_point_indicies) 
      upstream_points <- remaining_edges[next_points, 2]
      downstream_points <- remaining_edges[next_points, 1]
      additional_binary_bits <- c()
      counter <- 1
      while (counter <= length(downstream_points)) { # need to change with an index in base 8
        if (counter == length(downstream_points) || downstream_points[counter] !=  downstream_points[counter + 1]) {
          additional_binary_bits <- c(additional_binary_bits,  "1")
          counter <- counter + 1
        } else if (counter == (length(downstream_points) - 1) || downstream_points[counter] !=  downstream_points[counter + 2] ) {
          additional_binary_bits <- c(additional_binary_bits,  "1", "2")
          counter <- counter + 2
        } else if (counter == (length(downstream_points) - 2) || downstream_points[counter] !=  downstream_points[counter + 3]) {
          additional_binary_bits <- c(additional_binary_bits,  "1", "2", "3")
          counter <- counter + 3
        } else if (counter == (length(downstream_points) - 3) || downstream_points[counter] !=  downstream_points[counter + 4]) {
          additional_binary_bits <- c(additional_binary_bits,  "1", "2", "3","4")
          counter <- counter + 4
        } else if (counter == (length(downstream_points) - 4) || downstream_points[counter] !=  downstream_points[counter + 5]) {
          additional_binary_bits <- c(additional_binary_bits,  "1", "2", "3","4","5")
          counter <- counter + 5
        } else if (counter == (length(downstream_points) - 5) || downstream_points[counter] !=  downstream_points[counter + 6]) {
          additional_binary_bits <- c(additional_binary_bits,  "1", "2", "3","4","5","6")
          counter <- counter + 6  
        } else if (counter == (length(downstream_points) - 6) || downstream_points[counter] !=  downstream_points[counter + 7]) {
          additional_binary_bits <- c(additional_binary_bits,  "1", "2", "3","4","5","6","7")
          counter <- counter + 7  
        } else if (counter == (length(downstream_points) - 7) || downstream_points[counter] !=  downstream_points[counter + 8]) {
          additional_binary_bits <- c(additional_binary_bits,  "1", "2", "3","4","5","6","7","8")
          counter <- counter + 8}  
      }
      binaryid_indicies <- match(downstream_points, known_point_indicies)
      previous_binary_ids <- known_binaryids[binaryid_indicies]
      known_binaryids <- c(known_binaryids, paste(previous_binary_ids, additional_binary_bits, sep = ""))
      known_point_indicies <- c(known_point_indicies, upstream_points)
      known_rids <- c(known_rids, remaining_rids[next_points])
      remaining_edges <- remaining_edges[!next_points, , drop = FALSE]
      remaining_rids <- remaining_rids[!next_points]
      remaining_points <- remaining_points[remaining_points %in%  remaining_edges[, 2]]
      if (length(remaining_rids) == 0) 
        break
    }
    known_binaryids <- known_binaryids[-1]
    binary_ids <- known_binaryids[order(known_rids)]
    binary_ids_tables[[netid]] <- data.frame(rid = rids[[netid]],  binaryID = binary_ids)
    if (length(unique(binary_ids_tables[[netid]]$binaryID)) !=  length(binary_ids_tables[[netid]]$binaryID)) 
      stop("Internal error")
    write.table(binary_ids_tables[[netid]], file = file.path(path, paste("netID", netid, ".dat", sep = "")),  col.names = T, sep = ",", row.names = FALSE)
  }
  setwd(path)
  distance_matrices <- list()
  for (netid in 1:n_networks) {
    edge_updist_this_network <- edge_updist[[netid]]
    edges_this_network <- edges[[netid]]
    binary_id_table <- binary_ids_tables[[netid]]
    distance_matrix <- matrix(0, nrow(binary_id_table), nrow(binary_id_table))
    colnames(distance_matrix) <- rownames(distance_matrix) <- binary_id_table$rid
    partial_match_function <- function(binary_id1, binary_id2) {
      min_len <- min(nchar(binary_id1), nchar(binary_id2))
      for (j in 1:min_len) {
        if (substr(binary_id1, j, j) != substr(binary_id2,j, j)) 
          return(j - 1)
      }
      return(min_len)
    } # number of figures in the binary id (starting from left) that are identical (works also for non-binary strings)
    character_binary_ids <- as.character(binary_id_table$binaryID)
    for (i in 1:nrow(binary_id_table)) {
      current_binary_id <- binary_id_table$binaryID[i]
      current_rid <- as.character(binary_id_table$rid[i])
      current_updist <- edge_updist_this_network[current_rid]
      matching_characters <- sapply(character_binary_ids, partial_match_function, as.character(current_binary_id))
      matching_substring <- substr(binary_id_table$binaryID, 1, matching_characters)
      indices <- match(matching_substring, binary_id_table$binaryID)
      if (any(is.na(indices))) 
        stop("Internal Error")
      downstream_rids <- binary_id_table$rid[indices]
      downstream_updists <- edge_updist_this_network[as.character(downstream_rids)]
      distance_matrix[as.character(current_rid), ] <- pmax(current_updist -  downstream_updists, rep(0, length(downstream_updists)))
    }
    reindex <- match(binary_id_table$rid, rids[[netid]])
    distance_matrix <- distance_matrix[reindex, reindex]
    distance_matrices[[netid]] <- distance_matrix + t(distance_matrix)
  }
  
  # observation sites
  if (missing(obsSites)){
  obs_sites <- obsDesign(tree.graphs, edge_lengths, locations,  edge_updist, distance_matrices)
  } else { obs_sites <- list(n_networks); k <- 0
    for (i in 1:n_networks){
      nodes <- obsSites[which(sub_OCN$toCM[obsSites]==i)]
      # identify corresponding nodes in igraph
      for (j in 1:length(nodes)){nodes[j] <- which(locations[[i]][,1]==sub_OCN$X[nodes[j]] & locations[[i]][,2]==sub_OCN$Y[nodes[j]])[1]}
      edge_vec  <- numeric(0)
      for (s in nodes){
        edge_vec <- c(edge_vec, rids[[i]][which(edges[[i]][,2]==s)])
        k <- k+1}
      if (randomAllocation){ratio_vec <- runif(length(edge_vec))} else {ratio_vec <- 1+numeric(length(edge_vec))}
      obs_sites[[i]] <- data.frame(edge=edge_vec, ratio=ratio_vec, locID=((k-length(edge_vec)+1):k))
    }}
  if (missing(predSites)){
    pred_sites <- predDesign(tree.graphs, edge_lengths, locations,  edge_updist, distance_matrices)
  } else { pred_sites <- list(n_networks); k <- 0
  for (i in 1:n_networks){
    nodes <- predSites[which(sub_OCN$toCM[predSites]==i)]
    # identify corresponding nodes in igraph
    for (j in 1:length(nodes)){nodes[j] <- which(locations[[i]][,1]==sub_OCN$X[nodes[j]] & locations[[i]][,2]==sub_OCN$Y[nodes[j]])[1]}
    edge_vec  <- numeric(0)
    for (s in nodes){
      edge_vec <- c(edge_vec, rids[[i]][which(edges[[i]][,2]==s)])
      k <- k+1}
    if (randomAllocation){ratio_vec <- runif(length(edge_vec))} else {ratio_vec <- 1+numeric(length(edge_vec))}
    pred_sites[[i]] <- data.frame(edge=edge_vec, ratio=ratio_vec, locID=((k-length(edge_vec)+1):k))
  }}
  max_observed_locID <- max(unlist(lapply(obs_sites, function(x) max(x$locID))))
  for (i in 1:length(pred_sites)) {
    pred_sites[[i]]$locID <- pred_sites[[i]]$locID + max_observed_locID
  }
  n_obs_sites <- unlist(lapply(obs_sites, function(x) return(dim(x)[1])))
  n_pred_sites <- unlist(lapply(pred_sites, function(x) return(dim(x)[1])))
  sites_data <- data.frame()
  combined_site_location_data <- c()
  pred_data <- data.frame()
  combined_pred_location_data <- c()
  cumulative_pids <- 0
  
  for (netid in 1:n_networks) {
    edge_lengths_this_network <- edge_lengths[[netid]]
    edges_this_network <- edges[[netid]]
    edge_updist_this_network <- edge_updist[[netid]]
    locations_this_network <- locations[[netid]]
    n_locations_this_network <- n_obs_sites[netid] + n_pred_sites[netid]
    rids_this_network <- rids[[netid]]
    pred_sites_this_network <- pred_sites[[netid]]
    obs_sites_this_network <- obs_sites[[netid]]
    f <- function(row) {
      rid <- as.character(row[1])
      edge_id <- match(rid, rids_this_network)
      proportion <- as.numeric(row[2])
      downstream_point <- edges_this_network[edge_id, 1]
      upstream_point <- edges_this_network[edge_id, 2]
      downstream_location <- locations_this_network[downstream_point, 
                                                    ]
      upstream_location <- locations_this_network[upstream_point, 
                                                  ]
      location <- downstream_location + proportion * (upstream_location - 
                                                        downstream_location)
      ret <- c(location, (sqrt(sum((location - downstream_location)^2)) + 
                            edge_updist_this_network[rid] - edge_lengths_this_network[rid]))
      names(ret) <- c("NEAR_X", "NEAR_Y", "upDist")
      return(ret)
    }
    obs_location_data <- data.frame(rid = obs_sites_this_network$edge, 
                                    ratio = obs_sites_this_network$ratio, locID = obs_sites_this_network$locID, 
                                    stringsAsFactors = FALSE)
    pred_location_data <- data.frame(rid = pred_sites_this_network$edge, 
                                     ratio = pred_sites_this_network$ratio, locID = pred_sites_this_network$locID, 
                                     stringsAsFactors = FALSE)
    if (n_locations_this_network > 0) {
      pred_location_data_this_network <- t(apply(pred_location_data, 
                                                 1, f))
      if (length(pred_location_data_this_network) > 0) {
        colnames(pred_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
      }
      else {
        pred_location_data_this_network <- matrix(0, 
                                                  0, 3)
        colnames(pred_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
      }
      obs_location_data_this_network <- t(apply(obs_location_data, 
                                                1, f))
      if (length(obs_location_data_this_network) > 0) {
        colnames(obs_location_data_this_network) <- c("NEAR_X",  "NEAR_Y", "upDist")
      }
      else {
        obs_location_data_this_network <- matrix(0, 0, 
                                                 3)
        colnames(obs_location_data_this_network) <- c("NEAR_X", "NEAR_Y", "upDist")
      }
      if (n_obs_sites[netid] > 0) {
        obs_pids <- (1:n_obs_sites[netid]) + cumulative_pids
      }
      else obs_pids <- integer(0)
      if (n_pred_sites[netid] > 0) {
        pred_pids <- n_obs_sites[netid] + (1:n_pred_sites[netid]) +  cumulative_pids
      }
      else pred_pids <- integer(0)
    }
    else {
      obs_location_data_this_network <- pred_location_data_this_network <- data.frame(NEAR_X = numeric(0), NEAR_Y = numeric(0), upDist = numeric(0))
      obs_pids <- integer(0)
      pred_pids <- integer(0)
    }
    cumulative_pids <- cumulative_pids + n_locations_this_network
    obs_data_this_network <- data.frame(locID = obs_location_data[, 
                                                                  "locID"], upDist = obs_location_data_this_network[, 
                                                                                                                    "upDist"], pid = obs_pids, netID = rep(netid, 
                                                                                                                                                           length(obs_pids)), rid = obs_location_data[, "rid"], 
                                        ratio = obs_location_data[, "ratio"], shreve = line_data[match(obs_location_data[, 
                                                                                                                         "rid"], line_data[, "rid"]), "shreve"], 
                                        addfunccol = line_data[match(obs_location_data[, 
                                                                                       "rid"], line_data[, "rid"]), "addfunccol"], 
                                        stringsAsFactors = FALSE)
    if (ncol(obs_sites_this_network) > 3) {
      obs_data_this_network <- cbind(obs_data_this_network, 
                                     obs_sites_this_network[, -match(c("edge", 
                                                                       "ratio", "locID"), colnames(obs_sites_this_network)), 
                                                            drop = FALSE])
    }
    rownames(obs_data_this_network) <- obs_pids
    rownames(obs_location_data_this_network) <- obs_pids
    pred_data_this_network <- data.frame(locID = pred_location_data[, 
                                                                    "locID"], upDist = pred_location_data_this_network[, 
                                                                                                                       "upDist"], pid = pred_pids, netID = rep(netid, 
                                                                                                                                                               length(pred_pids)), rid = pred_location_data[, "rid"], 
                                         ratio = pred_location_data[, "ratio"], shreve = line_data[match(pred_location_data[, 
                                                                                                                            "rid"], line_data[, "rid"]), "shreve"], 
                                         addfunccol = line_data[match(pred_location_data[, 
                                                                                         "rid"], line_data[, "rid"]), "addfunccol"], 
                                         stringsAsFactors = FALSE)
    if (ncol(pred_sites_this_network) > 3) {
      pred_data_this_network <- cbind(pred_data_this_network, 
                                      pred_sites_this_network[, -match(c("edge", 
                                                                         "ratio", "locID"), colnames(pred_sites_this_network)), 
                                                              drop = FALSE])
    }
    rownames(pred_data_this_network) <- pred_pids
    rownames(pred_location_data_this_network) <- pred_pids
    if (n_obs_sites[netid] > 0) {
      sites_data <- rbind(obs_data_this_network[1:n_obs_sites[netid], 
                                                , drop = FALSE], sites_data)
      combined_site_location_data <- rbind(obs_location_data_this_network[1:n_obs_sites[netid], 
                                                                          , drop = FALSE], combined_site_location_data)
    }
    if (n_pred_sites[netid] > 0) {
      pred_data <- rbind(pred_data_this_network[(1:n_pred_sites[netid]), 
                                                , drop = FALSE], pred_data)
      combined_pred_location_data <- rbind(pred_location_data_this_network[(1:n_pred_sites[netid]), 
                                                                           , drop = FALSE], combined_pred_location_data)
    }
  }
  
  if (length(combined_site_location_data) == 0) 
    stop("At least one observation site must be present")
  sites <- SpatialPointsDataFrame(combined_site_location_data[, c("NEAR_X", "NEAR_Y"), drop = FALSE], sites_data, match.ID = TRUE)
  writeOGR(sites, ".", "sites", verbose = FALSE, 
           driver = "ESRI Shapefile")
  if (length(combined_pred_location_data) > 0) {
    preds <- SpatialPointsDataFrame(combined_pred_location_data[, 
                                                                c("NEAR_X", "NEAR_Y"), drop = FALSE], 
                                    pred_data, match.ID = TRUE)
    writeOGR(preds, ".", "preds", verbose = FALSE, 
             driver = "ESRI Shapefile")
  }
  setwd(old_wd)
  
  
  if (importToR) {
    if (sum(n_pred_sites) > 0) 
      return(importSSN(path, predpts = "preds", o.write = TRUE))
    else return(importSSN(path, o.write = TRUE))
  }
  else return(invisible(NULL))
}