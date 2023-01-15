

get_csvs <- function(path) {
    # get list of discretized CSVs recursively
    files <- list.files(path, pattern="\\.csv$", include.dirs=TRUE, recursive=TRUE)
}


concat_path <- function(path, filename) {
  # make sure path has file separator at the end
  n <- nchar(path)
  if (substr(path, n, n) != .Platform$file.sep) {
    path <- paste0(path, .Platform$file.sep)
  }

  # return correct joined path to file
  return(paste0(path, filename))
}

#plot persistence diagrams
plot_diagram <- function(pairs, dgm_max, dgm_min){
  finite_points <- matrix(pairs[pairs[,2] != Inf], ncol=2)
  #line below is needed if we compute homology in degree 0 
  #because we might have path connected components that die at infinity
  infinite_points <- matrix(pairs[pairs[,2] == Inf], ncol=2) 
  if (missing(dgm_min)){
    dgm_min <- min(pairs)
  }
  if (missing(dgm_max)){
    dgm_max <- max(pairs[pairs[,] != Inf])
    dgm_max <- dgm_max + 0.05*(dgm_max-dgm_min)
  }
  
  if (nrow(finite_points) > 0){
    plot(finite_points, col='orange', bty="o", asp=1, xlab='', ylab='', xlim=c(dgm_min,dgm_max), ylim=c(dgm_min,dgm_max), cex=1)
  } else {
    plot(c(dgm_min,dgm_min), col='white', bty="o", asp=1, xlab='', ylab='', xlim=c(dgm_min,dgm_max), ylim=c(dgm_min,dgm_max), cex=1)
  }
  points(infinite_points[,1],rep(dgm_max,nrow(infinite_points)), col='red', pch=20, cex=1)
  abline(0,1)
}

#compute PL from persistence diagram
# returns zero PL if persistence diagram is empty (instead of giving an error)
landscape0 <- function(data, degree, exact=FALSE, dx, min_x, max_x){
  if (length(data)==0) { # empty persistence diagram
    tdatools::landscape(matrix(0, nrow=1,ncol=2), degree=degree, exact=exact, dx=dx, min_x=min_x, max_x=max_x)
  } else {
    tdatools::landscape(data, degree=degree, exact=exact, dx=dx, min_x=min_x, max_x=max_x)
  }
}

#plot PL using color scheme
plot_landscape <- function(landscape, x_max, y_max){
  internal <- landscape$getInternal()
  infinity_sub <- -1
  line_width <- 1
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  depth_1 <- accessLevel(internal,1)
  
  if ((missing(x_max) | missing(y_max))){ 
    plot(depth_1[,1], depth_1[,2], type='l', xlab='', ylab='', ann=FALSE, bty="o", col=mycolors[1], lwd=line_width)
  } 
  #plot PL with specific limits of x-axis and y-axis
  else { 
    plot(depth_1[,1],depth_1[,2], xlim=c(0, x_max) , ylim=c(0, y_max), type='l', ann=FALSE, bty="o",col=mycolors[1], lwd=line_width)
  }
  
  #plot PL at other depths
  if (numLevels(internal) >1 ){
    for(depth in 2:numLevels(internal)){
      depth_k <- accessLevel(internal, depth)
      lines(depth_k[,1], depth_k[,2], col=mycolors[depth %% 15], lwd=line_width)
    }
  } 
}

#plot PL from depth d1 to depth d2 using color scheme
plot_landscape_depths <- function(landscape, d1, d2, x_max, y_max){
  internal <- landscape$getInternal()
  infinity_sub <- -1
  line_width <- 1
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  depth_1 <- accessLevel(internal,d1)
  
  #if d1 is larger than the max depth of PL then plot zero PL
  if (d1 >= numLevels(internal)){
    plot(NULL,  xlab='', ylab='', ann=FALSE, bty="o", col=mycolors[1], lwd=line_width)
  } 
  #if d1 is smaller than the max depth of PL then plot PL at depth d1
  else if ((missing(x_max) | missing(y_max)) & d1 < numLevels(internal)){ 
    plot(depth_1[,1], depth_1[,2], type='l', xlab='', ylab='', ann=FALSE, bty="o", col=mycolors[1], lwd=line_width)
  } 
  #plot PL with specific limits of x-axis and y-axis
  else if (d1 < numLevels(internal)){ 
    plot(depth_1[,1],depth_1[,2], xlim=c(0, x_max) , ylim=c(0, y_max), type='l', ann=FALSE, bty="o",col=mycolors[1], lwd=line_width)
  }
  
  #plot PL at other depths
  if (d1 < numLevels(internal) & d2 <= numLevels(internal)){
    for(depth in d1:d2){
      depth_k <- accessLevel(internal, depth)
      lines(depth_k[,1], depth_k[,2], col=mycolors[depth %% 15], lwd=line_width)
    }
  } 
  #if d2 is larger than the max depth of PL then plot PL at all depth starting at d1   
  if (d1 < numLevels(internal) & d2 >= numLevels(internal)){
    for(depth in d1:numLevels(internal)){
      depth_k <- accessLevel(internal, depth)
      lines(depth_k[,1], depth_k[,2], col=mycolors[depth %% 15], lwd=line_width)
    }
  }
}

#convert PL to a single vector
vectorize_landscapes <- function(PL_list, depth_cap=0){
  # input: list of PLs and a highest depth of PL to include in vector
  # output: a matrix where each row is the concatenation 
  # of the y values of a PL at each depth (up to depth_cap). 
  if (depth_cap == 0){ # no depth cap -- set to max depth
    max_depth <- 0
    for (i in 1:length(PL_list)){
      max_depth <- max(max_depth, dim(PL_list[[i]]$getInternal())[1])
    }
    depth_cap <- max_depth
  }
  
  vect_length <- depth_cap*dim(PL_list[[1]]$getInternal())[2]
  vect_PLs <- matrix(0, nrow=length(PL_list), ncol=vect_length)
  for (i in 1:length(PL_list)){
    if (dim(PL_list[[i]]$getInternal())[1] < depth_cap) {
      # need the transpose because as.vector takes columns of a matrix, not rows
      temp_vec <- as.vector(t(PL_list[[i]]$getInternal()[,,2]))
      vect_PLs[i,1:length(temp_vec)] <- temp_vec
    } else {
      temp_vec <- as.vector(t(PL_list[[i]]$getInternal()[1:depth_cap,,2]))
      vect_PLs[i,] <- temp_vec
    }
  }
  return(vect_PLs)
}

# convert persistence diagrams from data structure GUDHI (Dionysus) to tda-tools (Ripser)
gudhi2tdatools <- function(gudhi_dgm) {
  # extract first homology from gudhi diagram
  if (length(which(gudhi_dgm[,1]==1)) >1){
    h1 <- gudhi_dgm[which(gudhi_dgm[,1]==1),][,2:3] # size is dim(h1) x 2
  }
  if (length(which(gudhi_dgm[,1]==1))==1){
    h <- gudhi_dgm[which(gudhi_dgm[,1]==1),2:3] # size is dim(h1) x 2
    h1 <- rbind(rep(0,2),h)
  }
  if (length(which(gudhi_dgm[,1]==1))==0){
    h1 <- rbind(rep(0,2),rep(0,2))
  }
  # extract zero homology from gudhi diagram
  if (length(which(gudhi_dgm[,1]==0)) >1){
    h0 <- gudhi_dgm[which(gudhi_dgm[,1]==0),][,2:3] # size is dim(h1) x 2
  }
  if (length(which(gudhi_dgm[,1]==0))==1){
    h <- gudhi_dgm[which(gudhi_dgm[,1]==0),2:3] # size is dim(h1) x 2
    h0 <- rbind(rep(0,2),h)
  }
  if (length(which(gudhi_dgm[,1]==0))==0){
    h0 <- rbind(rep(0,2),rep(0,2))
  }
  # create tdatools diagram; load with zero and first homology
  tdatools_dgm <- diagram(0,'point-cloud', dim_max = 1)
  tdatools_dgm$pairs[[2]] <- h1
  tdatools_dgm$pairs[[1]] <- h0
  
  return(tdatools_dgm)
}

#Euclidean distance between two vectors 
euclidean.distance <- function(u, v) sqrt(sum((u - v) ^ 2))

#permutation test for two group of data that are saved as matrices
permutation_test <- function(group1, group2, nrepeats = 10000){
  M <- rbind(group1, group2)
  n <- nrow(M) # number of total data points
  k <- nrow(group1) # number of data points in first group
  observation <- euclidean.distance(colMeans(group1),colMeans(group2))
  ngreater_distances <- 0
  for (i in 1:nrepeats){
    permutation <- sample(1:n)
    distance <- euclidean.distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
    if (distance >= observation)
      ngreater_distances <- ngreater_distances + 1
  }
  return(ngreater_distances/nrepeats)
}
