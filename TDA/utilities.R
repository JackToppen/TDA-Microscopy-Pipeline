
#Some functions here are from Peter Bubenik's labs that can be found at 
#https://people.clas.ufl.edu/peterbubenik/intro-to-tda/

#Also, some functions are adapted from 
#https://github.com/althomas/tda-for-worm-behavior


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

############tda functions for general_pipeline.Rmd#######################################

# remove dimension zero points from a persistence diagram
remove_dimzero_from_diagrams <- function(PD){
  zeroes <- which(PD[, 1] == 0)
  new_PD <- PD[-zeroes,] 
  return(new_PD)
}

#compute and save a list of persistence diagrams
diagrams_list <- function(data_files, cell_types, csv_files_path, save_filename, save_file_location){
  PDs <- list()
  max_birth <- list() 
  max_death <- list() 
  for (i in 1:length(data_files)){
    print(sprintf("Processing file %s", data_files[i]))
    path <- concat_path(csv_files_path, data_files[i])
    cells <- read.csv(path)
    cells <- cells[which(is.element(cells[,8], cell_types)),]
    #compute persistence homology using Delaunay complex filtration (also known as Alpha complex filtration)
    PH <-  alphaComplexDiag(cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), location = TRUE)
    PDs[[i]] <- PH[["diagram"]]
    PDs[[i]] <- remove_dimzero_from_diagrams(PDs[[i]])
    PDs[[i]][,2:3] <- sqrt(PDs[[i]][,2:3])
    #max birth of 1-degree persistence diagrams 
    max_birth[[i]] <- max(PDs[[i]][,2])
    #max death of 1-degree persistence diagrams 
    max_death[[i]] <- max(PDs[[i]][,3])
  }
  
  max_birth <- max(unlist(max_birth))
  max_death <- max(unlist(max_death))
  
  #save the list of persistence diagrams, the max birth and max death values
  save(PDs, max_birth, max_death, file=paste0(save_file_location,save_filename))
  
  return(PDs)
}

#plot persistence diagrams from a list
plot_diagrams_from_list <- function(PD_list, number_of_files, birth, death){
  par(pty="s")
  max_radius <- max(birth, death)
  for (i in 1:number_of_files){
    plot(PD_list[[i]][,2:3], asp=1, xlim=c(0, max_radius) , ylim=c(0, max_radius), xlab='', ylab='', col='orange', bty="o", pch=20, cex=1)
    abline(0,1)
  }
}

#plot representative cycles that persist (live) over a certain threshold
plot_representative_cycles <- function(data_files, cell_types, csv_files_path, threshold){
  par(pty="s")
  for (i in 1:length(data_files)){
    print(sprintf("Processing file %s", data_files[i]))
    path <- concat_path(csv_files_path, data_files[i])
    cells <- read.csv(path)
    cells <- cells[which(is.element(cells[,8], cell_types)),]
    #compute persistence homology using Delaunay complex filtration (also known as Alpha complex filtration)
    filtration <- alphaComplexFiltration(cells[,1:2], printProgress = TRUE)
    PH <-  alphaComplexDiag(cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), location = TRUE)
    PD <- PH[["diagram"]]
    #plot cycles that persist over specific persistence_param
    ones <- which(PD[,1] == 1)
    if (length(ones) > 1 ){
      plot(filtration[["coordinates"]], pch = 19, cex=0.1, ylab="", xlab="", axes = FALSE)
      for (m in ones[1]:(length(ones)+ones[1]-1)){
        cycles <- PH[["cycleLocation"]][m]
        if ((sqrt(PD[m,3]) - sqrt(PD[m,2])  > threshold) ){
          for (s in 1:length(cycles)){
            for (l in 1:dim(cycles[[s]])[1]){
              lines(cycles[[s]][l,,], col="blue", lwd=1)
            }
          } 
        } 
      }
    }
  }
}

#compute and save a list of persistence landscapes
landscapes_list <- function(PD_list, number_of_files, birth, death, discr_step, save_filename, save_file_location){
  max_x <- (death+birth)/2 + 5
  radius_values <- seq(0, max_x, discr_step)
  PLs <- list()
  for (i in 1:number_of_files){ 
    PLs[[i]] <- t(landscape(PD_list[[i]],dimension=1,KK=1:1000,radius_values))
  }
  
  #save the list of persistence diagrams, the max birth and max death values
  save(PDs, radius_values, file=paste0(save_file_location,save_filename))
  
  return(PLs)
}

#plot a list of persistence landscapes
plot_landscapes_from_list <- function(PL_list,number_of_files, birth, death, discr_step){
  par(pty="s")
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  max_x <- (death+birth)/2 + 5
  max_height <- death/2 
  radius_values <- seq(0, max_x, discr_step)
  for (i in 1:number_of_files){ 
    plot(radius_values, PL_list[[i]][1,], type="l", ylab="", xlab="", xlim=c(0, max_x) , ylim=c(0, max_height), ann=FALSE, bty="o",col=mycolors[1])
    for(k in 2:dim(PL_list[[i]])[1]){
      lines(radius_values,PL_list[[i]][k,],type="l",col=mycolors[k %% 15])
    }
  }
}

# convert a vector to a persistence landscape
landscape_from_vector <- function(PL_vector, radius_values){
  m <- length(radius_values)
  K <- length(PL_vector)/m
  PL <- Matrix(0, nrow = K, ncol=m, sparse = TRUE)
  for (i in 1:K){
    PL[i,1:m] <- PL_vector[(1+(i-1)*m):(i*m)]
  }
  return(PL)
}

# Matrix of persistence landscape row vectors from list of persistence landscapes
landscape_matrix_from_list <- function(PL_list){
  n <- length(PL_list)
  m <- ncol(PL_list[[1]])
  max_depth <- integer(n)
  for (i in 1:n)
    max_depth[i] <- nrow(PL_list[[i]])
  K <- max(max_depth)
  PL_matrix <- Matrix(0, nrow = n, ncol = m*K, sparse = TRUE)
  for (i in 1:n)
    for (j in 1:max_depth[i])
      PL_matrix[i,(1+(j-1)*m):(j*m)] <- PL_list[[i]][j,]
  return(PL_matrix)
}

#compute and save the average persistence landscape
average_persistence_landscape <- function(PL_list, birth, death, save_filename, save_file_location){
  PL_matrix <- landscape_matrix_from_list(PL_list)
  average_PL_vector <- colMeans(PL_matrix, sparseResult = TRUE)
  average_PL <- landscape_from_vector(average_PL_vector, radius_values)
  
  save(average_PL, birth, death, file=paste0(save_file_location,save_filename))
  
  return(average_PL)
}

#plot average persistence landscape
plot_average_persistence_landscape <- function(average_PL, birth, death, discr_step){
  max_x <- (death+birth)/2 + 5
  radius_values <- seq(0, max_x, discr_step)
  
  
  par(pty="s")
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  max_height <- death/2 
  plot(radius_values, average_PL[1,], type="l", ylab="", xlab="", xlim=c(0, max_x) , ylim=c(0, max_height), ann=FALSE, bty="o",col=mycolors[1])
  for(k in 2:dim(average_PL)[1]){
    lines(radius_values,average_PL[k,],type="l",col=mycolors[k %% 15])
  }
}

#compute the difference between two average persistence landscapes 
plot_difference_average_PL <- function(average_PL1, average_PL2, birth1, death1, birth2, death2, discr_step){
  if (ncol(average_PL1) < ncol(average_PL2)){
    average_PL1 <- cbind(average_PL1, matrix(0, nrow(average_PL1), ncol(average_PL2)-ncol(average_PL1)))
    difference_matrix <- average_PL1 - average_PL2
    max_x <- (max_birth2+max_death2)/2 + 5
    max_height <- death2/2 
  } else{
    average_PL2 <- cbind(average_PL2, matrix(0, nrow(average_PL2), ncol(average_PL2)-ncol(average_PL2)))
    difference_matrix <- average_PL1 - average_PL2
    max_x <- (max_birth1+max_death1)/2 + 5
    max_height <- death1/2 
  }
  
  radius_values <- seq(0, max_x, discr_step)
  
  par(pty="s")
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  
  plot(radius_values, difference_matrix[2,], type="l", ylab="", xlab="", xlim=c(0, max_x), ylim=c(-max_height, max_height), ann=FALSE, bty="o",col=mycolors[1])
  for(k in 2:dim(difference_matrix)[1]){
    lines(radius_values, difference_matrix[k,],type="l",col=mycolors[k %% 15])
  }
}

#permutation test on two groups of PLs saved as lists
permutation_test_for_PLs <- function(PL1, PL2, nrepeats = 10000){
  PL1_matrix <- landscape_matrix_from_list(PL1)
  PL2_matrix <- landscape_matrix_from_list(PL2)
  # append zeros if necessary so that the matrices have the same number of columns
  num_columns <- max(ncol(PL1_matrix),ncol(PL2_matrix))
  PL1_matrix <- cbind(PL1_matrix, Matrix(0,nrow=nrow(PL1_matrix),ncol=num_columns-ncol(PL1_matrix)))
  PL2_matrix <- cbind(PL2_matrix, Matrix(0,nrow=nrow(PL2_matrix),ncol=num_columns-ncol(PL2_matrix)))
  
  M <- rbind(PL1_matrix, PL2_matrix)
  n <- nrow(M) # number of total data points
  k <- nrow(PL1_matrix) # number of data points in first group
  observation <- euclidean_distance(colMeans(PL1_matrix),colMeans(PL2_matrix))
  ngreater_distances <- 0
  for (i in 1:nrepeats){
    permutation <- sample(1:n)
    distance <- euclidean_distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
    if (distance >= observation)
      ngreater_distances <- ngreater_distances + 1
  }
  return(print(sprintf("p-value %f:", ngreater_distances/nrepeats)))
}

############tda functions for hipsc_pipeline.Rmd ('tda-tools' package needed)###########

#plot persistence diagram
plot_diagram <- function(pairs, dgm_max){
  par(pty="s")
  finite_points <- matrix(pairs[pairs[,2] != Inf], ncol=2)
  #line below is needed if we compute homology in degree 0 
  #because we might have path connected components that die at infinity
  infinite_points <- matrix(pairs[pairs[,2] == Inf], ncol=2) 
  
  if (missing(dgm_max)){
    dgm_max <- max(pairs[pairs[,] != Inf])
    dgm_max <- dgm_max + 0.05*(dgm_max)
  }
  
  if (nrow(finite_points) > 0){
    plot(finite_points, col='orange', bty="o", asp=1, xlab='', ylab='', xlim=c(0,dgm_max), ylim=c(0,dgm_max), pch=20, cex=1)
  } else {
    plot(c(dgm_min,dgm_min), col='white', bty="o", asp=1, xlab='', ylab='', xlim=c(0,dgm_max), ylim=c(0,dgm_max), pch=20, cex=1)
  }
  points(infinite_points[,1],rep(dgm_max,nrow(infinite_points)), col='red', pch=20, cex=1)
  abline(0,1)
}

#compute persistence landscape from a persistence diagram
# returns zero persistence landscape if persistence diagram is empty (instead of giving an error)
landscape0 <- function(data, degree, exact=FALSE, dx, min_x, max_x){
  if (length(data)==0) { # empty persistence diagram
    tdatools::landscape(matrix(0, nrow=1,ncol=2), degree=degree, exact=exact, dx=dx, min_x=min_x, max_x=max_x)
  } else {
    tdatools::landscape(data, degree=degree, exact=exact, dx=dx, min_x=min_x, max_x=max_x)
  }
}

#plot persistence landscape using color scheme
plot_landscape <- function(landscape, x_max, y_min, y_max){
  par(pty="s")
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
  #plot persistence landscape with specific limits of x-axis and y-axis
  else { 
    plot(depth_1[,1],depth_1[,2], xlim=c(0, x_max) , ylim=c(y_min, y_max), type='l', ann=FALSE, bty="o",col=mycolors[1], lwd=line_width)
  }
  
  #plot persistence landscape at other depths
  if (numLevels(internal) >1 ){
    for(depth in 2:numLevels(internal)){
      depth_k <- accessLevel(internal, depth)
      lines(depth_k[,1], depth_k[,2], col=mycolors[depth %% 15], lwd=line_width)
    }
  } 
}

#convert persistence landscape to a single vector
vectorize_landscapes <- function(PL_list, depth_cap=0){
  # input: list of persistence landscapes and a highest depth of persistence landscape to include in vector
  # output: a matrix where each row is the concatenation 
  # of the y values of a persistence landscape at each depth (up to depth_cap). 
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

# convert persistence diagrams from data structure GUDHI to tdatools
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
euclidean_distance <- function(u, v) sqrt(sum((u - v) ^ 2))

#permutation test on two groups of data that are saved as matrices
permutation_test <- function(group1 , group2, num.repeats = 10000){
  # append zeros if necessary so that the matrices have the same number of columns
  num.columns <- max(ncol(group1),ncol(group2))
  group1 <- cbind(group1, Matrix(0,nrow=nrow(group1),ncol=num.columns-ncol(group1)))
  group2 <- cbind(group2, Matrix(0,nrow=nrow(group2),ncol=num.columns-ncol(group2)))
  t.obs <- euclidean.distance(colMeans(group1),colMeans(group2))
  k <- dim(group1)[1]
  M <- rbind(group1,group2)
  n <- dim(M)[1]
  count <- 0
  for (i in 1:num.repeats){
    permutation <- sample(1:n)
    t <- euclidean.distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
    if (t >= t.obs)
      count <- count + 1
  }
  return(count/num.repeats)
}

